/**
 * @file gvcf_compressor.cpp
 * @brief gVCF compression implementation
 */

#include "gvcf_compressor.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "../compression_strategy.h"
#include "../logger.h"
#include <algorithm>
#include <cstring>
#include <sstream>

namespace gvcf {

// Global compression strategy for gVCF module
static std::unique_ptr<CompressionStrategy> g_gvcf_strategy;
static compression_backend_t g_gvcf_backend = compression_backend_t::bsc;  // Default to bsc for best ratio/speed balance

// Set the compression backend (called before compression starts)
void SetGVCFBackend(compression_backend_t backend) {
    if (g_gvcf_backend != backend || !g_gvcf_strategy) {
        g_gvcf_backend = backend;
        g_gvcf_strategy.reset();  // Force re-creation
    }
}

compression_backend_t GetGVCFBackend() {
    return g_gvcf_backend;
}

static void EnsureGVCFStrategy() {
    if (!g_gvcf_strategy) {
        InitializeCompressionBackend(g_gvcf_backend);
        g_gvcf_strategy = MakeCompressionStrategy(g_gvcf_backend, p_bsc_fixed_fields);
    }
}

static bool GVCFCompressData(const uint8_t* data, size_t size, std::vector<uint8_t>& output) {
    EnsureGVCFStrategy();
    if (!g_gvcf_strategy) return false;
    std::vector<uint8_t> input(data, data + size);
    return g_gvcf_strategy->Compress(input, output);
}

static bool GVCFDecompressData(const uint8_t* data, size_t size, std::vector<uint8_t>& output) {
    EnsureGVCFStrategy();
    if (!g_gvcf_strategy) return false;
    std::vector<uint8_t> input(data, data + size);
    return g_gvcf_strategy->Decompress(input, output);
}

// ============================================================================
// GVCFCompressor Implementation
// ============================================================================

GVCFCompressor::GVCFCompressor(const GSC_Params& params)
    : params_(params)
    , vcf_file_(nullptr)
    , vcf_header_(nullptr)
    , vcf_record_(nullptr)
    , output_file_(nullptr)
    , current_block_id_(0)
    , num_samples_(0)
    , is_gvcf_(false)
    , has_end_field_(false)
    , has_min_dp_(false)
    , has_non_ref_(false)
{
    std::memset(&stats_, 0, sizeof(stats_));
    config_.block_size = 10000;
    config_.enable_adaptive = true;
    config_.mask_threshold = 0.7f;
}

GVCFCompressor::~GVCFCompressor() {
    CloseInput();
    if (output_file_) {
        fclose(output_file_);
        output_file_ = nullptr;
    }
}

bool GVCFCompressor::Compress() {
    auto logger = LogManager::Instance().Logger();

    // Initialize compression backend with user-specified or default backend
    SetGVCFBackend(params_.backend);
    InitializeCompressionBackend(params_.backend);
    backend_ = CreateBackend();

    logger->info("Using compression backend: {}",
        params_.backend == compression_backend_t::bsc ? "bsc" :
        params_.backend == compression_backend_t::zstd ? "zstd" : "brotli");

    // Open input file
    if (!OpenInput()) {
        logger->error("Failed to open input file: {}", params_.in_file_name);
        return false;
    }

    // Read header
    if (!ReadHeader()) {
        logger->error("Failed to read VCF header");
        return false;
    }

    // Check if it's gVCF format
    if (!DetectGVCFFormat()) {
        logger->warn("Input file does not appear to be gVCF format, proceeding anyway");
    }

    // Check sample count
    if (num_samples_ != 1) {
        logger->error("gVCF compression currently supports single-sample files only (found {} samples)",
                     num_samples_);
        return false;
    }

    logger->info("Starting gVCF compression: {} sample, {} format",
                num_samples_, is_gvcf_ ? "gVCF" : "VCF");

    // Open output file
    output_filename_ = params_.out_file_name;
    if (output_filename_.empty() || output_filename_ == "-") {
        logger->error("Output filename required for gVCF compression");
        return false;
    }

    output_file_ = fopen(output_filename_.c_str(), "wb");
    if (!output_file_) {
        logger->error("Failed to open output file: {}", output_filename_);
        return false;
    }

    // Write file header
    if (!WriteFileHeader()) {
        logger->error("Failed to write file header");
        return false;
    }

    // Process blocks
    GVCFBlock block;
    block.Reserve(config_.block_size);

    while (true) {
        block.Clear();

        if (!ReadBlock(block)) {
            break;
        }

        if (block.variant_count == 0) {
            break;
        }

        stats_.total_variants += block.variant_count;
        stats_.total_blocks++;

        if (!CompressAndWriteBlock(block)) {
            logger->error("Failed to compress block {}", current_block_id_);
            return false;
        }

        current_block_id_++;

        if (stats_.total_blocks % 100 == 0) {
            logger->debug("Processed {} blocks, {} variants",
                         stats_.total_blocks, stats_.total_variants);
        }
    }

    // Write file footer
    if (!WriteFileFooter()) {
        logger->error("Failed to write file footer");
        return false;
    }

    // Calculate final statistics
    fseek(output_file_, 0, SEEK_END);
    stats_.compressed_size = ftell(output_file_);
    if (stats_.original_size > 0) {
        stats_.compression_ratio = static_cast<float>(stats_.compressed_size) /
                                   static_cast<float>(stats_.original_size);
    }

    logger->info("gVCF compression complete: {} variants in {} blocks",
                stats_.total_variants, stats_.total_blocks);
    logger->info("Compression ratio: {:.2f}% ({} -> {} bytes)",
                stats_.compression_ratio * 100.0f,
                stats_.original_size, stats_.compressed_size);

    // Output per-field statistics
    uint64_t total_field = stats_.chrom_compressed + stats_.pos_compressed +
                           stats_.gt_compressed + stats_.dp_compressed +
                           stats_.gq_compressed + stats_.pl_compressed;
    if (total_field > 0) {
        logger->info("Field breakdown: CHROM={} POS={} GT={} DP={} GQ={} PL={}",
                    stats_.chrom_compressed, stats_.pos_compressed,
                    stats_.gt_compressed, stats_.dp_compressed,
                    stats_.gq_compressed, stats_.pl_compressed);
    }

    CloseInput();
    fclose(output_file_);
    output_file_ = nullptr;

    return true;
}

bool GVCFCompressor::OpenInput() {
    htsFile* fp = hts_open(params_.in_file_name.c_str(), "r");
    if (!fp) {
        return false;
    }
    vcf_file_ = fp;

    bcf1_t* rec = bcf_init1();
    if (!rec) {
        hts_close(fp);
        vcf_file_ = nullptr;
        return false;
    }
    vcf_record_ = rec;

    return true;
}

bool GVCFCompressor::ReadHeader() {
    htsFile* fp = static_cast<htsFile*>(vcf_file_);
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        return false;
    }
    vcf_header_ = hdr;

    // Get sample info
    num_samples_ = bcf_hdr_nsamples(hdr);
    sample_names_.clear();
    for (uint32_t i = 0; i < num_samples_; i++) {
        sample_names_.push_back(hdr->samples[i]);
    }

    // Get header text
    kstring_t str = {0, 0, nullptr};
    bcf_hdr_format(hdr, 0, &str);
    if (str.s) {
        header_text_ = std::string(str.s, str.l);
        free(str.s);
    }

    return true;
}

bool GVCFCompressor::DetectGVCFFormat() {
    bcf_hdr_t* hdr = static_cast<bcf_hdr_t*>(vcf_header_);

    // Check for INFO/END field
    has_end_field_ = (bcf_hdr_id2int(hdr, BCF_DT_ID, "END") >= 0);

    // Check for FORMAT/MIN_DP field
    has_min_dp_ = (bcf_hdr_id2int(hdr, BCF_DT_ID, "MIN_DP") >= 0);

    // Check for <NON_REF> in ALT header
    has_non_ref_ = (header_text_.find("<NON_REF>") != std::string::npos ||
                   header_text_.find("NON_REF") != std::string::npos);

    is_gvcf_ = has_end_field_ || has_min_dp_ || has_non_ref_;

    auto logger = LogManager::Instance().Logger();
    logger->debug("gVCF detection: END={}, MIN_DP={}, NON_REF={}, is_gvcf={}",
                 has_end_field_, has_min_dp_, has_non_ref_, is_gvcf_);

    return is_gvcf_;
}

bool GVCFCompressor::ReadBlock(GVCFBlock& block) {
    htsFile* fp = static_cast<htsFile*>(vcf_file_);
    bcf_hdr_t* hdr = static_cast<bcf_hdr_t*>(vcf_header_);
    bcf1_t* rec = static_cast<bcf1_t*>(vcf_record_);

    block.sample_count = 1;

    int* gt_arr = nullptr;
    int ngt_arr = 0;
    int* dp_arr = nullptr;
    int ndp_arr = 0;
    int* gq_arr = nullptr;
    int ngq_arr = 0;
    int* min_dp_arr = nullptr;
    int nmin_dp_arr = 0;
    int* pl_arr = nullptr;
    int npl_arr = 0;
    int* ad_arr = nullptr;
    int nad_arr = 0;
    int* end_arr = nullptr;
    int nend_arr = 0;

    while (block.variant_count < config_.block_size) {
        int ret = bcf_read(fp, hdr, rec);
        if (ret < 0) {
            break; // EOF or error
        }

        bcf_unpack(rec, BCF_UN_ALL);

        // CHROM
        const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
        block.position.chrom.push_back(chrom ? chrom : ".");

        // POS
        block.position.pos.push_back(rec->pos + 1); // Convert to 1-based

        // ID
        block.position.id.push_back(rec->d.id ? rec->d.id : ".");

        // REF
        block.sequence.ref.push_back(rec->d.allele[0] ? rec->d.allele[0] : ".");

        // ALT
        std::vector<std::string> alts;
        for (int i = 1; i < rec->n_allele; i++) {
            alts.push_back(rec->d.allele[i] ? rec->d.allele[i] : ".");
        }
        if (alts.empty()) {
            alts.push_back(".");
        }
        block.sequence.alt.push_back(alts);

        // QUAL
        block.quality.qual.push_back(bcf_float_is_missing(rec->qual) ? -1.0f : rec->qual);

        // FILTER
        if (rec->d.n_flt > 0) {
            std::string flt_str;
            for (int i = 0; i < rec->d.n_flt; i++) {
                if (i > 0) flt_str += ";";
                const char* flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
                flt_str += flt ? flt : ".";
            }
            block.quality.filter.push_back(flt_str);
        } else {
            block.quality.filter.push_back(".");
        }

        // INFO/END
        int ret_end = bcf_get_info_int32(hdr, rec, "END", &end_arr, &nend_arr);
        if (ret_end > 0 && end_arr) {
            block.info.end.push_back(end_arr[0]);
        } else {
            block.info.end.push_back(-1); // No END field
        }

        // GT
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt > 0 && gt_arr) {
            int8_t a1 = -1, a2 = -1;
            bool phased = false;

            if (ngt >= 1) {
                if (gt_arr[0] == bcf_gt_missing) {
                    a1 = -1;
                } else {
                    a1 = bcf_gt_allele(gt_arr[0]);
                }
            }
            if (ngt >= 2) {
                if (gt_arr[1] == bcf_gt_missing) {
                    a2 = -1;
                } else {
                    a2 = bcf_gt_allele(gt_arr[1]);
                    phased = bcf_gt_is_phased(gt_arr[1]);
                }
            }

            block.sample.gt.push_back(GenotypeData(a1, a2, phased));
        } else {
            block.sample.gt.push_back(GenotypeData(-1, -1, false));
        }

        // DP
        int ndp = bcf_get_format_int32(hdr, rec, "DP", &dp_arr, &ndp_arr);
        if (ndp > 0 && dp_arr && dp_arr[0] != bcf_int32_missing) {
            block.sample.dp.push_back(dp_arr[0]);
        } else {
            block.sample.dp.push_back(-1);
        }

        // GQ
        int ngq = bcf_get_format_int32(hdr, rec, "GQ", &gq_arr, &ngq_arr);
        if (ngq > 0 && gq_arr && gq_arr[0] != bcf_int32_missing) {
            block.sample.gq.push_back(gq_arr[0]);
        } else {
            block.sample.gq.push_back(-1);
        }

        // MIN_DP
        int nmin_dp = bcf_get_format_int32(hdr, rec, "MIN_DP", &min_dp_arr, &nmin_dp_arr);
        if (nmin_dp > 0 && min_dp_arr && min_dp_arr[0] != bcf_int32_missing) {
            block.sample.min_dp.push_back(min_dp_arr[0]);
        } else {
            block.sample.min_dp.push_back(-1);
        }

        // PL
        int npl = bcf_get_format_int32(hdr, rec, "PL", &pl_arr, &npl_arr);
        if (npl > 0 && pl_arr) {
            std::vector<int32_t> pl_vals;
            for (int i = 0; i < npl; i++) {
                if (pl_arr[i] != bcf_int32_missing && pl_arr[i] != bcf_int32_vector_end) {
                    pl_vals.push_back(pl_arr[i]);
                }
            }
            block.sample.pl.push_back(pl_vals);
        } else {
            block.sample.pl.push_back(std::vector<int32_t>());
        }

        // AD
        int nad = bcf_get_format_int32(hdr, rec, "AD", &ad_arr, &nad_arr);
        if (nad > 0 && ad_arr) {
            std::vector<int32_t> ad_vals;
            for (int i = 0; i < nad; i++) {
                if (ad_arr[i] != bcf_int32_missing && ad_arr[i] != bcf_int32_vector_end) {
                    ad_vals.push_back(ad_arr[i]);
                }
            }
            block.sample.ad.push_back(ad_vals);
        } else {
            block.sample.ad.push_back(std::vector<int32_t>());
        }

        block.variant_count++;
    }

    // Free allocated arrays
    if (gt_arr) free(gt_arr);
    if (dp_arr) free(dp_arr);
    if (gq_arr) free(gq_arr);
    if (min_dp_arr) free(min_dp_arr);
    if (pl_arr) free(pl_arr);
    if (ad_arr) free(ad_arr);
    if (end_arr) free(end_arr);

    return block.variant_count > 0;
}

void GVCFCompressor::CloseInput() {
    if (vcf_record_) {
        bcf_destroy1(static_cast<bcf1_t*>(vcf_record_));
        vcf_record_ = nullptr;
    }
    if (vcf_header_) {
        bcf_hdr_destroy(static_cast<bcf_hdr_t*>(vcf_header_));
        vcf_header_ = nullptr;
    }
    if (vcf_file_) {
        hts_close(static_cast<htsFile*>(vcf_file_));
        vcf_file_ = nullptr;
    }
}

std::shared_ptr<CompressionBackend> GVCFCompressor::CreateBackend() {
    // Use GSCBackendAdapter with compression strategy
    auto compress_fn = [](const std::vector<uint8_t>& input, std::vector<uint8_t>& output) -> bool {
        if (input.empty()) {
            output.clear();
            return true;
        }
        return GVCFCompressData(input.data(), input.size(), output);
    };
    auto decompress_fn = [](const std::vector<uint8_t>& input, std::vector<uint8_t>& output) -> bool {
        if (input.empty()) {
            output.clear();
            return true;
        }
        return GVCFDecompressData(input.data(), input.size(), output);
    };
    auto decompress_ptr_fn = [](const uint8_t* input, size_t size, std::vector<uint8_t>& output) -> bool {
        if (size == 0) {
            output.clear();
            return true;
        }
        return GVCFDecompressData(input, size, output);
    };
    return std::make_shared<GSCBackendAdapter>(compress_fn, decompress_fn, decompress_ptr_fn);
}

bool GVCFCompressor::WriteFileHeader() {
    // Write magic number
    fwrite(&GVCF_FILE_MAGIC, sizeof(uint32_t), 1, output_file_);

    // Write version (V3 = compressed header)
    uint32_t version = 3;
    fwrite(&version, sizeof(uint32_t), 1, output_file_);

    // Write compression backend (V2+)
    uint8_t backend_id = static_cast<uint8_t>(params_.backend);
    fwrite(&backend_id, sizeof(uint8_t), 1, output_file_);

    // Write number of samples
    fwrite(&num_samples_, sizeof(uint32_t), 1, output_file_);

    // Compress and write header text (V3+)
    std::vector<uint8_t> header_data(header_text_.begin(), header_text_.end());
    std::vector<uint8_t> compressed_header;

    EnsureGVCFStrategy();
    bool header_compressed = false;
    if (g_gvcf_strategy && g_gvcf_strategy->Compress(header_data, compressed_header)) {
        if (compressed_header.size() < header_data.size()) {
            header_compressed = true;
        }
    }

    // Write compression flag
    uint8_t header_flag = header_compressed ? 1 : 0;
    fwrite(&header_flag, sizeof(uint8_t), 1, output_file_);

    if (header_compressed) {
        // Write original size (for decompression buffer allocation)
        uint32_t original_size = static_cast<uint32_t>(header_data.size());
        fwrite(&original_size, sizeof(uint32_t), 1, output_file_);
        // Write compressed header
        uint32_t compressed_size = static_cast<uint32_t>(compressed_header.size());
        fwrite(&compressed_size, sizeof(uint32_t), 1, output_file_);
        fwrite(compressed_header.data(), 1, compressed_size, output_file_);
    } else {
        // Write uncompressed header (fallback)
        uint32_t header_size = static_cast<uint32_t>(header_text_.size());
        fwrite(&header_size, sizeof(uint32_t), 1, output_file_);
        fwrite(header_text_.data(), 1, header_size, output_file_);
    }

    // Write sample names
    for (const auto& name : sample_names_) {
        uint32_t name_size = static_cast<uint32_t>(name.size());
        fwrite(&name_size, sizeof(uint32_t), 1, output_file_);
        fwrite(name.data(), 1, name_size, output_file_);
    }

    // Reserve space for block count and offsets (will be written in footer)
    uint64_t placeholder = 0;
    fwrite(&placeholder, sizeof(uint64_t), 1, output_file_);

    return true;
}

bool GVCFCompressor::CompressAndWriteBlock(const GVCFBlock& block) {
    // Record block offset
    uint64_t offset = ftell(output_file_);
    block_offsets_.push_back(offset);

    // Compress block
    CompressedGVCFBlock compressed;
    GVCFBlockCompressor compressor(backend_, config_);

    if (!compressor.Compress(block, compressed)) {
        return false;
    }

    // Serialize to buffer
    std::vector<uint8_t> buffer;
    if (!compressed.Serialize(buffer)) {
        return false;
    }

    // Write block size and data
    uint32_t block_size = static_cast<uint32_t>(buffer.size());
    fwrite(&block_size, sizeof(uint32_t), 1, output_file_);
    fwrite(buffer.data(), 1, buffer.size(), output_file_);

    // Update statistics
    stats_.chrom_compressed += compressed.chrom.data.size();
    stats_.pos_compressed += compressed.pos.data.size();
    stats_.gt_compressed += compressed.gt_mask.data.size() + compressed.gt_patches.data.size();
    stats_.dp_compressed += compressed.dp.data.size();
    stats_.gq_compressed += compressed.gq.data.size();
    stats_.pl_compressed += compressed.pl.data.size();

    // Debug: print all field sizes
    auto logger = LogManager::Instance().Logger();
    logger->debug("Field sizes: CHROM={} POS={} ID={} REF={} ALT={} QUAL={} FILTER={} END={} GT={} DP={} GQ={} MIN_DP={} PL={} AD={}",
        compressed.chrom.data.size(), compressed.pos.data.size(),
        compressed.id.data.size(), compressed.ref.data.size(),
        compressed.alt.data.size(), compressed.qual.data.size(),
        compressed.filter.data.size(), compressed.info_end.data.size(),
        compressed.gt_mask.data.size() + compressed.gt_patches.data.size(),
        compressed.dp.data.size(), compressed.gq.data.size(),
        compressed.min_dp.data.size(), compressed.pl.data.size(),
        compressed.ad.data.size());

    return true;
}

bool GVCFCompressor::WriteFileFooter() {
    // Get footer position
    uint64_t footer_offset = ftell(output_file_);

    // Write block count
    uint32_t num_blocks = static_cast<uint32_t>(block_offsets_.size());
    fwrite(&num_blocks, sizeof(uint32_t), 1, output_file_);

    // Write block offsets
    fwrite(block_offsets_.data(), sizeof(uint64_t), block_offsets_.size(), output_file_);

    // Write total variant count
    fwrite(&stats_.total_variants, sizeof(uint64_t), 1, output_file_);

    // Write footer offset at the end
    fwrite(&footer_offset, sizeof(uint64_t), 1, output_file_);

    return true;
}

// ============================================================================
// GVCFDecompressor Implementation
// ============================================================================

GVCFDecompressor::GVCFDecompressor(const GSC_Params& params)
    : params_(params)
    , input_file_(nullptr)
    , output_file_(nullptr)
    , output_header_(nullptr)
    , num_blocks_(0)
{
}

GVCFDecompressor::~GVCFDecompressor() {
    CloseInput();
}

bool GVCFDecompressor::Decompress() {
    auto logger = LogManager::Instance().Logger();

    // Use GSCBackendAdapter with compression strategy (matches compressor)
    auto decompress_fn = [](const std::vector<uint8_t>& input, std::vector<uint8_t>& output) -> bool {
        if (input.empty()) {
            output.clear();
            return true;
        }
        return GVCFDecompressData(input.data(), input.size(), output);
    };
    auto decompress_ptr_fn = [](const uint8_t* input, size_t size, std::vector<uint8_t>& output) -> bool {
        if (size == 0) {
            output.clear();
            return true;
        }
        return GVCFDecompressData(input, size, output);
    };
    backend_ = std::make_shared<GSCBackendAdapter>(nullptr, decompress_fn, decompress_ptr_fn);

    // Open input file
    if (!OpenInput()) {
        logger->error("Failed to open input file: {}", params_.in_file_name);
        return false;
    }

    // Read file header
    if (!ReadFileHeader()) {
        logger->error("Failed to read file header");
        return false;
    }

    logger->info("Starting gVCF decompression: {} blocks", num_blocks_);

    // Determine output mode based on file extension
    std::string output_mode = "w";  // Default: plain VCF
    const std::string& out_name = params_.out_file_name;
    if (params_.out_type == file_type::BCF_File ||
        (out_name.size() > 4 && out_name.substr(out_name.size() - 4) == ".bcf")) {
        output_mode = "wb";  // BCF
    } else if ((out_name.size() > 7 && out_name.substr(out_name.size() - 7) == ".vcf.gz") ||
               (out_name.size() > 3 && out_name.substr(out_name.size() - 3) == ".gz")) {
        output_mode = "wz";  // gzip compressed VCF
    }
    logger->info("Output format: {}", output_mode == "wb" ? "BCF" : (output_mode == "wz" ? "VCF.GZ" : "VCF"));

    // Open output file
    htsFile* out_fp = hts_open(params_.out_file_name.c_str(), output_mode.c_str());
    if (!out_fp) {
        logger->error("Failed to open output file: {}", params_.out_file_name);
        return false;
    }
    output_file_ = out_fp;

    // Create and write header - parse header text line by line
    bcf_hdr_t* hdr = bcf_hdr_init("w");

    // Parse header lines
    std::istringstream iss(header_text_);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
        if (line[0] == '#' && line[1] == '#') {
            // Meta-info line
            bcf_hdr_append(hdr, line.c_str());
        } else if (line[0] == '#') {
            // Header line - extract sample names
            // #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE...
            // Skip standard columns and add samples
            std::istringstream header_iss(line);
            std::string col;
            int col_idx = 0;
            while (std::getline(header_iss, col, '\t')) {
                if (col_idx >= 9) {
                    // This is a sample name
                    bcf_hdr_add_sample(hdr, col.c_str());
                }
                col_idx++;
            }
        }
    }

    // Sync the header
    if (bcf_hdr_sync(hdr) < 0) {
        logger->error("Failed to sync header");
        bcf_hdr_destroy(hdr);
        return false;
    }

    output_header_ = hdr;

    if (bcf_hdr_write(out_fp, hdr) < 0) {
        logger->error("Failed to write VCF header");
        return false;
    }

    // Process blocks
    GVCFBlock block;
    // Note: decompressor is created inside ReadBlock, not here

    for (uint32_t i = 0; i < num_blocks_; i++) {
        logger->debug("Processing block {}/{}", i, num_blocks_);
        if (!ReadBlock(i, block)) {
            logger->error("Failed to read block {}", i);
            return false;
        }
        logger->debug("Block {} read successfully, writing {} variants", i, block.variant_count);

        // Write all variants in block
        for (uint32_t j = 0; j < block.variant_count; j++) {
            logger->debug("Writing variant {}/{}", j, block.variant_count);
            if (!WriteVCFRecord(block, j)) {
                logger->error("Failed to write variant {}", j);
                return false;
            }
        }
        logger->debug("Block {} complete", i);

        if ((i + 1) % 100 == 0) {
            logger->debug("Decompressed {} / {} blocks", i + 1, num_blocks_);
        }
    }

    logger->info("gVCF decompression complete");

    // Cleanup
    bcf_hdr_destroy(static_cast<bcf_hdr_t*>(output_header_));
    hts_close(static_cast<htsFile*>(output_file_));
    output_file_ = nullptr;
    output_header_ = nullptr;

    return true;
}

bool GVCFDecompressor::OpenInput() {
    input_file_ = fopen(params_.in_file_name.c_str(), "rb");
    return input_file_ != nullptr;
}

bool GVCFDecompressor::ReadFileHeader() {
    auto logger = LogManager::Instance().Logger();

    // Read and verify magic number
    uint32_t magic;
    if (fread(&magic, sizeof(uint32_t), 1, input_file_) != 1) {
        return false;
    }
    if (magic != GVCF_FILE_MAGIC) {
        return false;
    }

    // Read version
    uint32_t version;
    if (fread(&version, sizeof(uint32_t), 1, input_file_) != 1) {
        return false;
    }
    if (version > 3) {  // Support up to V3
        logger->error("Unsupported file version: {}", version);
        return false;
    }

    // Read compression backend (V2+)
    compression_backend_t file_backend = compression_backend_t::zstd;  // V1 default
    if (version >= 2) {
        uint8_t backend_id;
        if (fread(&backend_id, sizeof(uint8_t), 1, input_file_) != 1) {
            return false;
        }
        if (backend_id <= static_cast<uint8_t>(compression_backend_t::brotli)) {
            file_backend = static_cast<compression_backend_t>(backend_id);
        }
    }

    // Set backend for decompression
    SetGVCFBackend(file_backend);
    InitializeCompressionBackend(file_backend);
    logger->info("Using compression backend: {}",
        file_backend == compression_backend_t::bsc ? "bsc" :
        file_backend == compression_backend_t::zstd ? "zstd" : "brotli");

    // Read number of samples
    uint32_t num_samples;
    if (fread(&num_samples, sizeof(uint32_t), 1, input_file_) != 1) {
        return false;
    }

    // Read header text (V3: may be compressed)
    if (version >= 3) {
        // V3 format: [flag] [original_size] [compressed_size] [data] or [flag] [size] [data]
        uint8_t header_flag;
        if (fread(&header_flag, sizeof(uint8_t), 1, input_file_) != 1) {
            return false;
        }

        if (header_flag == 1) {
            // Compressed header
            uint32_t original_size;
            if (fread(&original_size, sizeof(uint32_t), 1, input_file_) != 1) {
                return false;
            }
            uint32_t compressed_size;
            if (fread(&compressed_size, sizeof(uint32_t), 1, input_file_) != 1) {
                return false;
            }
            std::vector<uint8_t> compressed_data(compressed_size);
            if (fread(compressed_data.data(), 1, compressed_size, input_file_) != compressed_size) {
                return false;
            }

            // Decompress
            EnsureGVCFStrategy();
            std::vector<uint8_t> decompressed_data;
            if (!g_gvcf_strategy || !g_gvcf_strategy->Decompress(compressed_data, decompressed_data)) {
                logger->error("Failed to decompress header");
                return false;
            }
            header_text_.assign(decompressed_data.begin(), decompressed_data.end());
        } else {
            // Uncompressed header
            uint32_t header_size;
            if (fread(&header_size, sizeof(uint32_t), 1, input_file_) != 1) {
                return false;
            }
            header_text_.resize(header_size);
            if (fread(&header_text_[0], 1, header_size, input_file_) != header_size) {
                return false;
            }
        }
    } else {
        // V1/V2 format: uncompressed header
        uint32_t header_size;
        if (fread(&header_size, sizeof(uint32_t), 1, input_file_) != 1) {
            return false;
        }
        header_text_.resize(header_size);
        if (fread(&header_text_[0], 1, header_size, input_file_) != header_size) {
            return false;
        }
    }

    // Read sample names
    sample_names_.clear();
    for (uint32_t i = 0; i < num_samples; i++) {
        uint32_t name_size;
        if (fread(&name_size, sizeof(uint32_t), 1, input_file_) != 1) {
            return false;
        }
        std::string name(name_size, '\0');
        if (fread(&name[0], 1, name_size, input_file_) != name_size) {
            return false;
        }
        sample_names_.push_back(name);
    }

    // Skip placeholder
    uint64_t placeholder;
    if (fread(&placeholder, sizeof(uint64_t), 1, input_file_) != 1) {
        return false;
    }

    // Read footer (seek to end, read footer offset, then read footer)
    fseek(input_file_, -static_cast<long>(sizeof(uint64_t)), SEEK_END);
    uint64_t footer_offset;
    if (fread(&footer_offset, sizeof(uint64_t), 1, input_file_) != 1) {
        return false;
    }

    fseek(input_file_, footer_offset, SEEK_SET);

    // Read block count
    if (fread(&num_blocks_, sizeof(uint32_t), 1, input_file_) != 1) {
        return false;
    }

    // Read block offsets
    block_offsets_.resize(num_blocks_);
    if (fread(block_offsets_.data(), sizeof(uint64_t), num_blocks_, input_file_) != num_blocks_) {
        return false;
    }

    return true;
}

bool GVCFDecompressor::ReadBlock(uint32_t block_id, GVCFBlock& block) {
    auto logger = LogManager::Instance().Logger();

    if (block_id >= block_offsets_.size()) {
        logger->error("Invalid block_id {} (num_blocks={})", block_id, block_offsets_.size());
        return false;
    }

    logger->debug("Reading block {} at offset {}", block_id, block_offsets_[block_id]);
    fseek(input_file_, block_offsets_[block_id], SEEK_SET);

    // Read block size
    uint32_t block_size;
    if (fread(&block_size, sizeof(uint32_t), 1, input_file_) != 1) {
        logger->error("Failed to read block size");
        return false;
    }
    logger->debug("Block size: {}", block_size);

    // Sanity check block size
    if (block_size > 100 * 1024 * 1024) { // 100MB limit
        logger->error("Block size too large: {}", block_size);
        return false;
    }

    // Read block data
    std::vector<uint8_t> buffer(block_size);
    if (fread(buffer.data(), 1, block_size, input_file_) != block_size) {
        logger->error("Failed to read block data");
        return false;
    }

    // Deserialize compressed block
    CompressedGVCFBlock compressed;
    logger->debug("Deserializing compressed block...");
    if (!compressed.Deserialize(buffer)) {
        logger->error("Failed to deserialize compressed block");
        return false;
    }
    logger->debug("Deserialized: {} variants", compressed.variant_count);

    // Decompress
    logger->debug("Decompressing block...");
    GVCFBlockDecompressor decompressor(backend_);
    bool result = decompressor.Decompress(compressed, block);
    if (!result) {
        logger->error("Failed to decompress block");
    } else {
        logger->debug("Block decompressed successfully, variant_count={}", block.variant_count);
    }
    return result;
}

bool GVCFDecompressor::WriteVCFRecord(const GVCFBlock& block, uint32_t idx) {
    auto logger = LogManager::Instance().Logger();

    htsFile* out_fp = static_cast<htsFile*>(output_file_);
    bcf_hdr_t* hdr = static_cast<bcf_hdr_t*>(output_header_);

    if (!out_fp || !hdr) {
        logger->error("Output file or header is null");
        return false;
    }

    bcf1_t* rec = bcf_init1();
    if (!rec) {
        logger->error("Failed to init bcf1_t");
        return false;
    }

    logger->debug("  Setting CHROM: {}", block.position.chrom[idx]);
    // Set CHROM and POS
    rec->rid = bcf_hdr_name2id(hdr, block.position.chrom[idx].c_str());
    if (rec->rid < 0) {
        logger->warn("Unknown CHROM: {}, adding to header", block.position.chrom[idx]);
        // Try to add the contig to header
        std::string contig_line = "##contig=<ID=" + block.position.chrom[idx] + ">";
        bcf_hdr_append(hdr, contig_line.c_str());
        bcf_hdr_sync(hdr);
        rec->rid = bcf_hdr_name2id(hdr, block.position.chrom[idx].c_str());
    }

    logger->debug("  Setting POS: {}", block.position.pos[idx]);
    rec->pos = block.position.pos[idx] - 1; // Convert to 0-based

    // Set ID
    logger->debug("  Setting ID: {}", block.position.id[idx]);
    bcf_update_id(hdr, rec, block.position.id[idx].c_str());

    // Set alleles (REF + ALT)
    logger->debug("  Setting REF: {}", block.sequence.ref[idx]);
    std::vector<const char*> alleles;
    alleles.push_back(block.sequence.ref[idx].c_str());
    logger->debug("  Setting ALT count: {}", block.sequence.alt[idx].size());
    for (const auto& alt : block.sequence.alt[idx]) {
        alleles.push_back(alt.c_str());
    }
    bcf_update_alleles(hdr, rec, alleles.data(), alleles.size());

    // Set QUAL
    if (block.quality.qual[idx] >= 0) {
        rec->qual = block.quality.qual[idx];
    }

    // Set FILTER
    if (!block.quality.filter[idx].empty() && block.quality.filter[idx] != ".") {
        int flt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, block.quality.filter[idx].c_str());
        if (flt_id >= 0) {
            bcf_update_filter(hdr, rec, &flt_id, 1);
        }
    }

    // Set INFO/END
    if (idx < block.info.end.size() && block.info.end[idx] > 0) {
        int32_t end_val = static_cast<int32_t>(block.info.end[idx]);
        bcf_update_info_int32(hdr, rec, "END", &end_val, 1);
    }

    // Set GT
    if (idx < block.sample.gt.size()) {
        const auto& gt = block.sample.gt[idx];
        int32_t gt_arr[2];
        if (gt.allele1 < 0) {
            gt_arr[0] = bcf_gt_missing;
        } else {
            gt_arr[0] = bcf_gt_unphased(gt.allele1);
        }
        if (gt.allele2 < 0) {
            gt_arr[1] = bcf_gt_missing;
        } else {
            gt_arr[1] = gt.phased ? bcf_gt_phased(gt.allele2) : bcf_gt_unphased(gt.allele2);
        }
        bcf_update_genotypes(hdr, rec, gt_arr, 2);
    }

    // Set DP
    if (idx < block.sample.dp.size() && block.sample.dp[idx] >= 0) {
        int32_t dp = block.sample.dp[idx];
        bcf_update_format_int32(hdr, rec, "DP", &dp, 1);
    }

    // Set GQ
    if (idx < block.sample.gq.size() && block.sample.gq[idx] >= 0) {
        int32_t gq = block.sample.gq[idx];
        bcf_update_format_int32(hdr, rec, "GQ", &gq, 1);
    }

    // Set MIN_DP
    if (idx < block.sample.min_dp.size() && block.sample.min_dp[idx] >= 0) {
        int32_t min_dp = block.sample.min_dp[idx];
        bcf_update_format_int32(hdr, rec, "MIN_DP", &min_dp, 1);
    }

    // Set PL
    if (idx < block.sample.pl.size() && !block.sample.pl[idx].empty()) {
        bcf_update_format_int32(hdr, rec, "PL",
                               block.sample.pl[idx].data(),
                               block.sample.pl[idx].size());
    }

    // Set AD
    if (idx < block.sample.ad.size() && !block.sample.ad[idx].empty()) {
        bcf_update_format_int32(hdr, rec, "AD",
                               block.sample.ad[idx].data(),
                               block.sample.ad[idx].size());
    }

    // Write record
    int ret = bcf_write(out_fp, hdr, rec);
    bcf_destroy1(rec);

    return ret >= 0;
}

void GVCFDecompressor::CloseInput() {
    if (input_file_) {
        fclose(input_file_);
        input_file_ = nullptr;
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

bool IsGVCFFile(const std::string& filename) {
    htsFile* fp = hts_open(filename.c_str(), "r");
    if (!fp) return false;

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        return false;
    }

    // Check for gVCF indicators
    bool has_end = (bcf_hdr_id2int(hdr, BCF_DT_ID, "END") >= 0);
    bool has_min_dp = (bcf_hdr_id2int(hdr, BCF_DT_ID, "MIN_DP") >= 0);

    // Check header for NON_REF
    kstring_t str = {0, 0, nullptr};
    bcf_hdr_format(hdr, 0, &str);
    bool has_non_ref = false;
    if (str.s) {
        has_non_ref = (strstr(str.s, "NON_REF") != nullptr);
        free(str.s);
    }

    bcf_hdr_destroy(hdr);
    hts_close(fp);

    return has_end || has_min_dp || has_non_ref;
}

bool IsGVCFCompressed(const std::string& filename) {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) return false;

    uint32_t magic;
    bool result = false;

    if (fread(&magic, sizeof(uint32_t), 1, fp) == 1) {
        result = (magic == GVCF_FILE_MAGIC);
    }

    fclose(fp);
    return result;
}

} // namespace gvcf
