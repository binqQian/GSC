#include "decompressor.h"
#include "logger.h"
#include <atomic>
#include <bitset>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <cstdlib>
#include <limits>
#include <memory>
#include <mutex>
#include <utility>
#include <zlib.h>

using namespace std::chrono;

namespace {

template <typename F>
class ScopeExit
{
    F fn_;
    bool active_ = true;

public:
    explicit ScopeExit(F fn) : fn_(std::move(fn)) {}
    ~ScopeExit()
    {
        if (active_)
            fn_();
    }
    void dismiss() { active_ = false; }
};

template <typename F>
ScopeExit<F> MakeScopeExit(F fn)
{
    return ScopeExit<F>(std::move(fn));
}

bool ends_with(const std::string &value, const std::string &suffix) {
    return value.size() >= suffix.size() &&
           std::equal(suffix.rbegin(), suffix.rend(), value.rbegin());
}

std::string vcf_write_mode(const std::string &path, char compression_level) {
    // Use bgzip when the caller asks for a gzipped VCF filename.
    bool wants_gzip = path != "-" &&
                      (ends_with(path, ".vcf.gz") || ends_with(path, ".vcf.bgz") ||
                       ends_with(path, ".vcf.bgzf") || ends_with(path, ".gz"));
    if (!wants_gzip)
        return "w";

    if (compression_level == 'u')
        return "w";

    std::string mode = "wz";
    if (compression_level >= '0' && compression_level <= '9')
        mode.push_back(compression_level);
    return mode;
}

bool use_bgzf_threads(const std::string& vcf_mode, uint32_t threads) {
    return threads > 1 && vcf_mode.rfind("wz", 0) == 0;
}

void append_u64(std::string& out, uint64_t value) {
    char buf[32];
    auto result = std::to_chars(buf, buf + sizeof(buf), value);
    if (result.ec == std::errc()) {
        out.append(buf, static_cast<size_t>(result.ptr - buf));
    }
}

void append_first_alt(std::string& out, const std::string& alt) {
    size_t comma_pos = alt.find(',');
    if (comma_pos == std::string::npos)
        out.append(alt);
    else
        out.append(alt.data(), comma_pos);
}

void build_bim_line(std::string& out, const variant_desc_t& desc, const std::string& genetic_distance) {
    out.clear();
    out.reserve(desc.chrom.size() + desc.id.size() + desc.ref.size() + desc.alt.size() +
                genetic_distance.size() + 32);
    out.append(desc.chrom);
    out.push_back('\t');
    out.append(desc.id);
    out.push_back('\t');
    out.append(genetic_distance);
    out.push_back('\t');
    append_u64(out, desc.pos);
    out.push_back('\t');
    append_first_alt(out, desc.alt);
    out.push_back('\t');
    out.append(desc.ref);
    out.push_back('\n');
}

void build_pvar_line(std::string& out, const variant_desc_t& desc) {
    out.clear();
    out.reserve(desc.chrom.size() + desc.id.size() + desc.ref.size() + desc.alt.size() + 32);
    out.append(desc.chrom);
    out.push_back('\t');
    append_u64(out, desc.pos);
    out.push_back('\t');
    if (desc.id.empty())
        out.push_back('.');
    else
        out.append(desc.id);
    out.push_back('\t');
    out.append(desc.ref);
    out.push_back('\t');
    out.append(desc.alt);
    out.push_back('\n');
}

bool patch_u32_le_at(const std::string& path, long offset, uint32_t value) {
    FILE* f = fopen(path.c_str(), "r+b");
    if (!f)
        return false;
    if (fseek(f, offset, SEEK_SET) != 0) {
        fclose(f);
        return false;
    }
    const size_t n = fwrite(&value, sizeof(value), 1, f);
    fclose(f);
    return n == 1;
}

bool is_biallelic_desc_for_export(const variant_desc_t& desc) {
    // Skip GSC multi-allelic split markers (<N>/<M>) and general multi-allelic ALT lists.
    // BGEN/PGEN outputs here are biallelic-only.
    if (desc.alt.find('<') != std::string::npos)
        return false;
    if (desc.alt.find(',') != std::string::npos)
        return false;
    return true;
}

template <typename T>
class BoundedQueue
{
    std::mutex mtx_;
    std::condition_variable cv_pop_;
    std::condition_variable cv_push_;
    std::deque<T> q_;
    size_t max_size_ = 0;
    bool done_ = false;

public:
    explicit BoundedQueue(size_t max_size) : max_size_(max_size) {}

    void Complete()
    {
        std::lock_guard<std::mutex> lock(mtx_);
        done_ = true;
        cv_pop_.notify_all();
        cv_push_.notify_all();
    }

    bool Push(T &&item)
    {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_push_.wait(lock, [&] { return done_ || q_.size() < max_size_; });
        if (done_)
            return false;
        q_.emplace_back(std::move(item));
        cv_pop_.notify_one();
        return true;
    }

    bool Pop(T &out)
    {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_pop_.wait(lock, [&] { return done_ || !q_.empty(); });
        if (q_.empty())
            return false;
        out = std::move(q_.front());
        q_.pop_front();
        cv_push_.notify_one();
        return true;
    }
};

struct DecompressedChunk
{
    uint32_t chunk_id = 0; // 0-based chunk id
    uint32_t start_actual_pos = 0;
    uint32_t end_actual_pos = 0;
    uint32_t chunk_variant_offset = 0;
    bool has_fixed_fields_rb_dir = false;
    uint32_t fixed_fields_chunk_version = 0;

    std::vector<block_t> fixed_variants_chunk;
    std::vector<std::vector<std::vector<uint32_t>>> sort_perm;
    std::vector<uint8_t> gt_indexes;
    std::vector<std::pair<uint32_t, uint32_t>> row_block_ranges;
    std::vector<std::vector<field_desc>> all_fields;
};

} // namespace




bool Decompressor::analyzeInputRange(uint32_t & start_chunk_id,uint32_t & end_chunk_id){

    size_t prev_chrom_size;
    bool query_flag = false;
    auto logger = LogManager::Instance().Logger();

    if (range == "")
    {
        // Use the actual chunk count from chunks_min_pos, which is read from the compressed file
        end_chunk_id = static_cast<uint32_t>(decompression_reader.chunks_min_pos.size());
    }
    else
    {
        std::regex pattern(R"(^([\w.]+)(?::(-?\d+))?(?:,(-?\d+))?:?,?$)");
        std::smatch matches;
        std::string cur_query_chrom;
        
        if (std::regex_match(range, matches, pattern)) {
            cur_query_chrom = matches[1].str();
            
            if(matches[2].matched){
                range_1 = std::stoi(matches[2].str());
                if(matches[3].matched)
                    range_2 = std::stoi(matches[3].str());
                else
                    range_2 = MAX;
            }else{
                
                if(matches[3].matched){
                    logger->error("Invalid input format.");
                    return 0;
                }else{
                    range_1 = 0;
                    range_2 = MAX;
                }
            }
            
            // range_1 = (matches[2].matched) ? std::stoi(matches[2].str()) : 0 ;
            // range_2 = (matches[3].matched) ? std::stoi(matches[3].str()) : MAX;

            // std::cerr << "Chromosome: " << cur_query_chrom << std::endl;
            // std::cerr << "Start position: " << range_1 << std::endl;
            // std::cerr << "End position: " << range_2 << std::endl;
        } else {
            
            logger->error("Invalid input format.");

        }
        
        // string cur_query_chrom = range.substr(0, range.find(':'));
        // auto curr_pos = range.find(':');
        // string_view cur_query_chrom(range.c_str(),curr_pos);
        
        // try {
        //     string cur_range = range.substr(curr_pos + 1);
            
        //     curr_pos = cur_range.find(',');
        //     if(curr_pos == string::npos){
                
        //         if(cur_range == ""){
        //             range_1 = 0;
        //             range_2 = MAX;
        //         }
        //         else{
        //             range_1 = stoll(cur_range);
        //             range_2 = MAX;
                    
        //         }
        //     }
        //     else{
        //         range_1 = stoll(cur_range.substr(0, curr_pos));
        //         cur_range = cur_range.substr(curr_pos + 1);
        //         if (cur_range == "")
        //             range_2 = MAX;    
        //         else
        //         {
        //             range_2 = stoll(cur_range);
        //         }    
                           
        //     }
        // } 
        // catch (const std::invalid_argument& e) {
        //     std::cerr << "Invalid argument: " << e.what() << std::endl;
        // }
        //  catch (const std::out_of_range& e) {
        //     std::cerr << "Out of range: " << e.what() << std::endl;
        // }
        if (range_2 < range_1)
        {
            logger->error("Invalid range selection.");
            return 0;
        }

        for (size_t i = 0; i < decompression_reader.d_where_chrom.size(); i++)
        {
            prev_chrom_size = i ? decompression_reader.d_where_chrom[i - 1].second : 0;
            
            if (cur_query_chrom == decompression_reader.d_where_chrom[i].first)
            {
                query_flag = true;
                start_chunk_id = end_chunk_id;
                end_chunk_id += (decompression_reader.d_where_chrom[i].second - prev_chrom_size + params.no_blocks-1) / params.no_blocks;

                break;
            }
            end_chunk_id += (decompression_reader.d_where_chrom[i].second - prev_chrom_size + params.no_blocks-1) / params.no_blocks;
            
        }
        
        if (!query_flag)
        {
            logger->error("The specified chromosome was not found.");
            return 0;
        }

        int64_t start_chunk = lower_bound(decompression_reader.chunks_min_pos.begin() + start_chunk_id, decompression_reader.chunks_min_pos.begin() + end_chunk_id, range_1) - decompression_reader.chunks_min_pos.begin() - 1;
        int64_t end_chunk = upper_bound(decompression_reader.chunks_min_pos.begin() + start_chunk_id, decompression_reader.chunks_min_pos.begin() + end_chunk_id, range_2) - decompression_reader.chunks_min_pos.begin() - 1;

        if (end_chunk < 0)
        {
            logger->warn("No variants in the requested range.");
            return 0;
        }
        start_chunk_id = start_chunk > (int)start_chunk_id ? start_chunk : start_chunk_id;
        end_chunk_id = end_chunk + 1;

       
    }
    return 1;
}
bool Decompressor::initDecompression(DecompressionReader &decompression_reader){
    auto logger = LogManager::Instance().Logger();
    logger->debug("initDecompression: entering");

    haplotype_count = decompression_reader.n_samples * (uint32_t)decompression_reader.ploidy;
    row_block_size = decompression_reader.max_block_rows ? decompression_reader.max_block_rows : haplotype_count;
    standard_block_size = row_block_size;
    if (!params.MB_memory)
    {
        max_stored_unique = 0;
    }
    else if (params.max_MB_memory)
    {
        if (decompression_reader.vec_len)
            max_stored_unique = (static_cast<uint64_t>(params.max_MB_memory) * 1000000ull) /
                                static_cast<uint64_t>(decompression_reader.vec_len);
        else
            max_stored_unique = 0;
    }
    else
    {
        max_stored_unique = standard_block_size * 2;
    }

    logger->debug("initDecompression: haplotype_count={}, row_block_size={}, useLegacyPath={}, vec_len={}",
                 haplotype_count, row_block_size, decompression_reader.useLegacyPath, decompression_reader.vec_len);

    if(row_block_size < 1024)
    {
        chunk_size = CHUNK_SIZE1;
    }
    else if(row_block_size < 4096 )
    {
        chunk_size = CHUNK_SIZE2;
    }
    else if(row_block_size < 8192)
    {
        chunk_size = CHUNK_SIZE3;
    }
    else
    {
        chunk_size = row_block_size;
    }
    params.no_blocks = chunk_size / row_block_size;
    rrr_rank_zeros_bit_vector[0] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_zeros_bit_vector[0]);
    rrr_rank_zeros_bit_vector[1] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_zeros_bit_vector[1]);
    rrr_rank_copy_bit_vector[0] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_copy_bit_vector[0]);
    rrr_rank_copy_bit_vector[1] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_copy_bit_vector[1]);

    // Lossless mode always uses full GT decoding; lossy sample/range paths may not.
    bool needs_full_decomp = (params.compress_mode == compress_mode_t::lossless_mode) || params.samples.empty();
    logger->debug("initDecompression: needs_full_decomp={}, params.samples empty={}", needs_full_decomp, params.samples.empty());
    if (needs_full_decomp){
        logger->debug("initDecompression: allocating decomp_data with size={}", decompression_reader.vec_len * 2);
        decomp_data = new uint8_t[decompression_reader.vec_len * 2];
        decomp_data_perm = new uint8_t[decompression_reader.vec_len * 2];
        zeros_only_vector = new uint8_t[decompression_reader.vec_len]();
        logger->debug("initDecompression: decomp_data allocated successfully");
    }
    int no_haplotypes = haplotype_count;
    if (params.samples != ""){
        no_haplotypes = smpl.no_samples * decompression_reader.ploidy;

        // For legacy path with sample subset, we need bit vectors
        if (decompression_reader.useLegacyPath) {
            uint64_t bv_size = decompression_reader.rrr_zeros_bit_vector[0].size();

            zeros_bit_vector[0] = sdsl::bit_vector(bv_size);
            zeros_bit_vector[1] = sdsl::bit_vector(bv_size);
            copy_bit_vector[0] = sdsl::bit_vector(bv_size);
            copy_bit_vector[1] = sdsl::bit_vector(bv_size);
            uint64_t v_pos;

            for (v_pos = 0; v_pos + 64 < bv_size; v_pos += 64)
            {
                zeros_bit_vector[0].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[0].get_int(v_pos, 64), 64);
                zeros_bit_vector[1].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[1].get_int(v_pos, 64), 64);
                copy_bit_vector[0].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[0].get_int(v_pos, 64), 64);
                copy_bit_vector[1].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[1].get_int(v_pos, 64), 64);
            }

            uint64_t tail_len = bv_size - v_pos;
            zeros_bit_vector[0].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[0].get_int(v_pos, tail_len), tail_len);
            zeros_bit_vector[1].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[1].get_int(v_pos, tail_len), tail_len);
            copy_bit_vector[0].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[0].get_int(v_pos, tail_len), tail_len);
            copy_bit_vector[1].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[1].get_int(v_pos, tail_len), tail_len);
        }
    }

    // Always set up full_byte_count and trailing_bits (needed for GT reconstruction in tiled mode)
    if(haplotype_count & 7)//%8)
    {
        full_byte_count = decompression_reader.vec_len - 1;
        trailing_bits = haplotype_count & 7;//%8;
    }
    else
    {
        full_byte_count = decompression_reader.vec_len;
        trailing_bits = 0;
    }

    // Allocate tmp_arr for full decompression (both full decomp and tiled sample subset)
    if (needs_full_decomp) {
        logger->debug("initDecompression: allocating tmp_arr with size={}", full_byte_count);
        tmp_arr = new long long[full_byte_count];
        logger->debug("initDecompression: tmp_arr allocated successfully");
    }
    logger->debug("initDecompression: calling initialXORLut");
    initialXORLut();
    logger->debug("initDecompression: calling initialLut");
    initialLut();
    logger->debug("initDecompression: lut initialization done");

    max_col_vec_len_ = 0;
    col_decomp_perm_buf_.clear();
    col_decomp_buf_.clear();
    if (!decompression_reader.useLegacyPath)
    {
        for (auto len : decompression_reader.col_block_vec_lens)
            if (len > max_col_vec_len_)
                max_col_vec_len_ = len;
        if (max_col_vec_len_ > 0)
        {
            const size_t buf_len = static_cast<size_t>(max_col_vec_len_) * 2;
            col_decomp_perm_buf_.resize(buf_len);
            col_decomp_buf_.resize(buf_len);
        }
    }

    if (out_type == file_type::PGEN_File && !gsc_to_pgen_lut_ready_)
        initGscToPgenLUT();
    if (out_type == file_type::BGEN_File && !gsc_to_bgen_lut_ready_)
        initGscToBgenLUT();

    if(out_type == file_type::VCF_File || out_type == file_type::BCF_File){
        int fmt_id = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"GT");
        bcf_enc_int1(&str, fmt_id);
        bcf_enc_size(&str, decompression_reader.ploidy, BCF_BT_INT8);
        char *tmp;
        str.m = str.l + no_haplotypes + 1;
            
        kroundup32(str.m);

        if ((tmp = (char*)realloc(str.s, str.m)))
            str.s = tmp;
        else
        {
            logger->error("initDecompression: realloc failed for GT header buffer");
            return false;
        }
        
        str.l = 3;
    }

    logger->debug("initDecompression: completed successfully");
    return true;
}
//*************************************************************************************************************************************
// Official Decompress Program Entry
bool Decompressor::decompressProcess()
{
    decompression_reader.SetNoThreads(params.no_threads);
    // unique_ptr<CompressedFileLoading> cfile(new CompressedFileLoading());

    if(params.compress_mode == compress_mode_t::lossless_mode){
        decompression_mode_type = true;
    }
    // unique_ptr<CompressedFileLoading> cfile(new CompressedFileLoading());
    if(!decompression_reader.OpenReading(in_file_name, decompression_mode_type))
        return false;

    auto cleanup_guard = MakeScopeExit([&] {
        decompression_reader.close();
        Close();
    });

    // Reject mode mismatches early: on-disk mode is stored in header as `mode_type`.
    if (decompression_reader.file_mode_type && params.compress_mode == compress_mode_t::lossly_mode)
    {
        auto logger = LogManager::Instance().Logger();
        logger->error("Input `.gsc` is lossless (Mode=true); do not use `-M/--mode_lossly` for decompression.");
        return false;
    }
    if (!decompression_reader.file_mode_type && params.compress_mode == compress_mode_t::lossless_mode)
    {
        auto logger = LogManager::Instance().Logger();
        logger->error("Input `.gsc` is lossy (Mode=false); please add `-M/--mode_lossly` for decompression.");
        return false;
    }
    if (params.compress_mode == compress_mode_t::lossless_mode && !params.samples.empty())
    {
        auto logger = LogManager::Instance().Logger();
        logger->info("Lossless sample query enabled: output fixed fields + GT only (skips other fields decoding).");
    }
    if(decompression_mode_type) {
        const bool start_threads = params.samples.empty();
        if(!decompression_reader.OpenReadingPart2(in_file_name, start_threads))
            return false;
        
    }    
    decompression_reader.decompress_meta(v_samples, header);
    
    if(analyzeInputSamples(v_samples)) // Retrieving sample name.
        return false;

    if (params.header_only)
    {
        if (params.split_flag)
        {
            auto logger = LogManager::Instance().Logger();
            logger->error("--header-only is not compatible with --split");
            return false;
        }
        if (!(out_type == file_type::VCF_File || out_type == file_type::BCF_File))
        {
            auto logger = LogManager::Instance().Logger();
            logger->error("--header-only is only supported for VCF/BCF output");
            return false;
        }
    }


    if(params.split_flag){

        if (!(out_type == file_type::VCF_File || out_type == file_type::BCF_File))
        {
            auto logger = LogManager::Instance().Logger();
            logger->error("Split output is only supported for VCF/BCF.");
            return false;
        }

        if(!splitFileWriting(static_cast<int>(decompression_reader.d_where_chrom.size())))
            return false;
        initOutSplitFile();
    }else{
        if(!OpenForWriting())
            return false;
        initOut();
    }

    if (params.header_only)
    {
        // Header was written by initOut(). Exit early without decoding records.
        return true;
    }

    if (!initDecompression(decompression_reader))
        return false;
    
    if(!analyzeInputRange(start_chunk_id,end_chunk_id)){
        return false;
    }
    
    if(!decompression_reader.setStartChunk(start_chunk_id)){
        return false;
    } 
    std::atomic<bool> abort{false};
    std::atomic<bool> ok{true};

    // Decompression pipeline:
    // - producer thread: readFixedFields + Decoder/DecoderByRange (+ optional other-fields)
    // - consumer thread: initalIndex + output (VCF/BCF/BED/PGEN/BGEN)
    //
    // This replaces the fragile 3-way barrier synchronization and avoids deadlocks on early-exit paths.
    const size_t pipeline_depth = kPrefetchDepth;
    BoundedQueue<std::unique_ptr<DecompressedChunk>> free_q(pipeline_depth);
    BoundedQueue<std::unique_ptr<DecompressedChunk>> ready_q(pipeline_depth);

    for (size_t i = 0; i < pipeline_depth; ++i)
        free_q.Push(std::make_unique<DecompressedChunk>());

    std::thread process_thread([&] {
        auto logger = LogManager::Instance().Logger();
        for (uint32_t chunk_id = start_chunk_id; chunk_id < end_chunk_id; ++chunk_id)
        {
            if (abort.load(std::memory_order_relaxed))
                break;

            std::unique_ptr<DecompressedChunk> buf;
            if (!free_q.Pop(buf))
                break;

            buf->chunk_id = chunk_id;
            buf->start_actual_pos = decompression_reader.getActualPos(chunk_id);
            buf->end_actual_pos = decompression_reader.getActualPos(chunk_id + 1);
            buf->chunk_variant_offset = 0;
            buf->fixed_variants_chunk.clear();
            buf->sort_perm.clear();
            buf->gt_indexes.clear();
            buf->row_block_ranges.clear();
            buf->all_fields.clear();

            if (!decompression_reader.readFixedFields())
            {
                logger->error("process_thread: readFixedFields FAILED at chunk_id={}", chunk_id);
                abort.store(true, std::memory_order_relaxed);
                ok.store(false, std::memory_order_relaxed);
                ready_q.Complete();
                free_q.Complete();
                break;
            }

            buf->has_fixed_fields_rb_dir = decompression_reader.has_fixed_fields_rb_dir;
            buf->fixed_fields_chunk_version = decompression_reader.fixed_fields_chunk_version;

            if (params.compress_mode == compress_mode_t::lossless_mode && params.samples.empty() &&
                (out_type == file_type::VCF_File || out_type == file_type::BCF_File))
            {
                if (chunk_id >= decompression_reader.actual_variants.size())
                {
                    logger->error("process_thread: actual_variants out of bounds (chunk_id={}, size={})",
                                  chunk_id, decompression_reader.actual_variants.size());
                    abort.store(true, std::memory_order_relaxed);
                    ok.store(false, std::memory_order_relaxed);
                    ready_q.Complete();
                    free_q.Complete();
                    break;
                }

                uint32_t no_actual_variants = decompression_reader.actual_variants[chunk_id];
                buf->all_fields.reserve(no_actual_variants);
                while (no_actual_variants--)
                {
                    buf->all_fields.emplace_back(vector<field_desc>(decompression_reader.keys.size()));
                    decompression_reader.GetVariants(buf->all_fields.back());
                }
            }

            uint32_t variants_before = 0;
            bool decode_ok = false;
            if (range != "" && buf->has_fixed_fields_rb_dir && !decompression_reader.useLegacyPath)
            {
                decode_ok = decompression_reader.DecoderByRange(buf->fixed_variants_chunk, buf->sort_perm, buf->gt_indexes,
                                                               chunk_id, range_1, range_2, variants_before,
                                                               buf->row_block_ranges);
            }
            else
            {
                decode_ok = decompression_reader.Decoder(buf->fixed_variants_chunk, buf->sort_perm, buf->gt_indexes, chunk_id);
            }

            if (!decode_ok)
            {
                logger->error("process_thread: Decoder FAILED at chunk_id={}", chunk_id);
                abort.store(true, std::memory_order_relaxed);
                ok.store(false, std::memory_order_relaxed);
                ready_q.Complete();
                free_q.Complete();
                break;
            }
            buf->chunk_variant_offset = variants_before;

            if (!ready_q.Push(std::move(buf)))
                break;
        }
        ready_q.Complete();
    });

    std::thread decompress_thread([&] {
        auto logger = LogManager::Instance().Logger();

        std::unique_ptr<DecompressedChunk> buf;
        while (ready_q.Pop(buf))
        {
            if (!buf)
                continue;
            if (abort.load(std::memory_order_relaxed))
                break;

            // Keep legacy semantics: cur_chunk_id is 1-based in the output stage
            // (code uses `cur_chunk_id-1` as the decoded chunk index).
            cur_chunk_id = buf->chunk_id + 1;
            start_chunk_actual_pos = buf->start_actual_pos;
            end_chunk_actual_pos = buf->end_actual_pos;
            chunk_variant_offset_io = buf->chunk_variant_offset;
            has_fixed_fields_rb_dir_io = buf->has_fixed_fields_rb_dir;
            fixed_fields_chunk_version_io = buf->fixed_fields_chunk_version;

            fixed_variants_chunk_io.swap(buf->fixed_variants_chunk);
            sort_perm_io.swap(buf->sort_perm);
            decompress_gt_indexes_io.swap(buf->gt_indexes);
            row_block_ranges_io.swap(buf->row_block_ranges);
            all_fields_io.swap(buf->all_fields);

            bool chunk_ok = true;
            if (!initalIndex())
            {
                logger->error("decompress_thread: initalIndex FAILED at cur_chunk_id={}", cur_chunk_id);
                chunk_ok = false;
            }
            else if (out_type == file_type::BED_File)
            {
                chunk_ok = (BedFormatDecompress() == 0);
            }
            else if (out_type == file_type::PGEN_File)
            {
                chunk_ok = (PgenFormatDecompress() == 0);
            }
            else if (out_type == file_type::BGEN_File)
            {
                chunk_ok = (BgenFormatDecompress() == 0);
            }
            else if (out_type == file_type::GDS_File)
            {
                // Not implemented; OpenForWriting() should have rejected this already.
                chunk_ok = false;
            }
            else
            {
                if (params.compress_mode == compress_mode_t::lossless_mode)
                {
                    if (params.samples.empty())
                    {
                        chunk_ok = (decompressAll() == 0);
                        all_fields_io.clear(); // release other-fields payloads for this chunk
                    }
                    else
                    {
                        // Lossless sample query: decode fixed fields + GT only.
                        chunk_ok = (decompressSampleSmart(range) == 0);
                    }
                }
                else
                {
                    if (params.samples == "")
                        chunk_ok = (decompressRange(range) == 0);
                    else
                        chunk_ok = (decompressSampleSmart(range) == 0);
                }
            }

            // Return buffers back to pool (keep capacity, drop content).
            fixed_variants_chunk_io.clear();
            sort_perm_io.clear();
            decompress_gt_indexes_io.clear();
            row_block_ranges_io.clear();
            all_fields_io.clear();

            buf->fixed_variants_chunk.swap(fixed_variants_chunk_io);
            buf->sort_perm.swap(sort_perm_io);
            buf->gt_indexes.swap(decompress_gt_indexes_io);
            buf->row_block_ranges.swap(row_block_ranges_io);
            buf->all_fields.swap(all_fields_io);
            buf->chunk_variant_offset = 0;
            buf->has_fixed_fields_rb_dir = false;
            buf->fixed_fields_chunk_version = 0;
            buf->start_actual_pos = 0;
            buf->end_actual_pos = 0;

            if (!free_q.Push(std::move(buf)))
                break;

            ClearUniqueCache();

            if (!chunk_ok)
            {
                abort.store(true, std::memory_order_relaxed);
                ok.store(false, std::memory_order_relaxed);
                ready_q.Complete();
                free_q.Complete();
                break;
            }
        }
    });

    process_thread.join();
    decompress_thread.join();
    return ok.load(std::memory_order_relaxed);
}

void Decompressor::ClearUniqueCache()
{
    if (done_unique.empty())
    {
        stored_unique.clear();
        return;
    }
    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();
    stored_unique.clear();
}

void Decompressor::finalizePgen()
{
    auto logger = LogManager::Instance().Logger();
    if (out_file_name == "-")
        return;
    const std::string path = out_file_name + ".pgen";
    if (!patch_u32_le_at(path, 3, export_variant_count_))
        logger->error("Failed to patch PGEN variant_count in {}", path);
}

void Decompressor::finalizeBgen()
{
    auto logger = LogManager::Instance().Logger();
    if (out_file_name == "-")
        return;
    const std::string path = out_file_name + ".bgen";
    if (!patch_u32_le_at(path, 8, export_variant_count_))
        logger->error("Failed to patch BGEN variant_count in {}", path);
}

bool Decompressor::Close(){
		if (in_file_name == "-") {
	        if(remove(decompression_reader.fname.c_str()) != 0)
			    perror("Error deleting temp file");

		}
    if(out_type == file_type::BED_File){
        out_fam.Close();
        out_bed.Close();
        out_bim.Close();
    }
    else if (out_type == file_type::PGEN_File)
    {
        out_pgen.Close();
        out_pvar.Close();
        out_psam.Close();
        finalizePgen();
    }
    else if (out_type == file_type::BGEN_File)
    {
        out_bgen.Close();
        out_sample.Close();
        finalizeBgen();
    }
    else{
        // Finalize parallel writer before closing output file
        if (parallel_writer_) {
            parallel_writer_->Finalize();
            delete parallel_writer_;
            parallel_writer_ = nullptr;
        }

        if(out_hdr){
            bcf_hdr_destroy(out_hdr);
            out_hdr = nullptr;
        }
        if (out_file)
        {
            hts_close(out_file);
            out_file = nullptr;
        }
        if(params.split_flag){
            if(split_files[cur_file]){
                hts_close(split_files[cur_file]);
                split_files[cur_file] = nullptr;
            }
        }

        if(rec){
            bcf_destroy1(rec);
        }
        if(str.m)
            free(str.s);
    }
    if (tmp_arr)
    {
        delete[] tmp_arr;
        tmp_arr = nullptr;
    }
    if (str.m)
    {
        str.s = nullptr;
        str.m = 0;
        str.l = 0;
    }
    return true;
}
// *****************************************************************************************************************
// bool Decompressor::createSamplesNameFile(vector<string> &v_samples)
// {
//     std::ofstream sn_file(out_samples_file_name + ".sn");
//     if (sn_file)
//     {
//         for (uint32_t i = 0; i < v_samples.size(); i++)
//         {
//             sn_file << v_samples[i] << std::endl;
//         }
//         std::cerr << "File with list of samples (" << out_samples_file_name + ".sn"
//                   << ") created." << std::endl;

//         sn_file.close();
//         return true;
//     }
//     else
//     {
//         std::cerr << "Could not open " << out_samples_file_name + ".sn"
//                   << "file with list of samples." << std::endl;
//         return false;
//     }
//     return true;
// }


// *****************************************************************************************************************
void Decompressor::initialXORLut(){

    uint8_t cur_xor_result = 0;
    uint8_t temp = 0;
    for (int i = 0; i < 256; ++i){
        cur_xor_result = i & perm_lut8[0];
        for(int j = 7; j > 0; --j){
            temp = ((i >> j) ^ (i >> (j-1))) & 1;
            cur_xor_result += temp * perm_lut8[8 - j];
        }
        map_t256[cur_xor_result] = i;
    }

}
// *****************************************************************************************************************
// void Decompressor::initialLut1()
// {
    
//     uint8_t mask;
//     for (int c = 0; c < 8; c++)
//     {
//         mask = 0x80 >> c;
//         for (int i = 0; i < 256; i++)
//         {
//             if (i & mask)
//             {
//                 for (int j = 0; j < 256; j++)
//                 {
//                     if (j & mask) // 11
//                     {
//                         gt_lookup_table[i][j][c] = '0';
//                     }
//                     else // 10
//                     {
//                         gt_lookup_table[i][j][c] = '.';
//                     }
//                 }
//             }
//             else
//             {
//                 for (int j = 0; j < 256; j++)
//                 {
//                     if (j & mask) // 01
//                     {
//                         gt_lookup_table[i][j][c] = '1';
//                     }
//                     else // 00
//                     {
//                         gt_lookup_table[i][j][c] = '0';
//                     }
//                 }
//             }
//         }
//     }
// }
// *****************************************************************************************************************
void Decompressor::initialLut()
{
    // xor_map_table();
    uint8_t mask;
    for (int c = 0; c < 8; c++)
    {
        mask = 0x80 >> c;
        for (int i = 0; i < 256; i++)
        {
            if (i & mask)
            {
                for (int j = 0; j < 256; j++)
                {
                    if (j & mask) // 11
                    {
                        gt_lookup_table[i][j][c] = '2';
                    }
                    else // 10
                    {
                        gt_lookup_table[i][j][c] = '.';
                    }
                }
            }
            else
            {
                for (int j = 0; j < 256; j++)
                {
                    if (j & mask) // 01
                    {
                        gt_lookup_table[i][j][c] = '1';
                    }
                    else // 00
                    {
                        gt_lookup_table[i][j][c] = '0';
                    }
                }
            }
        }
    }
}

void Decompressor::initGscToPgenLUT()
{
    for (int p0 = 0; p0 < 256; ++p0)
    {
        for (int p1 = 0; p1 < 256; ++p1)
        {
            uint8_t out_byte = 0;
            for (int sample_idx = 0; sample_idx < 4; ++sample_idx)
            {
                const char h0 = gt_lookup_table[p0][p1][sample_idx * 2];
                const char h1 = gt_lookup_table[p0][p1][sample_idx * 2 + 1];
                const bool missing = (h0 == '.' || h1 == '.' || h0 == '2' || h1 == '2');
                const uint8_t gt = missing ? 0b11 : static_cast<uint8_t>((h0 == '1') + (h1 == '1'));
                out_byte |= static_cast<uint8_t>(gt << (sample_idx * 2));
            }
            gsc_to_pgen_lut[p0][p1] = out_byte;
        }
    }
    gsc_to_pgen_lut_ready_ = true;
}

void Decompressor::initGscToBgenLUT()
{
    for (int p0 = 0; p0 < 256; ++p0)
    {
        for (int p1 = 0; p1 < 256; ++p1)
        {
            for (int sample_idx = 0; sample_idx < 4; ++sample_idx)
            {
                const char h0 = gt_lookup_table[p0][p1][sample_idx * 2];
                const char h1 = gt_lookup_table[p0][p1][sample_idx * 2 + 1];
                const bool missing = (h0 == '.' || h1 == '.' || h0 == '2' || h1 == '2');
                if (missing)
                {
                    gsc_to_bgen_ploidy[p0][p1][sample_idx] = 0x82;
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2] = 0;
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2 + 1] = 0;
                    continue;
                }

                gsc_to_bgen_ploidy[p0][p1][sample_idx] = 0x02;
                const int alt_count = (h0 == '1') + (h1 == '1');
                if (alt_count == 0)
                {
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2] = 255;
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2 + 1] = 0;
                }
                else if (alt_count == 1)
                {
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2] = 0;
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2 + 1] = 255;
                }
                else
                {
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2] = 0;
                    gsc_to_bgen_prob[p0][p1][sample_idx * 2 + 1] = 0;
                }
            }
        }
    }
    gsc_to_bgen_lut_ready_ = true;
}
void Decompressor::getRangeGT(uint8_t *a,const vector<uint32_t> &rev_perm,size_t no_rec, vector<uint8_t> &str){

    for(size_t i = 0;i<no_rec;i++){
        int cur_byte = rev_perm[i] >> 3;
        int cur_pos = rev_perm[i] % 8;
        int r1 = a[cur_byte];
        int r2 = a[decompression_reader.vec_len + cur_byte];
        str[i]= gt_lookup_table[r1][r2][cur_pos];
    }   
    
}
// // *****************************************************************************************************************
void Decompressor::getRangeGT(uint8_t *a, size_t no_rec, string &str)
{
    // int end,g;
    // if(no_rec & 7)//%8)
    // {
    //     end = decompression_reader.vec_len - 1;
    //     g = no_rec & 7;//%8;
    // }
    // else
    // {
    //     end = decompression_reader.vec_len;
    //     g = 0;
    // }
    // long long * lookup_table_ptr = nullptr;
    // long long *tmp_arr = new long long[end];
    // int vec1_start = 0,vec2_start=decompression_reader.vec_len;
    // for (vec1_start = 0; vec1_start < end; ++vec1_start)
    // {
    //     //memcpy(pt + (vec1_start << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]], 8);
    //     lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
    //     //*(tmp_arr+ vec1_start) = *lookup_table_ptr;
    //     tmp_arr[vec1_start] = *lookup_table_ptr;
    // }
    // memcpy(str.data(), tmp_arr, end << 3);
    // if(g)
    // {
    //     memcpy(str.data() + (end << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], g);
        
    // }   
    // if(tmp_arr)
    //     delete [] tmp_arr;  
}
// void Decompressor::getSamplesGT(vector<uint8_t> &str, uint8_t res1, uint8_t res2, int p)
// {

//     str += gt_lookup_table[res1][res2][p];

// }
// *****************************************************************************************************************
inline bool Decompressor::comp(const variant_desc_t &v1, const variant_desc_t &v2)
{
    return v1.pos <= v2.pos;
}
// *****************************************************************************************************************
inline bool Decompressor::comp1(const variant_desc_t &v1, const variant_desc_t &v2)
{
    return v1.pos < v2.pos;
}

// *****************************************************************************************************************
bool Decompressor::initalIndex()
{
    id_pos.clear();
    id_pos.reserve(no_variants_in_buf);

    const uint8_t *data = decompress_gt_indexes_io.data();
    const size_t size = decompress_gt_indexes_io.size();
    size_t cur = 0;
    while (cur < size)
    {
        const void *hit = memchr(data + cur, '\0', size - cur);
        if (!hit)
            break;
        const size_t idx = static_cast<const uint8_t *>(hit) - data;
        if (idx > std::numeric_limits<uint32_t>::max())
        {
            auto logger = LogManager::Instance().Logger();
            logger->error("initalIndex: gt_index offset {} exceeds uint32_t max", idx);
            return false;
        }
        id_pos.emplace_back(static_cast<uint32_t>(idx));
        cur = idx + 1;
    }
    return true;
}

// // *****************************************************************************************************************
void Decompressor::appendVCF(variant_desc_t &_desc, vector<uint8_t> &_my_str, size_t _no_haplotypes)
{

    bcf_clear(rec);
    record.clear();
    record.reserve(_desc.chrom.size() + _desc.id.size() + _desc.ref.size() +
                   _desc.alt.size() + _desc.qual.size() + _desc.filter.size() +
                   _desc.info.size() + 16);
    record.append(_desc.chrom);
    record.append("\t0\t");
    record.append(_desc.id);
    record.push_back('\t');
    record.append(_desc.ref);
    record.push_back('\t');
    record.append(_desc.alt);
    record.push_back('\t');
    record.append(_desc.qual);
    record.push_back('\t');
    record.append(_desc.filter);
    record.push_back('\t');
    record.append(_desc.info);
    kstring_t s;
    s.s = (char *)record.c_str();
    s.m = record.length();
    s.l = 0;
    vcf_parse(&s, out_hdr, rec);
    rec->pos = (int32_t)(_desc.pos - 1);
  
    if (out_genotypes)
    {
        if(count <= INT8_MAX && count > INT8_MIN+1){
            gt_arr_i8.resize(_no_haplotypes);

            int t = 0;
            for (auto cur_genotype : _my_str)
            {
                if (cur_genotype == '.')
                    gt_arr_i8[t++] = bcf_gt_missing;
                else
                    gt_arr_i8[t++] = bcf_gt_phased(int(cur_genotype - '0'));

            }
            
            str.l = 3;
            memcpy(str.s + str.l, gt_arr_i8.data(), _no_haplotypes * sizeof(int8_t));
            str.l += _no_haplotypes;
            str.s[str.l] = 0;
            // for(int i=0;i<str.l;i++)
            //     std::cerr<<(int)str.s[i]<;
    
            // bcf_update_genotypes(out_hdr, rec, gt_arr.data(), _no_haplotypes);
            bcf_update_genotypes_fast(out_hdr, rec,str);
        }
        else{

            gt_arr_i32.resize(_no_haplotypes);

            int t = 0;
            for (auto cur_genotype : _my_str)
            {
                if (cur_genotype == '.')
                    gt_arr_i32[t++] = bcf_gt_missing;
                else
                    gt_arr_i32[t++] = bcf_gt_phased(int(cur_genotype - '0'));

            }
            bcf_update_genotypes(out_hdr, rec, gt_arr_i32.data(), _no_haplotypes);
        }
    }
    if(params.split_flag){

        if(_desc.chrom != cur_chrom){
            if(cur_file != -1 ){
                if (split_files[cur_file])
                {
                    hts_close(split_files[cur_file]);
                    split_files[cur_file] = nullptr;
                }
            }
            cur_chrom = _desc.chrom;
            cur_file++;

        }

        bcf_write1(split_files[cur_file], out_hdr, rec);
    }else {
        // Use parallel writer if enabled
        if (use_parallel_writing_ && parallel_writer_) {
            // Reuse pooled record for parallel writing (writer takes ownership)
            bcf1_t* rec_copy = parallel_writer_->AcquireRecord();
            if (!rec_copy)
                rec_copy = bcf_dup(rec);
            else
                bcf_copy(rec_copy, rec);
            parallel_writer_->WriteRecord(rec_copy);
        } else {
            bcf_write1(out_file, out_hdr, rec);
        }
    }
    // bcf_write(out_file, out_hdr, rec);

}
// // *****************************************************************************************************************

void Decompressor::appendVCFToRec(variant_desc_t &_desc, vector<uint8_t> &_genotype, size_t _standard_block_size, vector<field_desc> &_fields, vector<key_desc> &_keys)
{
    auto logger = LogManager::Instance().Logger();
    logger->debug("appendVCFToRec: entering, _genotype.size()={}, _standard_block_size={}, _fields.size()={}, _keys.size()={}",
                  _genotype.size(), _standard_block_size, _fields.size(), _keys.size());
    if (_genotype.size() < _standard_block_size) {
        logger->error("appendVCFToRec: _genotype.size()={} < _standard_block_size={}", _genotype.size(), _standard_block_size);
    }
    if (_fields.size() != _keys.size()) {
        logger->error("appendVCFToRec: _fields.size()={} != _keys.size()={}", _fields.size(), _keys.size());
    }
    const bool use_parallel_record = (use_parallel_writing_ && parallel_writer_ && !params.split_flag);
    bcf1_t *target_rec = rec;
    if (use_parallel_record)
        target_rec = parallel_writer_->AcquireRecord();
    if (!target_rec)
    {
        logger->error("appendVCFToRec: failed to acquire output record");
        return;
    }

    bcf_clear(target_rec);
    record.clear();
    record.reserve(_desc.chrom.size() + _desc.id.size() + _desc.ref.size() +
                   _desc.alt.size() + _desc.qual.size() + 16);
    record.append(_desc.chrom);
    record.append("\t0\t");
    record.append(_desc.id);
    record.push_back('\t');
    record.append(_desc.ref);
    record.push_back('\t');
    record.append(_desc.alt);
    record.push_back('\t');
    record.append(_desc.qual);
    record.append("\t.\t.");
    kstring_t s;
    s.s = (char *)record.c_str();
    s.m = record.length();
    s.l = 0;
    vcf_parse(&s, out_hdr, target_rec);
    target_rec->pos = (int32_t)(_desc.pos - 1);
    int curr_size = 0;
    for (size_t i = 0; i < _fields.size(); i++)
    {
        if (i >= _keys.size()) {
            logger->error("appendVCFToRec: filter loop i={} >= _keys.size()={}", i, _keys.size());
            break;
        }
        uint32_t id = _keys[i].actual_field_id;
        if (id >= _keys.size() || id >= _fields.size()) {
            logger->error("appendVCFToRec: filter loop actual_field_id={} out of bounds (keys={}, fields={})", id, _keys.size(), _fields.size());
            continue;
        }
        if (_keys[id].keys_type == key_type_t::flt)
        {

            if (_fields[id].present)
                bcf_add_filter(out_hdr, target_rec, _keys[id].key_id);
        }
    }

    for (size_t i = 0; i < _keys.size(); i++)
    {
        uint32_t id = _keys[i].actual_field_id;
        if (id >= _keys.size() || id >= _fields.size()) {
            logger->error("appendVCFToRec: info loop actual_field_id={} out of bounds (keys={}, fields={})", id, _keys.size(), _fields.size());
            continue;
        }
        if (_keys[id].keys_type == key_type_t::info)
        {

            if (_fields[id].present)
                switch (_keys[id].type)
		        {
		        case BCF_HT_INT:
                    curr_size = _fields[id].data_size >> 2;
                    // std::cerr<<curr_size<<endl;
                    bcf_update_info_int32(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size);
			        break;
		        case BCF_HT_REAL:
                    curr_size = _fields[id].data_size >> 2;
                    bcf_update_info_float(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size);
			        break;
		        case BCF_HT_STR:
                    curr_size = _fields[id].data_size;
                    bcf_update_info(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_STR);
			        break;
		        case BCF_HT_FLAG:
                    curr_size = 1;
                    bcf_update_info_flag(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size);
			        break;
		        }

                
        }
    }

    // // FORMAT

    // field_desc gt_phased;
    for (size_t i = 0; i < _keys.size(); i++)
    {
        int id = _keys[i].actual_field_id;
        if (id < 0 || static_cast<size_t>(id) >= _fields.size()) {
            logger->error("appendVCFToRec: FORMAT loop id={} out of bounds (fields.size()={})", id, _fields.size());
            continue;
        }

        if(id == decompression_reader.key_gt_id)
        {
        //     // gt_phased.resize(_fields[id].data_size);
        // //    memcpy(gt_phased.data(), _fields[id].data, _fields[id].data_size);
        //     // gt_phased = move(_fields[id]);
            // Calculate actual sample count from _standard_block_size
            uint32_t actual_samples = _standard_block_size / decompression_reader.ploidy;
            bool is_sample_subset = (actual_samples != decompression_reader.n_samples);

            bool GT_NULL_flag = false;
            for(uint32_t i = 0; i < actual_samples; i++)
                if(_genotype[i] == '.'){
                    GT_NULL_flag = true;
                    break;
                }
            if(count <= INT8_MAX && count > INT8_MIN+1&&!GT_NULL_flag){

                gt_arr_u8.resize(_standard_block_size);
                uint32_t t = 0;
                int cur_gt = 0;
                for(uint32_t i = 0; i < actual_samples; i++){
                    cur_gt = i*decompression_reader.ploidy;
                    if(_genotype[cur_gt] == '.')
                        gt_arr_u8[cur_gt] = bcf_gt_missing;
                    else
                        gt_arr_u8[cur_gt] = bcf_gt_unphased((int)(_genotype[cur_gt]-'0'));
                    for(uint32_t j = 1;j< (uint32_t)decompression_reader.ploidy;j++){
                        cur_gt++;
                        // In sample subset mode, use phase info from the original sample
                        uint32_t orig_sample_idx = is_sample_subset ? sampleIDs[i] : i;
                        uint32_t phase_idx = orig_sample_idx * (decompression_reader.ploidy - 1) + (j - 1);
                        char phase_char = (phase_idx < _fields[id].data_size) ? _fields[id].data[phase_idx] : '|';
                        if(phase_char == '/'){
                            if(_genotype[cur_gt] == '.')
                                    gt_arr_u8[cur_gt] = bcf_gt_missing;
                                else
                                    gt_arr_u8[cur_gt] = bcf_gt_unphased((int)(_genotype[cur_gt]-'0'));

                        }
                        else if(phase_char == '|'){
                            if(_genotype[cur_gt] == '.')
                                gt_arr_u8[cur_gt] = bcf_next_gt_missing;
                            else
                                gt_arr_u8[cur_gt] = bcf_gt_phased((int)(_genotype[cur_gt]-'0'));
                        }
                        if (!is_sample_subset)
                            t++;

                    }
                }

                str.l = 3;

                memcpy(str.s + str.l, gt_arr_u8.data(), _standard_block_size * sizeof(uint8_t));
                str.l += _standard_block_size;
                str.s[str.l] = 0;

                bcf_update_genotypes_fast(out_hdr, target_rec, str);
            }
            else{

                gt_arr_i32.resize(_standard_block_size);
                uint32_t t = 0;
                int cur_gt = 0;
                for(uint32_t i = 0; i < actual_samples; i++){
                    cur_gt = i*decompression_reader.ploidy;
                    if(_genotype[cur_gt] == '.')
                        gt_arr_i32[cur_gt] = bcf_gt_missing;
                    else
                        gt_arr_i32[cur_gt] = bcf_gt_unphased(_genotype[cur_gt]-'0');
                    for(uint32_t j = 1;j< (uint32_t)decompression_reader.ploidy;j++){
                        cur_gt++;
                        // In sample subset mode, use phase info from the original sample
                        uint32_t orig_sample_idx = is_sample_subset ? sampleIDs[i] : i;
                        uint32_t phase_idx = orig_sample_idx * (decompression_reader.ploidy - 1) + (j - 1);
                        char phase_char = (phase_idx < _fields[id].data_size) ? _fields[id].data[phase_idx] : '|';
                        if(phase_char == '/'){
                        if(_genotype[cur_gt] == '.')
                                gt_arr_i32[cur_gt] = bcf_gt_missing;
                            else
                                gt_arr_i32[cur_gt] = bcf_gt_unphased(_genotype[cur_gt]-'0');

                        }
                        else if(phase_char == '|'){
                            if(_genotype[cur_gt] == '.')
                                gt_arr_i32[cur_gt] = bcf_next_gt_missing;
                            else
                                gt_arr_i32[cur_gt] = bcf_gt_phased(_genotype[cur_gt]-'0');
                        }
                        else{
                            gt_arr_i32[cur_gt] =  GT_NOT_CALL;
                        }
                        if (!is_sample_subset)
                            t++;
                    }
                }
                bcf_update_genotypes(out_hdr, target_rec, gt_arr_i32.data(), _standard_block_size);
            }
            continue;
        }
        if (_keys[id].keys_type == key_type_t::fmt)
        {
            // In sample subset mode, skip non-GT FORMAT fields as they contain data for all samples
            uint32_t actual_samples = _standard_block_size / decompression_reader.ploidy;
            bool is_sample_subset = (actual_samples != decompression_reader.n_samples);
            if (is_sample_subset)
                continue;

	            if (_fields[id].present)
	            {
	                switch (_keys[id].type)
			        {
			        case BCF_HT_INT:
	                    curr_size = _fields[id].data_size >> 2;
	                    if (_keys[id].name == "DP" || _keys[id].name == "GQ" || _keys[id].name == "MIN_DP" || _keys[id].name == "RGQ")
	                    {
	                        const int expected = (int)decompression_reader.n_samples;
	                        if (curr_size != expected)
	                        {
	                            auto logger = LogManager::Instance().Logger();
	                            logger->error("FMT {} size mismatch: have {} ints, expected {}", _keys[id].name, curr_size, expected);
	                            abort();
	                        }
	                    }
	                    else if (_keys[id].name == "SB")
	                    {
	                        const int expected = (int)(decompression_reader.n_samples * 4);
	                        if (curr_size != expected)
	                        {
	                            auto logger = LogManager::Instance().Logger();
	                            logger->error("FMT {} size mismatch: have {} ints, expected {}", _keys[id].name, curr_size, expected);
	                            abort();
	                        }
	                    }
	                    bcf_update_format(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_INT);
				        break;
		        case BCF_HT_REAL:
                    curr_size = _fields[id].data_size >> 2;
                    bcf_update_format(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_REAL);
			        break;
		        case BCF_HT_STR:
                    curr_size = _fields[id].data_size;
                    bcf_update_format(out_hdr, target_rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_STR);
			        break;
		        case BCF_HT_FLAG:
                    assert(0);
			        break;
		        }


            }
        }
    }
    if(params.split_flag){

        if(_desc.chrom != cur_chrom){
            if(cur_file != -1 ){
                if (split_files[cur_file])
                {
                    hts_close(split_files[cur_file]);
                    split_files[cur_file] = nullptr;
                }
            }
            cur_chrom = _desc.chrom;
            cur_file++;

        }

        bcf_write1(split_files[cur_file], out_hdr, target_rec);
    }else {
        if (use_parallel_record) {
            parallel_writer_->WriteRecord(target_rec);
        } else {
            bcf_write1(out_file, out_hdr, target_rec);
        }
    }
}
// 
bool Decompressor::SetVariantToRec(variant_desc_t &desc, vector<field_desc> &fields, vector<key_desc> &keys, vector<uint8_t> &_my_str, size_t _standard_block_size)
{
    
    if(desc.alt.find("<M>") == string::npos){

        if(desc.alt.find("<N>") == string::npos){
            if(count){
                
                appendVCFToRec(temp_desc, genotype, _standard_block_size, temp_fields, keys);  
            }
            
            appendVCFToRec(desc, _my_str, _standard_block_size, fields, keys);
            count = 0;
        
        }
        else{
            if(count){
                appendVCFToRec(temp_desc, genotype, _standard_block_size, temp_fields, keys);  
            }
            
            genotype = _my_str;
            // temp_fields.clear();
            // std::cerr<<"start move"<<endl;
            temp_fields.resize(fields.size());
            // temp_fields = std::move(fields);
            for(size_t i = 0;i<fields.size();i++)
                temp_fields[i] = std::move(fields[i]);
            // std::cerr<<"end move"<<endl;
            temp_desc = desc;
            temp_desc.alt = desc.alt.substr(0, desc.alt.find_first_of(','));
            count = 1;

        }
        fields_pos++;
    }
    else{
        if(count == 0){
            appendVCFToRec(desc, _my_str, _standard_block_size, fields, keys);
            count = 0;
            fields_pos++;
        }else{
            temp_desc.alt.append(",");
            temp_desc.alt += desc.alt.substr(0, desc.alt.find_first_of(','));
            uint8_t target = uint8_t('2' + count-1);
            for (size_t i = 0; i < _my_str.size(); i++)
            {
                if (genotype[i] == target && _my_str[i] == '2')
                    genotype[i]++;
            }
            count++;
        }
       
    }
    
    return true;
}
// // // **********************************************************************************************************************************************************

bool Decompressor::SetVariant(variant_desc_t &desc, vector<uint8_t> &_my_str, size_t _standard_block_size)
{
    if(desc.alt.find("<M>") == string::npos){
        if(desc.alt.find("<N>") == string::npos){
            if(count){
                appendVCF(temp_desc, genotype, _standard_block_size);
                
            }
            appendVCF(desc, _my_str, _standard_block_size);
            count = 0;
        }
        else{
            if(count)
                appendVCF(temp_desc, genotype, _standard_block_size);
            
            genotype = _my_str;
            temp_desc = desc;
            temp_desc.alt = desc.alt.substr(0, desc.alt.find_first_of(','));
            count = 1;
        }
    }
    else{
        if(count == 0){
            appendVCF(desc, _my_str, _standard_block_size);
            count = 0;
        }else{
            temp_desc.alt.append(",");
            temp_desc.alt += desc.alt.substr(0, desc.alt.find_first_of(','));
            char target = char('2' + count-1);
            for (size_t i = 0; i < _my_str.size(); i++)
            {

                if (genotype[i] == target && _my_str[i] == '2')
                    genotype[i]++;
            }
            count++;
        }
    }
    
    return true;
}
//*****************************************************************************************************************
int Decompressor::BedFormatDecompress(){

    auto logger = LogManager::Instance().Logger();
    if (!decompression_reader.useLegacyPath)
    {
        return BedFormatDecompressTiled();
    }
    done_unique.clear();
    stored_unique.clear();
    uint32_t cur_block_id = 0;
    uint32_t c_out_line = 0;
    uint32_t no_var = 0;
    uint32_t start_var = 0;
    int vec1_start,vec2_start = decompression_reader.vec_len;
    // fields_pos  = 0;
    vector<uint8_t> my_str(standard_block_size);
    std::string bim_line;
    // vector<uint32_t> rev_perm(standard_block_size);  //2024.1.16注释
    
    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos * 2 - rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) -rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);

    no_var = end_chunk_actual_pos - start_chunk_actual_pos;

    size_t cur_var;

    for (cur_var = start_var; cur_var + standard_block_size <= no_var; cur_var += standard_block_size )
    {
        
        cur_block_id = cur_var / standard_block_size;
        for(size_t i = 0; i < standard_block_size; i++){
            
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);

            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
            
            
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            }   



            // 优化1: 避免重复搜索
            // size_t comma_pos = desc.alt.find(",");
            // std::string desc_alt = (comma_pos == std::string::npos) ? desc.alt : desc.alt.substr(0, comma_pos);


            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[i];
            build_bim_line(bim_line, desc, genetic_distance);
            out_bim.Write(bim_line.c_str(), bim_line.size());
            for(int i = 0; i < (int)standard_block_size; i += 2){
                char first = my_str[i];
                char second = my_str[i + 1];
                out_bed.PutBit(bits_lut[first][second][0]);
                out_bed.PutBit(bits_lut[first][second][1]);
                // if(first == '0' &&  second == '0'){
                    
                //     out_bed.PutBit(1);
                //     out_bed.PutBit(1);
                // }
                // else if((first == '0' &&  second == '1') || (first == '1' &&  second == '0')){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(1);
                // }
                // else if(first == '1' &&  second == '1'){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(0);
                // }
                // else if(first == '.' ||  second == '.'){
                    
                //     out_bed.PutBit(1);
                    
                //     out_bed.PutBit(0);
                // }
            } 
            out_bed.FlushPartialByteBuffer();
          
            
        
        }
    }
    if(no_var % standard_block_size)
    {
        cur_block_id = cur_var / standard_block_size;            
        // reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size);  //2024.1.16
        
        for(;cur_var < no_var;++cur_var){
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
            // if(!(no_var%standard_block_size))
            //     decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
            // else
            //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);

            // for(int i = 0;i<rev_perm.size();i++)
            //     std::cerr<<rev_perm[i]<<" ";
            // std::cerr<<endl;
            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);//2024.1.16 rev_perm ---> sort_perm_io[cur_block_id][0]
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            } 
                
            c_out_line = cur_var % standard_block_size;

            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            build_bim_line(bim_line, desc, genetic_distance);
            out_bim.Write(bim_line.c_str(), bim_line.size());
            for(int i = 0; i < (int)standard_block_size; i += 2){
                char first = my_str[i];
                char second = my_str[i + 1];
                out_bed.PutBit(bits_lut[first][second][0]);
                out_bed.PutBit(bits_lut[first][second][1]);
                // if(first == '0' &&  second == '0'){
                //     out_bed.PutBit(1);
                //     out_bed.PutBit(1);
                // }
                // else if(first == '0' &&  second == '1' || first == '1' &&  second == '0'){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(1);
                // }
                // else if(first == '1' &&  second == '1'){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(0);
                // }
                // else if(first == '.' ||  second == '.'){
                //     out_bed.PutBit(1);
                //     out_bed.PutBit(0);
                // }
            } 
            out_bed.FlushPartialByteBuffer();
            // SetVariantToRec(desc, all_fields_io[fields_pos], decompression_reader.keys, my_str, standard_block_size);
        }
    }
    // if(cur_chunk_id == end_chunk_id && count){
      
    //     appendVCFToRec(temp_desc, genotype, static_cast<uint32_t>(standard_block_size), temp_fields, decompression_reader.keys);
    // }

    logger->info("Processed chunk {}", cur_chunk_id);

    for (auto &it : done_unique)
        delete[] it.second;

    done_unique.clear();

    return 0;

}

//*****************************************************************************************************************
int Decompressor::BedFormatDecompressTiled()
{
    auto logger = LogManager::Instance().Logger();
    logger->debug("BedFormatDecompressTiled: entering");
    done_unique.clear();
    stored_unique.clear();

    // Disable vector caching for tiled mode since different column blocks
    // have different vector lengths, which causes memory issues when reusing cached vectors
    uint64_t saved_max_stored_unique = max_stored_unique;
    max_stored_unique = 0;

    const uint32_t n_col_blocks = decompression_reader.n_col_blocks;
    const uint32_t full_vec_len = decompression_reader.vec_len;
    const uint32_t row_block_variants = row_block_size ? row_block_size : haplotype_count;

    logger->debug("BedFormatDecompressTiled: n_col_blocks={}, full_vec_len={}, row_block_variants={}, haplotype_count={}",
                 n_col_blocks, full_vec_len, row_block_variants, haplotype_count);

    if (!n_col_blocks || !row_block_variants)
    {
        logger->error("Invalid tiling metadata (n_col_blocks={}, row_block_variants={})", n_col_blocks, row_block_variants);
        return 1;
    }

    if (decompression_reader.col_block_ranges.size() != n_col_blocks ||
        decompression_reader.col_block_vec_lens.size() != n_col_blocks)
    {
        logger->error("Mismatch in column block metadata sizes");
        return 1;
    }

    vector<uint8_t> my_str(haplotype_count);
    std::string bim_line;

    uint64_t start_pair_id = static_cast<uint64_t>(start_chunk_actual_pos) * n_col_blocks;
    uint64_t gt_index_base_pair_id = start_pair_id;
    if (has_fixed_fields_rb_dir_io &&
        fixed_fields_chunk_version_io >= GSC_FIXED_FIELDS_RB_VERSION_V2 &&
        chunk_variant_offset_io)
    {
        gt_index_base_pair_id += static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    }
    uint64_t curr_non_copy_vec_id_offset = gt_index_base_pair_id * 2 - rrr_rank_zeros_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_zeros_bit_vector[1](gt_index_base_pair_id) - rrr_rank_copy_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_copy_bit_vector[1](gt_index_base_pair_id);

    uint64_t max_col_vec_len = 0;
    for (auto len : decompression_reader.col_block_vec_lens)
        if (len > max_col_vec_len)
            max_col_vec_len = len;

    if (max_col_vec_len == 0)
    {
        logger->error("max_col_vec_len is 0!");
        return 1;
    }

    vector<uint8_t> col_decomp_perm(max_col_vec_len * 2);
    vector<uint8_t> col_decomp(max_col_vec_len * 2);

    uint64_t pair_base = start_pair_id + static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;

    for (size_t block_id = 0; block_id < fixed_variants_chunk_io.size(); ++block_id)
    {
        const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());

        if (block_id >= sort_perm_io.size())
        {
            logger->error("block_id {} >= sort_perm_io.size() {}", block_id, sort_perm_io.size());
            return 1;
        }
        if (sort_perm_io[block_id].size() != n_col_blocks)
        {
            logger->error("sort_perm_io[{}].size()={} != n_col_blocks={}", block_id, sort_perm_io[block_id].size(), n_col_blocks);
            return 1;
        }

        logger->debug("BedFormatDecompressTiled: processing block_id={}, block_variants={}", block_id, block_variants);

        for (uint32_t var_in_block = 0; var_in_block < block_variants; ++var_in_block)
        {
            variant_desc_t desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];

            // Skip multi-allelic continuation records for BED output
            if (desc.alt.find("<M>") != string::npos)
                continue;

            // Decode GT data for all column blocks
            fill_n(decomp_data, full_vec_len * 2, 0);
            for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
            {
                const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
                const uint32_t col_vec_len = static_cast<uint32_t>(decompression_reader.col_block_vec_lens[cb]);
                const uint64_t pair_index = pair_base + static_cast<uint64_t>(cb) * block_variants + var_in_block;

                fill_n(col_decomp_perm.data(), col_vec_len * 2, 0);
                decoded_vector_row(pair_index * 2, 0, curr_non_copy_vec_id_offset, col_vec_len, 0, col_decomp_perm.data());
                decoded_vector_row(pair_index * 2 + 1, 0, curr_non_copy_vec_id_offset, col_vec_len, col_vec_len, col_decomp_perm.data());

                fill_n(col_decomp.data(), col_vec_len * 2, 0);
                decode_perm_rev(col_vec_len, col_vec_len, col_block_size, sort_perm_io[block_id][cb],
                                col_decomp_perm.data(), col_decomp.data());

                const uint32_t start_hap = decompression_reader.col_block_ranges[cb].first;
                insert_block_bits(decomp_data, col_decomp.data(), full_vec_len, start_hap, col_block_size);
                insert_block_bits(decomp_data + full_vec_len, col_decomp.data() + col_vec_len, full_vec_len, start_hap, col_block_size);
            }

            // Convert to character genotypes
            int vec1_start;
            int vec2_start = full_vec_len;
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if (trailing_bits)
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);

            // Write BIM line
            build_bim_line(bim_line, desc, genetic_distance);
            out_bim.Write(bim_line.c_str(), bim_line.size());

            // Write BED data
            for (int i = 0; i < (int)haplotype_count; i += 2)
            {
                char first = my_str[i];
                char second = my_str[i + 1];
                out_bed.PutBit(bits_lut[first][second][0]);
                out_bed.PutBit(bits_lut[first][second][1]);
            }
            out_bed.FlushPartialByteBuffer();
        }
        pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
    }

    logger->info("Processed chunk {}", cur_chunk_id);

    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();

    // Restore max_stored_unique
    max_stored_unique = saved_max_stored_unique;

    return 0;
}

int Decompressor::PgenFormatDecompress()
{
    auto logger = LogManager::Instance().Logger();
    if (!params.samples.empty())
    {
        logger->error("PGEN output does not currently support --samples subset.");
        return 1;
    }
    if (decompression_reader.useLegacyPath)
    {
        logger->error("PGEN output requires tiled .gsc input (recompress with non-legacy block settings).");
        return 1;
    }
    return PgenFormatDecompressTiled();
}

int Decompressor::PgenFormatDecompressTiled()
{
    auto logger = LogManager::Instance().Logger();
    done_unique.clear();
    stored_unique.clear();

    // Disable vector caching for tiled mode (matches BED tiled implementation)
    const uint64_t saved_max_stored_unique = max_stored_unique;
    max_stored_unique = 0;

    const uint32_t n_col_blocks = decompression_reader.n_col_blocks;
    const uint32_t full_vec_len = decompression_reader.vec_len;
    const uint32_t n_samples = decompression_reader.n_samples;
    const uint32_t row_block_variants = row_block_size ? row_block_size : haplotype_count;

    if (!n_col_blocks || !row_block_variants)
    {
        logger->error("Invalid tiling metadata for PGEN output");
        max_stored_unique = saved_max_stored_unique;
        return 1;
    }
    if (decompression_reader.col_block_ranges.size() != n_col_blocks ||
        decompression_reader.col_block_vec_lens.size() != n_col_blocks)
    {
        logger->error("Mismatch in column block metadata sizes for PGEN output");
        max_stored_unique = saved_max_stored_unique;
        return 1;
    }

    uint64_t start_pair_id = static_cast<uint64_t>(start_chunk_actual_pos) * n_col_blocks;
    uint64_t gt_index_base_pair_id = start_pair_id;
    if (has_fixed_fields_rb_dir_io &&
        fixed_fields_chunk_version_io >= GSC_FIXED_FIELDS_RB_VERSION_V2 &&
        chunk_variant_offset_io)
    {
        gt_index_base_pair_id += static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    }
    uint64_t curr_non_copy_vec_id_offset = gt_index_base_pair_id * 2 - rrr_rank_zeros_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_zeros_bit_vector[1](gt_index_base_pair_id) - rrr_rank_copy_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_copy_bit_vector[1](gt_index_base_pair_id);

    uint64_t max_col_vec_len = 0;
    for (auto len : decompression_reader.col_block_vec_lens)
        if (len > max_col_vec_len)
            max_col_vec_len = len;
    if (max_col_vec_len == 0)
    {
        logger->error("max_col_vec_len is 0 for PGEN output");
        max_stored_unique = saved_max_stored_unique;
        return 1;
    }
    vector<uint8_t> col_decomp_perm(max_col_vec_len * 2);
    vector<uint8_t> col_decomp(max_col_vec_len * 2);
    std::string pvar_line;

    uint64_t pair_base = start_pair_id + static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;

    for (size_t block_id = 0; block_id < fixed_variants_chunk_io.size(); ++block_id)
    {
        const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());
        if (block_id >= sort_perm_io.size() || sort_perm_io[block_id].size() != n_col_blocks)
        {
            logger->error("sort_perm_io mismatch for PGEN output (block_id={}, sort_perm_io.size()={}, expected n_col_blocks={})",
                          block_id, sort_perm_io.size(), n_col_blocks);
            max_stored_unique = saved_max_stored_unique;
            return 1;
        }
        for (uint32_t var_in_block = 0; var_in_block < block_variants; ++var_in_block)
        {
            const variant_desc_t &desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];
            if (!is_biallelic_desc_for_export(desc))
                continue;

            fill_n(decomp_data, full_vec_len * 2, 0);
            for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
            {
                const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
                const uint32_t col_vec_len = static_cast<uint32_t>(decompression_reader.col_block_vec_lens[cb]);
                const uint64_t pair_index = pair_base + static_cast<uint64_t>(cb) * block_variants + var_in_block;

                fill_n(col_decomp_perm.data(), col_vec_len * 2, 0);
                decoded_vector_row(pair_index * 2, 0, curr_non_copy_vec_id_offset, col_vec_len, 0, col_decomp_perm.data());
                decoded_vector_row(pair_index * 2 + 1, 0, curr_non_copy_vec_id_offset, col_vec_len, col_vec_len, col_decomp_perm.data());

                fill_n(col_decomp.data(), col_vec_len * 2, 0);
                decode_perm_rev(col_vec_len, col_vec_len, col_block_size, sort_perm_io[block_id][cb],
                                col_decomp_perm.data(), col_decomp.data());

                const uint32_t start_hap = decompression_reader.col_block_ranges[cb].first;
                insert_block_bits(decomp_data, col_decomp.data(), full_vec_len, start_hap, col_block_size);
                insert_block_bits(decomp_data + full_vec_len, col_decomp.data() + col_vec_len, full_vec_len, start_hap, col_block_size);
            }

            // Write PVAR row
            build_pvar_line(pvar_line, desc);
            out_pvar.Write(pvar_line.c_str(), pvar_line.size());

            // Write packed hardcalls to PGEN (4 samples per byte)
            const uint8_t* bvec0 = decomp_data;
            const uint8_t* bvec1 = decomp_data + full_vec_len;
            const uint32_t full_bytes = n_samples / 4;
            const uint32_t rem = n_samples % 4;
            for (uint32_t i = 0; i < full_bytes; ++i)
                out_pgen.PutByte(gsc_to_pgen_lut[bvec0[i]][bvec1[i]]);
            if (rem)
            {
                uint8_t out_byte = gsc_to_pgen_lut[bvec0[full_bytes]][bvec1[full_bytes]];
                out_byte &= static_cast<uint8_t>(0xFFu >> ((4 - rem) * 2));
                out_pgen.PutByte(out_byte);
            }

            ++export_variant_count_;
        }
        pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
    }

    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();

    max_stored_unique = saved_max_stored_unique;
    return 0;
}

int Decompressor::BgenFormatDecompress()
{
    auto logger = LogManager::Instance().Logger();
    if (!params.samples.empty())
    {
        logger->error("BGEN output does not currently support --samples subset.");
        return 1;
    }
    if (decompression_reader.useLegacyPath)
    {
        logger->error("BGEN output requires tiled .gsc input (recompress with non-legacy block settings).");
        return 1;
    }
    return BgenFormatDecompressTiled();
}

int Decompressor::BgenFormatDecompressTiled()
{
    auto logger = LogManager::Instance().Logger();
    done_unique.clear();
    stored_unique.clear();

    const uint64_t saved_max_stored_unique = max_stored_unique;
    max_stored_unique = 0;

    const uint32_t n_col_blocks = decompression_reader.n_col_blocks;
    const uint32_t full_vec_len = decompression_reader.vec_len;
    const uint32_t n_samples = decompression_reader.n_samples;
    const uint32_t row_block_variants = row_block_size ? row_block_size : haplotype_count;

    if (!n_col_blocks || !row_block_variants)
    {
        logger->error("Invalid tiling metadata for BGEN output");
        max_stored_unique = saved_max_stored_unique;
        return 1;
    }
    if (decompression_reader.col_block_ranges.size() != n_col_blocks ||
        decompression_reader.col_block_vec_lens.size() != n_col_blocks)
    {
        logger->error("Mismatch in column block metadata sizes for BGEN output");
        max_stored_unique = saved_max_stored_unique;
        return 1;
    }

    uint64_t start_pair_id = static_cast<uint64_t>(start_chunk_actual_pos) * n_col_blocks;
    uint64_t gt_index_base_pair_id = start_pair_id;
    if (has_fixed_fields_rb_dir_io &&
        fixed_fields_chunk_version_io >= GSC_FIXED_FIELDS_RB_VERSION_V2 &&
        chunk_variant_offset_io)
    {
        gt_index_base_pair_id += static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    }
    uint64_t curr_non_copy_vec_id_offset = gt_index_base_pair_id * 2 - rrr_rank_zeros_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_zeros_bit_vector[1](gt_index_base_pair_id) - rrr_rank_copy_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_copy_bit_vector[1](gt_index_base_pair_id);

    uint64_t max_col_vec_len = 0;
    for (auto len : decompression_reader.col_block_vec_lens)
        if (len > max_col_vec_len)
            max_col_vec_len = len;
    if (max_col_vec_len == 0)
    {
        logger->error("max_col_vec_len is 0 for BGEN output");
        max_stored_unique = saved_max_stored_unique;
        return 1;
    }
    vector<uint8_t> col_decomp_perm(max_col_vec_len * 2);
    vector<uint8_t> col_decomp(max_col_vec_len * 2);

    const uint32_t uncompressed_size = 4 + n_samples + 2 + 2 * n_samples; // N + ploidy + phased/B + probs
    vector<uint8_t> uncompressed;
    uncompressed.resize(uncompressed_size);

    uint64_t pair_base = start_pair_id + static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;

    for (size_t block_id = 0; block_id < fixed_variants_chunk_io.size(); ++block_id)
    {
        const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());
        if (block_id >= sort_perm_io.size() || sort_perm_io[block_id].size() != n_col_blocks)
        {
            logger->error("sort_perm_io mismatch for BGEN output (block_id={}, sort_perm_io.size()={}, expected n_col_blocks={})",
                          block_id, sort_perm_io.size(), n_col_blocks);
            max_stored_unique = saved_max_stored_unique;
            return 1;
        }
        for (uint32_t var_in_block = 0; var_in_block < block_variants; ++var_in_block)
        {
            const variant_desc_t &desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];
            if (!is_biallelic_desc_for_export(desc))
                continue;

            fill_n(decomp_data, full_vec_len * 2, 0);
            for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
            {
                const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
                const uint32_t col_vec_len = static_cast<uint32_t>(decompression_reader.col_block_vec_lens[cb]);
                const uint64_t pair_index = pair_base + static_cast<uint64_t>(cb) * block_variants + var_in_block;

                fill_n(col_decomp_perm.data(), col_vec_len * 2, 0);
                decoded_vector_row(pair_index * 2, 0, curr_non_copy_vec_id_offset, col_vec_len, 0, col_decomp_perm.data());
                decoded_vector_row(pair_index * 2 + 1, 0, curr_non_copy_vec_id_offset, col_vec_len, col_vec_len, col_decomp_perm.data());

                fill_n(col_decomp.data(), col_vec_len * 2, 0);
                decode_perm_rev(col_vec_len, col_vec_len, col_block_size, sort_perm_io[block_id][cb],
                                col_decomp_perm.data(), col_decomp.data());

                const uint32_t start_hap = decompression_reader.col_block_ranges[cb].first;
                insert_block_bits(decomp_data, col_decomp.data(), full_vec_len, start_hap, col_block_size);
                insert_block_bits(decomp_data + full_vec_len, col_decomp.data() + col_vec_len, full_vec_len, start_hap, col_block_size);
            }

            // Identifying data
            const std::string varid = desc.id.empty() ? "." : desc.id;
            const std::string rsid = varid;
            const std::string chrom = desc.chrom;
            const uint32_t pos = desc.pos;
            const std::string ref = desc.ref;
            const std::string alt = desc.alt;

            out_bgen.WriteUInt(static_cast<uint16_t>(varid.size()), 2);
            out_bgen.Write(varid.c_str(), varid.size());
            out_bgen.WriteUInt(static_cast<uint16_t>(rsid.size()), 2);
            out_bgen.Write(rsid.c_str(), rsid.size());
            out_bgen.WriteUInt(static_cast<uint16_t>(chrom.size()), 2);
            out_bgen.Write(chrom.c_str(), chrom.size());
            out_bgen.WriteUInt(pos, 4);

            out_bgen.WriteUInt(static_cast<uint16_t>(2), 2); // allele_count
            out_bgen.WriteUInt(static_cast<uint32_t>(ref.size()), 4);
            out_bgen.Write(ref.c_str(), ref.size());
            out_bgen.WriteUInt(static_cast<uint32_t>(alt.size()), 4);
            out_bgen.Write(alt.c_str(), alt.size());

            // Uncompressed genotype data (Layout 2)
            uncompressed[0] = static_cast<uint8_t>(n_samples & 0xFFu);
            uncompressed[1] = static_cast<uint8_t>((n_samples >> 8) & 0xFFu);
            uncompressed[2] = static_cast<uint8_t>((n_samples >> 16) & 0xFFu);
            uncompressed[3] = static_cast<uint8_t>((n_samples >> 24) & 0xFFu);

            const size_t ploidy_off = 4;
            const size_t phased_off = ploidy_off + n_samples;
            const size_t b_off = phased_off + 1;
            const size_t prob_off = b_off + 1;

            uncompressed[phased_off] = 0;
            uncompressed[b_off] = 8;

            const uint8_t* bvec0 = decomp_data;
            const uint8_t* bvec1 = decomp_data + full_vec_len;
            const uint32_t full_bytes = n_samples / 4;
            const uint32_t rem = n_samples % 4;

            for (uint32_t i = 0; i < full_bytes; ++i)
            {
                const uint8_t p0 = bvec0[i];
                const uint8_t p1 = bvec1[i];
                memcpy(uncompressed.data() + ploidy_off + i * 4, gsc_to_bgen_ploidy[p0][p1], 4);
                memcpy(uncompressed.data() + prob_off + i * 8, gsc_to_bgen_prob[p0][p1], 8);
            }
            if (rem)
            {
                const uint8_t p0 = bvec0[full_bytes];
                const uint8_t p1 = bvec1[full_bytes];
                memcpy(uncompressed.data() + ploidy_off + full_bytes * 4, gsc_to_bgen_ploidy[p0][p1], rem);
                memcpy(uncompressed.data() + prob_off + full_bytes * 8, gsc_to_bgen_prob[p0][p1], rem * 2);
            }

            // Compress with zlib
            vector<uint8_t> compressed;
            compressed.resize(compressBound(uncompressed.size()));
            uLongf compressed_size = compressed.size();
            const int zres = compress2(compressed.data(), &compressed_size,
                                       uncompressed.data(), uncompressed.size(),
                                       Z_BEST_SPEED);
            if (zres != Z_OK)
            {
                logger->error("BGEN zlib compression failed (zres={})", zres);
                max_stored_unique = saved_max_stored_unique;
                return 1;
            }

            const uint32_t comp_len = static_cast<uint32_t>(compressed_size) + 4;
            const uint32_t uncomp_len = static_cast<uint32_t>(uncompressed.size());
            out_bgen.WriteUInt(comp_len, 4);
            out_bgen.WriteUInt(uncomp_len, 4);
            out_bgen.Write(compressed.data(), static_cast<size_t>(compressed_size));

            ++export_variant_count_;
        }
        pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
    }

    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();

    max_stored_unique = saved_max_stored_unique;
    return 0;
}

//*****************************************************************************************************************
void Decompressor::insert_block_bits(uint8_t *dest, const uint8_t *src, uint32_t dest_vec_len, uint32_t start_hap, uint32_t block_hap)
{
    if (!block_hap)
        return;
    uint32_t src_vec_len = (block_hap + 7) / 8;
    uint32_t byte_offset = start_hap / 8;
    uint32_t bit_offset = start_hap % 8;
    if (bit_offset == 0)
    {
        memcpy(dest + byte_offset, src, src_vec_len);
        return;
    }
    for (uint32_t bit = 0; bit < block_hap; ++bit)
    {
        uint32_t src_byte = bit / 8;
        uint32_t dest_bit = start_hap + bit;
        if (src[src_byte] & perm_lut8[bit % 8])
            dest[dest_bit / 8] |= perm_lut8[dest_bit % 8];
    }
}

//*****************************************************************************************************************
int Decompressor::decompressAllTiled()
{
    auto logger = LogManager::Instance().Logger();
    logger->debug("decompressAllTiled: entering");
    done_unique.clear();
    stored_unique.clear();
    fields_pos = 0;

    const uint32_t n_col_blocks = decompression_reader.n_col_blocks;
    const uint32_t full_vec_len = decompression_reader.vec_len;
    const uint32_t row_block_variants = row_block_size ? row_block_size : haplotype_count;

    logger->debug("decompressAllTiled: n_col_blocks={}, full_vec_len={}, row_block_variants={}, haplotype_count={}",
                 n_col_blocks, full_vec_len, row_block_variants, haplotype_count);

    if (!n_col_blocks || !row_block_variants)
    {
        logger->error("Invalid tiling metadata (n_col_blocks={}, row_block_variants={})", n_col_blocks, row_block_variants);
        return 1;
    }

    logger->debug("decompressAllTiled: col_block_ranges.size()={}, col_block_vec_lens.size()={}",
                 decompression_reader.col_block_ranges.size(), decompression_reader.col_block_vec_lens.size());

    if (decompression_reader.col_block_ranges.size() != n_col_blocks ||
        decompression_reader.col_block_vec_lens.size() != n_col_blocks)
    {
        logger->error("Mismatch in column block metadata sizes");
        return 1;
    }

    vector<uint8_t> my_str(haplotype_count);
    logger->debug("decompressAllTiled: my_str allocated, size={}", my_str.size());

    uint64_t start_pair_id = static_cast<uint64_t>(start_chunk_actual_pos) * n_col_blocks;
    uint64_t gt_index_base_pair_id = start_pair_id;
    if (has_fixed_fields_rb_dir_io &&
        fixed_fields_chunk_version_io >= GSC_FIXED_FIELDS_RB_VERSION_V2 &&
        chunk_variant_offset_io)
    {
        gt_index_base_pair_id += static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    }
    uint64_t curr_non_copy_vec_id_offset = gt_index_base_pair_id * 2 - rrr_rank_zeros_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_zeros_bit_vector[1](gt_index_base_pair_id) - rrr_rank_copy_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_copy_bit_vector[1](gt_index_base_pair_id);

    logger->debug("decompressAllTiled: start_pair_id={}, curr_non_copy_vec_id_offset={}", start_pair_id, curr_non_copy_vec_id_offset);

    uint64_t max_col_vec_len = 0;
    for (auto len : decompression_reader.col_block_vec_lens)
        if (len > max_col_vec_len)
            max_col_vec_len = len;

    logger->debug("decompressAllTiled: max_col_vec_len={}", max_col_vec_len);
    if (max_col_vec_len == 0)
    {
        logger->error("max_col_vec_len is 0!");
        return 1;
    }

    vector<uint8_t> col_decomp_perm(max_col_vec_len * 2);
    vector<uint8_t> col_decomp(max_col_vec_len * 2);

    const uint32_t actual_variants = static_cast<uint32_t>(all_fields_io.size());
    uint32_t total_variants = 0;
    for (const auto &block : fixed_variants_chunk_io)
    {
        total_variants += static_cast<uint32_t>(block.data_compress.size());
    }

    if (!total_variants || !actual_variants)
    {
        logger->error("Invalid variant counts (total_variants={}, actual_variants={})", total_variants, actual_variants);
        return 1;
    }

    logger->debug("decompressAllTiled: fixed_variants_chunk_io.size()={}, sort_perm_io.size()={}, total_variants={}, actual_variants={}",
                 fixed_variants_chunk_io.size(), sort_perm_io.size(), total_variants, actual_variants);

    // Calculate range filtering boundaries if range is specified
    uint32_t no_var = total_variants;
    uint32_t start_var = 0;
    bool has_range = (range != "");

    if (has_range)
    {
        int start_block, end_block;
        int start_position, end_position;

        if (cur_chunk_id - 1 == start_chunk_id)
        {
            calculate_start_position(start_block, start_position);
            start_var = chunk_variant_offset_io + uint32_t(start_block * row_block_variants + start_position);
        }
        if (cur_chunk_id == end_chunk_id)
        {
            calculate_end_position(end_block, end_position);
            no_var = chunk_variant_offset_io + uint32_t(end_block * row_block_variants + end_position);
        }
        logger->debug("decompressAllTiled: range filtering - start_var={}, no_var={}", start_var, no_var);
    }

    // Check if sample subset is specified
    bool has_sample_subset = (params.samples != "");
    uint32_t output_haplotypes = haplotype_count;
    uint32_t output_samples = decompression_reader.n_samples;
    vector<uint8_t> subset_str;
    vector<vector<pair<uint32_t, uint32_t>>> selected_haps_by_cb;
    if (has_sample_subset)
    {
        output_samples = smpl.no_samples;
        output_haplotypes = output_samples * decompression_reader.ploidy;
        subset_str.resize(output_haplotypes);
        logger->debug("decompressAllTiled: sample subset - output_samples={}, output_haplotypes={}", output_samples, output_haplotypes);

        selected_haps_by_cb.resize(n_col_blocks);
        uint32_t assigned = 0;
        for (uint32_t g = 0; g < output_samples; ++g)
        {
            for (uint32_t p = 0; p < decompression_reader.ploidy; ++p)
            {
                uint32_t out_hap = g * decompression_reader.ploidy + p;
                uint32_t orig_hap = sampleIDs[g] * decompression_reader.ploidy + p;
                int found_cb = -1;
                for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
                {
                    uint32_t start = decompression_reader.col_block_ranges[cb].first;
                    uint32_t size = decompression_reader.col_block_ranges[cb].second;
                    if (orig_hap >= start && orig_hap < start + size)
                    {
                        found_cb = static_cast<int>(cb);
                        uint32_t local_orig = orig_hap - start;
                        selected_haps_by_cb[cb].push_back(make_pair(out_hap, local_orig));
                        assigned++;
                        break;
                    }
                }
                if (found_cb < 0)
                {
                    logger->error("decompressAllTiled: failed to map orig_hap={} to any column block", orig_hap);
                    return 1;
                }
            }
        }
        if (assigned != output_haplotypes)
        {
            logger->error("decompressAllTiled: assigned_haplotypes={} != output_haplotypes={}", assigned, output_haplotypes);
            return 1;
        }
    }

    uint64_t pair_base = start_pair_id + static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    size_t processed_variants = 0;
    uint32_t global_var_idx = chunk_variant_offset_io;

    for (size_t block_id = 0; block_id < fixed_variants_chunk_io.size(); ++block_id)
    {
        const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());

        // Validate sort_perm_io access
        if (block_id >= sort_perm_io.size())
        {
            logger->error("block_id {} >= sort_perm_io.size() {}", block_id, sort_perm_io.size());
            return 1;
        }
        if (sort_perm_io[block_id].size() != n_col_blocks)
        {
            logger->error("sort_perm_io[{}].size()={} != n_col_blocks={}", block_id, sort_perm_io[block_id].size(), n_col_blocks);
            return 1;
        }

        struct ColSubsetCtx
        {
            uint32_t cb = 0;
            uint32_t local_count = 0;
            vector<uint32_t> out_hap;
            vector<uint8_t> perm_bit;
            vector<pair<uint32_t, uint32_t>> byte_where; // (byte_no, local_index)
            vector<ByteGroup> groups;
            vector<uint8_t> bytes; // [local_count*2] for parity0 + parity1
        };

        vector<ColSubsetCtx> subset_ctxs;
        if (has_sample_subset)
        {
            subset_ctxs.reserve(n_col_blocks);
            for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
            {
                if (selected_haps_by_cb[cb].empty())
                    continue;

                const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
                vector<uint32_t> orig_to_perm(col_block_size);
                reverse_perm(sort_perm_io[block_id][cb], orig_to_perm, static_cast<int>(col_block_size));

                ColSubsetCtx ctx;
                ctx.cb = cb;
                ctx.local_count = static_cast<uint32_t>(selected_haps_by_cb[cb].size());
                ctx.out_hap.resize(ctx.local_count);
                ctx.perm_bit.resize(ctx.local_count);
                ctx.byte_where.reserve(ctx.local_count);
                ctx.bytes.resize(static_cast<size_t>(ctx.local_count) * 2);

                for (uint32_t i = 0; i < ctx.local_count; ++i)
                {
                    uint32_t out_h = selected_haps_by_cb[cb][i].first;
                    uint32_t local_orig = selected_haps_by_cb[cb][i].second;
                    uint32_t perm_idx = orig_to_perm[local_orig];
                    ctx.out_hap[i] = out_h;
                    ctx.perm_bit[i] = static_cast<uint8_t>(perm_idx & 7);
                    ctx.byte_where.emplace_back(static_cast<uint32_t>(perm_idx >> 3), i);
                }

                sort(ctx.byte_where.begin(), ctx.byte_where.end(),
                     [](const auto &a, const auto &b) { return a.first < b.first || (a.first == b.first && a.second < b.second); });

                for (uint32_t i = 0; i < ctx.local_count;)
                {
                    uint32_t byte_no = ctx.byte_where[i].first;
                    uint32_t j = i + 1;
                    while (j < ctx.local_count && ctx.byte_where[j].first == byte_no)
                        ++j;
                    ctx.groups.push_back(ByteGroup{byte_no, i, j});
                    i = j;
                }

                subset_ctxs.emplace_back(std::move(ctx));
            }
        }

        logger->debug("decompressAllTiled: processing block_id={}, block_variants={}, processed_variants={} (subset_ctxs={})",
                     block_id, block_variants, processed_variants, subset_ctxs.size());

        for (uint32_t var_in_block = 0; var_in_block < block_variants; ++var_in_block, ++processed_variants, ++global_var_idx)
        {
            // Skip variants outside the range
            if (has_range && (global_var_idx < start_var || global_var_idx >= no_var))
            {
                // Still need to advance fields_pos for <N> records
                variant_desc_t desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];
                if (desc.alt.find("<M>") == string::npos)
                    fields_pos++;
                continue;
            }

            variant_desc_t desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];

            // Use the previous fields entry for multi-allelic continuation records (<M>) since fields_pos
            // was already advanced when the <N> record was seen.
            int field_index = fields_pos;
            if (desc.alt.find("<M>") != string::npos && count > 0)
                field_index = fields_pos - 1;

            if (field_index < 0 || static_cast<uint32_t>(field_index) >= actual_variants) {
                logger->error("decompressAllTiled: field_index={} outside actual_variants={} (processed_variants={}, fields_pos={}, alt={})",
                              field_index, actual_variants, processed_variants, fields_pos, desc.alt);
                return 1;
            }

            if (has_sample_subset)
            {
                // Tiled sample subset: decode only the requested bytes per column-block and directly produce GTs.
                for (auto &ctx : subset_ctxs)
                {
                    const uint64_t pair_index = pair_base + static_cast<uint64_t>(ctx.cb) * block_variants + var_in_block;
                    if (!decode_vector_row_partial_bytes(pair_index * 2, curr_non_copy_vec_id_offset,
                                                         ctx.groups, ctx.byte_where,
                                                         ctx.bytes.data(), ctx.local_count))
                        return 1;
                    if (!decode_vector_row_partial_bytes(pair_index * 2 + 1, curr_non_copy_vec_id_offset,
                                                         ctx.groups, ctx.byte_where,
                                                         ctx.bytes.data() + ctx.local_count, ctx.local_count))
                        return 1;

                    for (uint32_t i = 0; i < ctx.local_count; ++i)
                    {
                        uint8_t b0 = map_t256[ctx.bytes[i]];
                        uint8_t b1 = map_t256[ctx.bytes[i + ctx.local_count]];
                        subset_str[ctx.out_hap[i]] = gt_lookup_table[b0][b1][ctx.perm_bit[i]];
                    }
                }
                SetVariantToRec(desc, all_fields_io[field_index], decompression_reader.keys, subset_str, output_haplotypes);
            }
            else
            {
                fill_n(decomp_data, full_vec_len * 2, 0);
                for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
                {
                    const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
                    const uint32_t col_vec_len = static_cast<uint32_t>(decompression_reader.col_block_vec_lens[cb]);
                    const uint64_t pair_index = pair_base + static_cast<uint64_t>(cb) * block_variants + var_in_block;

                    fill_n(col_decomp_perm.data(), col_vec_len * 2, 0);
                    decoded_vector_row(pair_index * 2, 0, curr_non_copy_vec_id_offset, col_vec_len, 0, col_decomp_perm.data());
                    decoded_vector_row(pair_index * 2 + 1, 0, curr_non_copy_vec_id_offset, col_vec_len, col_vec_len, col_decomp_perm.data());

                    fill_n(col_decomp.data(), col_vec_len * 2, 0);
                    decode_perm_rev(col_vec_len, col_vec_len, col_block_size, sort_perm_io[block_id][cb],
                                    col_decomp_perm.data(), col_decomp.data());

                    const uint32_t start_hap = decompression_reader.col_block_ranges[cb].first;
                    insert_block_bits(decomp_data, col_decomp.data(), full_vec_len, start_hap, col_block_size);
                    insert_block_bits(decomp_data + full_vec_len, col_decomp.data() + col_vec_len, full_vec_len, start_hap, col_block_size);
                }

                int vec1_start;
                int vec2_start = full_vec_len;
                for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
                {
                    lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                    tmp_arr[vec1_start] = *lookup_table_ptr;
                }
                memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
                if (trailing_bits)
                    memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);

                SetVariantToRec(desc, all_fields_io[field_index], decompression_reader.keys, my_str, haplotype_count);
            }
        }
        pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
    }

    if (cur_chunk_id == end_chunk_id && count)
    {
        logger->debug("decompressAllTiled: appending final multi-allelic record, count={}, has_sample_subset={}", count, has_sample_subset);
        if (has_sample_subset)
            appendVCFToRec(temp_desc, genotype, output_haplotypes, temp_fields, decompression_reader.keys);
        else
            appendVCFToRec(temp_desc, genotype, static_cast<uint32_t>(haplotype_count), temp_fields, decompression_reader.keys);
        logger->debug("decompressAllTiled: final multi-allelic record appended");
    }

    logger->info("Processed chunk {}", cur_chunk_id);

    logger->debug("decompressAllTiled: cleaning done_unique, size={}", done_unique.size());
    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();
    logger->debug("decompressAllTiled: done_unique cleared, returning");

    return 0;
}

//*****************************************************************************************************************
// Helper function to decode one variant's GT data in tiled mode
void Decompressor::decodeVariantGtTiled(uint32_t block_id, uint32_t var_in_block, uint32_t n_col_blocks,
                                        uint64_t pair_base, uint64_t curr_non_copy_vec_id_offset,
                                        uint32_t full_vec_len, vector<uint8_t> &my_str)
{
    const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());

    if (max_col_vec_len_ == 0)
    {
        for (auto len : decompression_reader.col_block_vec_lens)
            if (len > max_col_vec_len_)
                max_col_vec_len_ = len;
    }
    if (max_col_vec_len_ > 0)
    {
        const size_t need = static_cast<size_t>(max_col_vec_len_) * 2;
        if (col_decomp_perm_buf_.size() < need)
            col_decomp_perm_buf_.resize(need);
        if (col_decomp_buf_.size() < need)
            col_decomp_buf_.resize(need);
    }
    uint8_t *col_decomp_perm = col_decomp_perm_buf_.data();
    uint8_t *col_decomp = col_decomp_buf_.data();

    // Clear full decomp_data buffer
    fill_n(decomp_data, full_vec_len * 2, 0);

    // Process each column block
    for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
    {
        const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
        const uint32_t col_vec_len = static_cast<uint32_t>(decompression_reader.col_block_vec_lens[cb]);
        const uint64_t pair_index = pair_base + static_cast<uint64_t>(cb) * block_variants + var_in_block;

        fill_n(col_decomp_perm, col_vec_len * 2, 0);
        decoded_vector_row(pair_index * 2, 0, curr_non_copy_vec_id_offset, col_vec_len, 0, col_decomp_perm);
        decoded_vector_row(pair_index * 2 + 1, 0, curr_non_copy_vec_id_offset, col_vec_len, col_vec_len, col_decomp_perm);

        fill_n(col_decomp, col_vec_len * 2, 0);
        decode_perm_rev(col_vec_len, col_vec_len, col_block_size, sort_perm_io[block_id][cb],
                        col_decomp_perm, col_decomp);

        const uint32_t start_hap = decompression_reader.col_block_ranges[cb].first;
        insert_block_bits(decomp_data, col_decomp, full_vec_len, start_hap, col_block_size);
        insert_block_bits(decomp_data + full_vec_len, col_decomp + col_vec_len, full_vec_len, start_hap, col_block_size);
    }

    // Convert bit vectors to GT characters using lookup table
    int vec1_start;
    int vec2_start = full_vec_len;
    for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
    {
        lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
        tmp_arr[vec1_start] = *lookup_table_ptr;
    }
    memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
    if (trailing_bits)
        memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
}

//*****************************************************************************************************************
// Range decompression for tiled GT blocks
int Decompressor::decompressRangeTiled(const string &range)
{
    auto logger = LogManager::Instance().Logger();
    done_unique.clear();
    stored_unique.clear();

    const uint32_t n_col_blocks = decompression_reader.n_col_blocks;
    const uint32_t full_vec_len = decompression_reader.vec_len;
    const uint32_t row_block_variants = row_block_size ? row_block_size : haplotype_count;

    if (!n_col_blocks || !row_block_variants)
    {
        logger->error("Invalid tiling metadata (n_col_blocks={}, row_block_variants={})", n_col_blocks, row_block_variants);
        return 1;
    }

    if (fixed_variants_chunk_io.empty())
        return 0;

    vector<uint8_t> my_str(haplotype_count);
    bool skip_processing = false;

    // Calculate start pair_id for this chunk
    uint64_t start_pair_id = static_cast<uint64_t>(start_chunk_actual_pos) * n_col_blocks;
    uint64_t gt_index_base_pair_id = start_pair_id;
    if (has_fixed_fields_rb_dir_io &&
        fixed_fields_chunk_version_io >= GSC_FIXED_FIELDS_RB_VERSION_V2 &&
        chunk_variant_offset_io)
    {
        gt_index_base_pair_id += static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    }
    uint64_t curr_non_copy_vec_id_offset = gt_index_base_pair_id * 2 - rrr_rank_zeros_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_zeros_bit_vector[1](gt_index_base_pair_id) - rrr_rank_copy_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_copy_bit_vector[1](gt_index_base_pair_id);

    const bool has_range = (range != "");
    int64_t start_pos = std::numeric_limits<int64_t>::min();
    int64_t end_pos = std::numeric_limits<int64_t>::max();
    if (has_range)
    {
        if (cur_chunk_id - 1 == start_chunk_id)
            start_pos = range_1;
        if (cur_chunk_id == end_chunk_id)
            end_pos = range_2;
    }

    // Calculate pair_base for each row_block
    uint64_t pair_base = start_pair_id + static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    bool stop = false;

    const bool use_row_ranges = has_range && !row_block_ranges_io.empty();
    for (size_t block_id = 0; block_id < fixed_variants_chunk_io.size(); ++block_id)
    {
        const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());
        uint32_t start_idx = 0;
        uint32_t end_idx = block_variants;
        if (use_row_ranges)
        {
            const auto &range_pair = row_block_ranges_io[block_id];
            start_idx = std::min(range_pair.first, block_variants);
            end_idx = std::min(range_pair.second, block_variants);
            if (start_idx >= end_idx)
            {
                pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
                continue;
            }
        }

        for (uint32_t var_in_block = start_idx; var_in_block < end_idx; ++var_in_block)
        {
            variant_desc_t desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];
            if (has_range && !use_row_ranges)
            {
                const int64_t pos = static_cast<int64_t>(desc.pos);
                if (pos < start_pos)
                    continue;
                if (pos > end_pos)
                {
                    if (cur_chunk_id == end_chunk_id)
                        stop = true;
                    break;
                }
            }

            if (!params.out_AC_AN)
            {
                const bool qual_id_reject =
                    (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual ||
                     !(params.out_id == desc.id || params.out_id.empty()));
                if (qual_id_reject)
                    continue;
            }

            // Decode this variant
            decodeVariantGtTiled(block_id, var_in_block, n_col_blocks, pair_base,
                                 curr_non_copy_vec_id_offset, full_vec_len, my_str);

            // Apply filters if needed
            if (params.out_AC_AN)
            {
                uint32_t AN = haplotype_count;
                uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                    [](char c) { return c == '1' || c == '2' || c == '.'; });

                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;
                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }

            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual ||
                              !(params.out_id == desc.id || params.out_id.empty()));
            if (skip_processing)
                continue;

            SetVariant(desc, my_str, haplotype_count);
        }
        if (stop)
            break;
        pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
    }

    if (cur_chunk_id == end_chunk_id && count)
        appendVCF(temp_desc, genotype, haplotype_count);

    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();

    logger->info("Processed chunk {}", cur_chunk_id);
    return 0;
}

//*****************************************************************************************************************
// Sample subset decompression for tiled GT blocks
int Decompressor::decompressSampleSmartTiled(const string &range)
{
    auto logger = LogManager::Instance().Logger();
    done_unique.clear();
    stored_unique.clear();

    const uint32_t n_col_blocks = decompression_reader.n_col_blocks;
    const uint32_t row_block_variants = row_block_size ? row_block_size : haplotype_count;

    if (!n_col_blocks || !row_block_variants)
    {
        logger->error("Invalid tiling metadata");
        return 1;
    }

    if (fixed_variants_chunk_io.empty())
        return 0;

    uint32_t no_haplotypes = smpl.no_samples * decompression_reader.ploidy;
    vector<uint8_t> gt_variant_data(no_haplotypes);
    bool skip_processing = false;

    // Calculate start pair_id for this chunk
    uint64_t start_pair_id = static_cast<uint64_t>(start_chunk_actual_pos) * n_col_blocks;
    uint64_t gt_index_base_pair_id = start_pair_id;
    if (has_fixed_fields_rb_dir_io &&
        fixed_fields_chunk_version_io >= GSC_FIXED_FIELDS_RB_VERSION_V2 &&
        chunk_variant_offset_io)
    {
        gt_index_base_pair_id += static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;
    }
    uint64_t curr_non_copy_vec_id_offset = gt_index_base_pair_id * 2 - rrr_rank_zeros_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_zeros_bit_vector[1](gt_index_base_pair_id) - rrr_rank_copy_bit_vector[0](gt_index_base_pair_id) -
                                           rrr_rank_copy_bit_vector[1](gt_index_base_pair_id);

    const bool has_range = (range != "");
    int64_t start_pos = std::numeric_limits<int64_t>::min();
    int64_t end_pos = std::numeric_limits<int64_t>::max();
    if (has_range)
    {
        if (cur_chunk_id - 1 == start_chunk_id)
            start_pos = range_1;
        if (cur_chunk_id == end_chunk_id)
            end_pos = range_2;
    }

    vector<vector<pair<uint32_t, uint32_t>>> selected_haps_by_cb(n_col_blocks);
    uint32_t assigned = 0;
    for (uint32_t g = 0; g < smpl.no_samples; ++g)
    {
        for (uint32_t p = 0; p < decompression_reader.ploidy; ++p)
        {
            uint32_t out_hap = g * decompression_reader.ploidy + p;
            uint32_t orig_hap = sampleIDs[g] * decompression_reader.ploidy + p;
            int found_cb = -1;
            for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
            {
                uint32_t start = decompression_reader.col_block_ranges[cb].first;
                uint32_t size = decompression_reader.col_block_ranges[cb].second;
                if (orig_hap >= start && orig_hap < start + size)
                {
                    found_cb = static_cast<int>(cb);
                    uint32_t local_orig = orig_hap - start;
                    selected_haps_by_cb[cb].push_back(make_pair(out_hap, local_orig));
                    assigned++;
                    break;
                }
            }
            if (found_cb < 0)
            {
                logger->error("decompressSampleSmartTiled: failed to map orig_hap={} to any column block", orig_hap);
                return 1;
            }
        }
    }
    if (assigned != no_haplotypes)
    {
        logger->error("decompressSampleSmartTiled: assigned_haplotypes={} != no_haplotypes={}", assigned, no_haplotypes);
        return 1;
    }

    struct ColSubsetCtx
    {
        uint32_t cb = 0;
        uint32_t local_count = 0;
        vector<uint32_t> out_hap;
        vector<uint8_t> perm_bit;
        vector<pair<uint32_t, uint32_t>> byte_where; // (byte_no, local_index)
        vector<ByteGroup> groups;
        vector<uint8_t> bytes; // [local_count*2] for parity0 + parity1
    };

    // Calculate pair_base for each row_block
    uint64_t pair_base = start_pair_id + static_cast<uint64_t>(chunk_variant_offset_io) * n_col_blocks;

    bool stop = false;
    for (size_t block_id = 0; block_id < fixed_variants_chunk_io.size(); ++block_id)
    {
        const uint32_t block_variants = static_cast<uint32_t>(fixed_variants_chunk_io[block_id].data_compress.size());

        // Validate sort_perm_io access
        if (block_id >= sort_perm_io.size())
        {
            logger->error("block_id {} >= sort_perm_io.size() {}", block_id, sort_perm_io.size());
            return 1;
        }
        if (sort_perm_io[block_id].size() != n_col_blocks)
        {
            logger->error("sort_perm_io[{}].size()={} != n_col_blocks={}", block_id, sort_perm_io[block_id].size(), n_col_blocks);
            return 1;
        }

        vector<ColSubsetCtx> subset_ctxs;
        subset_ctxs.reserve(n_col_blocks);
        for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
        {
            if (selected_haps_by_cb[cb].empty())
                continue;

            const uint32_t col_block_size = decompression_reader.col_block_ranges[cb].second;
            vector<uint32_t> orig_to_perm(col_block_size);
            reverse_perm(sort_perm_io[block_id][cb], orig_to_perm, static_cast<int>(col_block_size));

            ColSubsetCtx ctx;
            ctx.cb = cb;
            ctx.local_count = static_cast<uint32_t>(selected_haps_by_cb[cb].size());
            ctx.out_hap.resize(ctx.local_count);
            ctx.perm_bit.resize(ctx.local_count);
            ctx.byte_where.reserve(ctx.local_count);
            ctx.bytes.resize(static_cast<size_t>(ctx.local_count) * 2);

            for (uint32_t i = 0; i < ctx.local_count; ++i)
            {
                uint32_t out_h = selected_haps_by_cb[cb][i].first;
                uint32_t local_orig = selected_haps_by_cb[cb][i].second;
                uint32_t perm_idx = orig_to_perm[local_orig];
                ctx.out_hap[i] = out_h;
                ctx.perm_bit[i] = static_cast<uint8_t>(perm_idx & 7);
                ctx.byte_where.emplace_back(static_cast<uint32_t>(perm_idx >> 3), i);
            }

            sort(ctx.byte_where.begin(), ctx.byte_where.end(),
                 [](const auto &a, const auto &b) { return a.first < b.first || (a.first == b.first && a.second < b.second); });

            for (uint32_t i = 0; i < ctx.local_count;)
            {
                uint32_t byte_no = ctx.byte_where[i].first;
                uint32_t j = i + 1;
                while (j < ctx.local_count && ctx.byte_where[j].first == byte_no)
                    ++j;
                ctx.groups.push_back(ByteGroup{byte_no, i, j});
                i = j;
            }

            subset_ctxs.emplace_back(std::move(ctx));
        }

        for (uint32_t var_in_block = 0; var_in_block < block_variants; ++var_in_block)
        {
            variant_desc_t desc = fixed_variants_chunk_io[block_id].data_compress[var_in_block];
            if (has_range)
            {
                const int64_t pos = static_cast<int64_t>(desc.pos);
                if (pos < start_pos)
                    continue;
                if (pos > end_pos)
                {
                    if (cur_chunk_id == end_chunk_id)
                        stop = true;
                    break;
                }
            }

            // Decode only the requested bytes per column-block and directly produce GTs.
            for (auto &ctx : subset_ctxs)
            {
                const uint64_t pair_index = pair_base + static_cast<uint64_t>(ctx.cb) * block_variants + var_in_block;
                if (!decode_vector_row_partial_bytes(pair_index * 2, curr_non_copy_vec_id_offset,
                                                     ctx.groups, ctx.byte_where,
                                                     ctx.bytes.data(), ctx.local_count))
                    return 1;
                if (!decode_vector_row_partial_bytes(pair_index * 2 + 1, curr_non_copy_vec_id_offset,
                                                     ctx.groups, ctx.byte_where,
                                                     ctx.bytes.data() + ctx.local_count, ctx.local_count))
                    return 1;

                for (uint32_t i = 0; i < ctx.local_count; ++i)
                {
                    uint8_t b0 = map_t256[ctx.bytes[i]];
                    uint8_t b1 = map_t256[ctx.bytes[i + ctx.local_count]];
                    gt_variant_data[ctx.out_hap[i]] = gt_lookup_table[b0][b1][ctx.perm_bit[i]];
                }
            }

            // Apply filters if needed
            if (params.out_AC_AN)
            {
                uint32_t AN = no_haplotypes;
                uint32_t AC = std::count_if(gt_variant_data.begin(), gt_variant_data.end(),
                                    [](char c) { return c == '1' || c == '2' || c == '.'; });

                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;
                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }

            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual ||
                              !(params.out_id == desc.id || params.out_id.empty()));
            if (skip_processing)
                continue;

            SetVariant(desc, gt_variant_data, no_haplotypes);
        }

        pair_base += static_cast<uint64_t>(block_variants) * n_col_blocks;
        if (stop)
            break;
    }

    if (cur_chunk_id == end_chunk_id && count)
        appendVCF(temp_desc, genotype, no_haplotypes);

    for (auto &it : done_unique)
        delete[] it.second;
    done_unique.clear();

    logger->info("Processed chunk {}", cur_chunk_id);
    return 0;
}

//*****************************************************************************************************************
int Decompressor::decompressAll(){

    if (!decompression_reader.useLegacyPath)
        return decompressAllTiled();

    auto logger = LogManager::Instance().Logger();
    done_unique.clear();
    stored_unique.clear();
    uint32_t cur_block_id = 0;
    uint32_t c_out_line = 0;
    uint32_t no_var = 0;
    uint32_t start_var = 0;
    int vec1_start,vec2_start=decompression_reader.vec_len;
    fields_pos  = 0;
    vector<uint8_t> my_str(standard_block_size);
    // vector<uint32_t> rev_perm(standard_block_size);  
    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos * 2 - rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) -rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);

    no_var = end_chunk_actual_pos - start_chunk_actual_pos;

    size_t cur_var;

    for (cur_var = start_var; cur_var + standard_block_size <= no_var; cur_var += standard_block_size )
    {
        
        cur_block_id = cur_var / standard_block_size;
        for(size_t i = 0; i < standard_block_size; i++){
            
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);

            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
            
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            }   
            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[i];

            SetVariantToRec(desc, all_fields_io[fields_pos], decompression_reader.keys, my_str,standard_block_size);
            
        
        }
    }
    if(no_var % standard_block_size)
    {
        cur_block_id = cur_var / standard_block_size;            
        // reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size); //2024.1.16
        
        for(;cur_var < no_var;++cur_var){
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
            // if(!(no_var%standard_block_size))
            //     decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
            // else
            //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);

            // for(int i = 0;i<rev_perm.size();i++)
            //     std::cerr<<rev_perm[i]<<" ";
            // std::cerr<<endl;
            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data); //2024.1.16 rev_perm ---> sort_perm_io[cur_block_id][0]
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            } 
                
            c_out_line = cur_var % standard_block_size;

            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            
            SetVariantToRec(desc, all_fields_io[fields_pos], decompression_reader.keys, my_str, standard_block_size);
        }
    }
    
    if(cur_chunk_id == end_chunk_id && count){
        
        appendVCFToRec(temp_desc, genotype, static_cast<uint32_t>(standard_block_size), temp_fields, decompression_reader.keys);
    }
    
    logger->info("Processed chunk {}", cur_chunk_id);

    for (auto &it : done_unique)
        delete[] it.second;

    done_unique.clear();

    return 0;
    
}
// // *****************************************************************************************************************
// Decompress by range
int Decompressor::decompressRange(const string &range)
{
    auto logger = LogManager::Instance().Logger();
    if (!decompression_reader.useLegacyPath)
    {
        return decompressRangeTiled(range);
    }
    done_unique.clear();
    stored_unique.clear();
    uint32_t cur_block_id = 0;
    uint32_t c_out_line = 0;
    uint32_t no_var = 0;
    uint32_t start_var = 0;
    int vec1_start,vec2_start;
    bool skip_processing = false;
    // vector<uint32_t> rev_perm(standard_block_size);  
    vector<uint8_t> my_str(standard_block_size);

    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos * 2 - rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) -rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);
    
    if (range == "")
    {
        no_var = end_chunk_actual_pos - start_chunk_actual_pos;

        size_t cur_var;
        for (cur_var = start_var; cur_var + standard_block_size <= no_var; cur_var += standard_block_size )
        {
            cur_block_id = cur_var / standard_block_size;
            // reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size);
            for(size_t i = 0; i < standard_block_size; i++){
                variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[i];
                const bool needs_ac_an = params.out_AC_AN;
                if (!needs_ac_an)
                {
                    const bool qual_id_reject =
                        (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual ||
                         !(params.out_id == desc.id || params.out_id.empty()));
                    if (qual_id_reject)
                    {
                        cur_vec_id += 2;
                        continue;
                    }
                }

                vec2_start = decompression_reader.vec_len;
                fill_n(decomp_data, decompression_reader.vec_len*2, 0);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
                decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
            
                for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
                {
                    lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);

                    tmp_arr[vec1_start] = *lookup_table_ptr;
                }
                memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
                if(trailing_bits)
                {
                    memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                    
                }
            
                // string my_str = "";
                // getRangeGT(decomp_data, standard_block_size, my_str);
                // for(int i = 0;i<(int)decompression_reader.vec_len*2;i++)
                //     decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
                // getRangeGT(decomp_data_perm,rev_perm,standard_block_size, my_str);    
                // c_out_line = i;

                if(params.out_AC_AN){
                    
                    uint32_t AN = standard_block_size;
                    uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                        [](char c)
                                        { return c == '1' || c == '2' || c == '.'; });
                    
                    skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);

                    if (skip_processing)
                        continue;

                    desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);

                    
                } 
                if (params.out_AC_AN)
                    skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

                if (skip_processing)
                    continue;  

                SetVariant(desc, my_str, standard_block_size);
            
            }
        }
        if(no_var % standard_block_size)
        {
            cur_block_id = cur_var / standard_block_size;
                           
            // reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size); 
    
            for(;cur_var < no_var;++cur_var){
                c_out_line = cur_var % standard_block_size;
                variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
                const bool needs_ac_an = params.out_AC_AN;
                if (!needs_ac_an)
                {
                    const bool qual_id_reject =
                        (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual ||
                         !(params.out_id == desc.id || params.out_id.empty()));
                    if (qual_id_reject)
                    {
                        cur_vec_id += 2;
                        continue;
                    }
                }

                vec2_start = decompression_reader.vec_len;
                fill_n(decomp_data, decompression_reader.vec_len*2, 0);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
                // if(!(no_var%standard_block_size))
                //     decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
                // else
                //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);

                decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data); //2024.1.16 rev_perm ---> sort_perm_io[cur_block_id][0]
                for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
                {
                    lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                    tmp_arr[vec1_start] = *lookup_table_ptr;
                }
                memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
                if(trailing_bits)
                {
                    memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                    
                } 
                // if(!(no_var%standard_block_size))
                //     for(int i = 0;i< (int)decompression_reader.vec_len*2;i++)
                //         decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
                        
                // // string my_str = "";
                // // vector<uint8_t> my_str(standard_block_size);
                // // getRangeGT(decomp_data_perm, standard_block_size, my_str);
                // getRangeGT(decomp_data_perm,rev_perm,standard_block_size, my_str);
                    
                if(params.out_AC_AN){
                    
                    uint32_t AN = standard_block_size;
                    uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                        [](char c)
                                        { return c == '1' || c == '2' || c == '.'; });
                    
                    skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);

                    if (skip_processing)
                        continue;

                    desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);

                    
                }   
                
                if (params.out_AC_AN)
                    skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

                if (skip_processing)
                    continue;
                // for(int i = 0; i < standard_block_size; i++){
                //     std::cerr << my_str[i] << " ";
                // }
                // std::cerr << endl; 
                SetVariant(desc, my_str, standard_block_size);
            }
        }

        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, standard_block_size);
        for (auto &it : done_unique)
            delete[] it.second;

        done_unique.clear();
    }
    else
    {
       
        int start_block, end_block;
        int start_position, end_position;
        
        if(cur_chunk_id-1 == start_chunk_id){
            
            calculate_start_position(start_block,start_position);
            start_var = uint32_t(start_block * standard_block_size + start_position);
            
        }else{
            start_var = 0;
        }
        if(cur_chunk_id == end_chunk_id){
            
            calculate_end_position(end_block,end_position);
            
            no_var = uint32_t(end_block * standard_block_size + end_position);
        }else{

            no_var = end_chunk_actual_pos - start_chunk_actual_pos;
        }       

        cur_vec_id = (start_chunk_actual_pos + start_var) * 2;
        // bool  tail_flag = true; 
        for (size_t cur_var = start_var;cur_var < no_var; cur_var++)
        {

            cur_block_id = cur_var / standard_block_size;
            // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size && tail_flag){
            //     reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size);
            //     sort_perm_io[cur_block_id][0] = rev_perm;
            //     tail_flag = false;
            // }                                
            
            c_out_line = cur_var % standard_block_size;
            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            const bool needs_ac_an = params.out_AC_AN;
            if (!needs_ac_an)
            {
                const bool qual_id_reject =
                    (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual ||
                     !(params.out_id == desc.id || params.out_id.empty()));
                if (qual_id_reject)
                {
                    cur_vec_id += 2;
                    continue;
                }
            }

            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
            // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() == standard_block_size)
            //     
            // else
            //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);
            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id][0], decomp_data_perm, decomp_data);
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            } 
                
            // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() == standard_block_size){
            //     for(int i = 0;i< (int)decompression_reader.vec_len*2;i++)
            //         decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
            // }
            // getRangeGT(decomp_data_perm,rev_perm,standard_block_size, my_str);

            // // string my_str = "";

            // // getRangeGT(decomp_data, standard_block_size, my_str);

            if(params.out_AC_AN){
                uint32_t AN = standard_block_size;

                uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                    [](char c)
                                    { return c == '1' || c == '2' || c == '.'; });
                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;
                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }

            if (params.out_AC_AN)
                skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

            if (skip_processing)
                continue;
            SetVariant(desc, my_str, standard_block_size);
            
        }
        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, standard_block_size);
        for (auto &it : done_unique)
            delete[] it.second;
        done_unique.clear();
    }
    logger->info("Processed chunk {}", cur_chunk_id);
    // if (decomp_data)
    //     delete[] decomp_data;
    // if (decomp_data_perm)
    //     delete[] decomp_data_perm;
    // if (new_decomp_data_perm)
    //     delete[] new_decomp_data_perm;

    return 0;
}
// // *****************************************************************************************************************
int Decompressor::decompressSampleSmart(const string &range)
{
    auto logger = LogManager::Instance().Logger();
    if (!decompression_reader.useLegacyPath)
    {
        return decompressSampleSmartTiled(range);
    }
    
    uint32_t no_var;

    size_t cur_block_id;

    size_t c_out_line;

    if (smpl.no_samples > NO_SAMPLE_THRESHOLD) // || range != "")
    {
        return decompressRangeSample(range);
    }
    else if (range != "")
    {
        ;
    }
    uint32_t no_haplotypes = smpl.no_samples * decompression_reader.ploidy;

    std::vector<std::pair<uint32_t, uint32_t>> whichByte_whereInRes(no_haplotypes + 1); // +1 for guard

    bool is_unique = false;
    
    bool skip_processing = false;

    uint32_t unique_pos = 0, unique_pos_first_in_block = 0;

    uint64_t curr_zeros = 0, curr_copy = 0;

    vector<uint32_t> rev_perm(standard_block_size);

    vector<uint8_t> gt_variant_data(no_haplotypes);

    uint32_t  prev_block_id = 0xFFFF;

    uint32_t where; // uchar where;

    uint32_t ind_id_orig;

    unique_pos = 0;

    curr_zeros = 0, curr_copy = 0;

    uint32_t written_records = 0;

    uint8_t *resUnique = nullptr;

    resUnique = new uint8_t[standard_block_size*2 * no_haplotypes];

    uint8_t *resAll = nullptr;

    resAll = new uint8_t[2 * no_haplotypes];

    uint32_t first_vec_in_block = 0;

    // int var_block_num=fixed_variants_chunk_io.size();
    uint64_t prev_chr_zeros_copy = rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) + rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) +
                                       rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) + rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);
    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos*2 - prev_chr_zeros_copy;
    
    if (range != "")
    {

        int start_block, end_block;
        int start_position, end_position;
        uint32_t start_var;
        if(cur_chunk_id-1 == start_chunk_id){
            calculate_start_position(start_block,start_position);
            start_var = uint32_t(start_block * standard_block_size + start_position);
            
        }else{
            start_var = 0;
        }
        if(cur_chunk_id == end_chunk_id){
            calculate_end_position(end_block,end_position);
            no_var = uint32_t(end_block * standard_block_size + end_position);
        }else{
            no_var = end_chunk_actual_pos - start_chunk_actual_pos;
        }
        
        cur_vec_id = (start_var + start_chunk_actual_pos) * 2;
  
        for (size_t cur_var = start_var;cur_var < no_var ; cur_var++)
        {
            cur_block_id = cur_var / standard_block_size;

            if (cur_block_id != prev_block_id) // Get perm and find out which bytes need decoding
            {
                first_vec_in_block = start_chunk_actual_pos*2 + cur_block_id * standard_block_size*2;
                unique_pos = 0;
                if (prev_block_id == 0xFFFF) // to set curr_zeros, curr_copy (first processed block)
                {
                    uint8_t parity = first_vec_in_block & 1;
                    uint32_t id = first_vec_in_block >> 1;
                    curr_zeros = rrr_rank_zeros_bit_vector[0](id + ((parity))) + rrr_rank_zeros_bit_vector[1](id);
                    curr_copy = rrr_rank_copy_bit_vector[0](id + ((parity))) + rrr_rank_copy_bit_vector[1](id);
                }
                unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
                // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size) 
                //     rev_perm = sort_perm_io[cur_block_id][0];
                // else
                reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size);

            //                 for(int i = 0;i<rev_perm.size();i++)
            //     std::cerr<<rev_perm[i]<<" ";
            // std::cerr<<endl;
                for (uint32_t s = 0; s < smpl.no_samples; s++)
                {
                    for (uint32_t p = 0; p < decompression_reader.ploidy; p++)
                    {
                        ind_id_orig = sampleIDs[s] * decompression_reader.ploidy + p;
                        where = s * decompression_reader.ploidy + p;
                        // std::cerr<<perm[ind_id_orig]<<endl;
                        whichByte_whereInRes[where] = std::make_pair(rev_perm[ind_id_orig] >> 3, where); // now original_id_only
                    }
                }
                whichByte_whereInRes[no_haplotypes] = std::make_pair(0xFFFFFFFF, no_haplotypes); // guard
                // Sort by byte_no
                std::sort(whichByte_whereInRes.begin(), whichByte_whereInRes.end());
                prev_block_id = cur_block_id;

                // Get vectors from all block
                // copies only within the same block, so only part of resUnique (within the same block and using the same perm) will be used - no need to fix previous unique bytes (these are in different blocks, with diff perm)
                
                for (uint64_t i = first_vec_in_block; i < cur_vec_id; i++)
                {

                    extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy, curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block, false);
                    if (is_unique)
                    {
                        memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                        unique_pos++;
                    }
                }
            }
        //     //   else
            {

        //         // previous vectors in block are decompressed already
        
                for (uint64_t i = cur_vec_id; i <= cur_vec_id + 1; i++)
                {
                    extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy,curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block);

                    if (is_unique)
                    {
                        
                        //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
                        {
                            memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                        }
                        unique_pos++;
      
                    }
                }
            }
            written_records++;
            prev_block_id = cur_block_id;
            cur_vec_id += 2;
            int r_i = 0;
            // string my_str;

            size_t c_out_line = cur_var % standard_block_size;

            // if (fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size){
            //     for (int g = 0; g < (int)smpl.no_samples; g++)
            //     {
            //         int samples_start = sampleIDs[g] * decompression_reader.ploidy;
            //         for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
            //         {
            //             uint32_t gt_pos = rev_perm[samples_start+p] % 8;
            //             gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
            //             // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos,standard_block_size); 

            //         }
            //     }
            // }
            // else{
            for (int g = 0; g < (int)smpl.no_samples; g++)
            {
                int samples_start = sampleIDs[g] * decompression_reader.ploidy;
                for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
                {
                    uint32_t gt_pos = rev_perm[samples_start+ p] % 8;
                    resAll[r_i] = map_t256[resAll[r_i]];
                    resAll[r_i + no_haplotypes] = map_t256[resAll[r_i + no_haplotypes]];
                    gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
                    // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos, standard_block_size); 

                }
            }
            // }
            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            if(params.out_AC_AN){
                uint32_t AN = no_haplotypes;

                uint32_t AC = std::count_if(gt_variant_data.begin(), gt_variant_data.end(),
                                    [](char c)
                                    { return c == '1' || c == '2' || c == '.'; });
                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;

                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }
            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

            if (skip_processing)
                continue;          
            SetVariant(desc, gt_variant_data, static_cast<size_t> (no_haplotypes));
            
        }
        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, static_cast<size_t> (no_haplotypes));

    }
    else
    {
        no_var = end_chunk_actual_pos - start_chunk_actual_pos;

        for (size_t cur_var = 0; cur_var < no_var; cur_var++)
        {
        
            // Permutations
            cur_block_id = cur_var  / standard_block_size;

            if (cur_block_id != prev_block_id) // get perm and find out which bytes need decoding
            {
                first_vec_in_block = start_chunk_actual_pos*2 +cur_block_id * standard_block_size*2;
                unique_pos = 0;
                if (prev_block_id == 0xFFFF) // to set curr_zeros, curr_copy (first processed block)
                {
                    uint8_t parity = first_vec_in_block & 1;
                    uint32_t id = first_vec_in_block >> 1;
                    curr_zeros = rrr_rank_zeros_bit_vector[0](id + ((parity))) + rrr_rank_zeros_bit_vector[1](id);
                    curr_copy = rrr_rank_copy_bit_vector[0](id + ((parity))) + rrr_rank_copy_bit_vector[1](id);
                }

                unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
                // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size)  
                //     rev_perm = sort_perm_io[cur_block_id][0];                                                  
                // else
                reverse_perm(sort_perm_io[cur_block_id][0], rev_perm, standard_block_size);

                for (uint32_t s = 0; s < smpl.no_samples; s++)
                {
                    for (uint32_t p = 0; p < decompression_reader.ploidy; p++)
                    {
                        ind_id_orig = sampleIDs[s] * decompression_reader.ploidy + p;
                        where = s * decompression_reader.ploidy + p;

                        whichByte_whereInRes[where] = std::make_pair(rev_perm[ind_id_orig] >> 3, where); // now original_id_only
                    }
                }
                whichByte_whereInRes[no_haplotypes] = std::make_pair(0xFFFFFFFF, no_haplotypes); // guard

                std::sort(whichByte_whereInRes.begin(), whichByte_whereInRes.end());
                // std::cerr<<"prev_block_id:"<<prev_block_id<<endl;
                // std::cerr<<"start_chunk_actual_pos:"<<start_chunk_actual_pos<<endl;
                prev_block_id = cur_block_id;
                // Get vectors from all block
                
                // for (uint64_t i = cur_vec_id; i <= cur_vec_id + 1; i++)
                // {
                //     extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy,curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block);

                //     if (is_unique)
                //     {
                //         memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                //         unique_pos++;
                //     }
                // }
            }
            // else
            {
                for (uint64_t i = cur_vec_id; i <= cur_vec_id + 1; i++)
                {

                    extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy,curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block);

                    if (is_unique)
                    {
                        memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                        unique_pos++;
                        
                    }
                }
            }
            cur_vec_id += 2;
            prev_block_id = cur_block_id;
            c_out_line = cur_var % standard_block_size;
            int r_i = 0;
            // if (fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size){
            //     for (int g = 0; g < (int)smpl.no_samples; g++)
            //     {
            //         int samples_start = sampleIDs[g] * decompression_reader.ploidy;
            //         for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
            //         {
            //             uint32_t gt_pos = rev_perm[samples_start+ p] % 8;
            //             gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
            //             // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos,standard_block_size); 

            //         }
            //     }
            // }
            // else{
            for (int g = 0; g < (int)smpl.no_samples; g++)
            {
                int samples_start = sampleIDs[g] * decompression_reader.ploidy;
                for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
                {
                    uint32_t gt_pos = rev_perm[samples_start+ p] % 8;
                    resAll[r_i] = map_t256[resAll[r_i]];
                    resAll[r_i + no_haplotypes] = map_t256[resAll[r_i + no_haplotypes]];
                    gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
                    // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos, standard_block_size); 

                }
            }
            // }

            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            
            if(params.out_AC_AN){

                uint32_t AN = no_haplotypes;

                uint32_t AC = std::count_if(gt_variant_data.begin(), gt_variant_data.end(),
                                    [](char c)
                                    { return c == '1' || c == '2' || c == '.'; });
                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;

                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }
            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

            if (skip_processing)
                continue;         
            SetVariant(desc, gt_variant_data, static_cast<size_t> (no_haplotypes));

        }
        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, static_cast<size_t> (no_haplotypes));
    }
    logger->info("Processed chunk {}", cur_chunk_id);
    if (resUnique)
        delete[] resUnique;
    if (resAll)
        delete[] resAll;

    return 0;
}
int Decompressor::decompressRangeSample(const string &range)
{
    return 0;
}
// // ***************************************************************************************************************************************************************
// // Get the byte corresponding to the sample genotype in the matrix
uint8_t Decompressor::extract_partial_bytes(uint64_t vec_id, std::vector<std::pair<uint32_t, uint32_t>> &whichByte_whereInRes, uint8_t *resUnique, bool &is_uniqe_id, uint64_t &curr_zeros, uint64_t &curr_copy, uint64_t curr_non_copy_vec_id_offset, uint8_t *resAll, uint32_t unique_pos_first_in_block, bool full_decode) // uint8_t Decompressor::extract_partial_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uint8_t> > & whichByte_whereInRes, uint8_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy, uint8_t * resAll)
{
    uint32_t next_haplotype = 0;
    uint32_t last_byte = whichByte_whereInRes.back().first; // Last byte to decode
    uint32_t tmp = 0;
    uint8_t parity = vec_id & 1;
    uint64_t vector = vec_id >> 1, curr_non_copy_vec_id;
    uint64_t no_haplotypes = whichByte_whereInRes.size() - 1;
    uint64_t vec_start = parity * no_haplotypes; //(vec_id - first_vec_in_block) * no_haplotypes;
    // size_t index_pos_size;
    size_t id_pos_size;
    uint64_t index_start = 0;
    size_t cur_i = 0;
    size_t count = 0;
    // uint8_t prev_index1 = 0, prev_index2 = 0;
    uint32_t cur_index;
    size_t prev_i = 999999;
    uint32_t cur_id = 0;
    int cur_pos = 0;
    if (copy_bit_vector[parity][vector]) // If vector is a copy od other vector (certainly it is placed within the same block)
    {
        if (!full_decode)
        {
            is_uniqe_id = false;
            curr_copy++;
            return 0;
        }
        unsigned long long bit_pos = (curr_copy)*decompression_reader.used_bits_cp;
        decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);

        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)

        curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;

        curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
        // std::cerr<<"unique_pos_first_in_block:"<<curr_non_copy_vec_id<<" "<<unique_pos_first_in_block<<endl;
        is_uniqe_id = false;
        curr_copy++;
        memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id - unique_pos_first_in_block) * no_haplotypes, no_haplotypes);
        

        return 0;
    }
    else if (zeros_bit_vector[parity][vector]) 
    {
        is_uniqe_id = false;
        curr_zeros++;
        if (!full_decode)
            return 0;
        fill_n(resAll + vec_start, no_haplotypes, 0);
        return 0;
    }
    else
    {
        curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
        is_uniqe_id = true;
    }
    fill_n(resAll + no_haplotypes*parity, no_haplotypes, 0);

    if (curr_non_copy_vec_id == 0)
        index_start = 0;
    else
        index_start = id_pos[curr_non_copy_vec_id - 1] + 1;
    const size_t index_end = static_cast<size_t>(id_pos[curr_non_copy_vec_id]);
    if (index_end < index_start || index_end > decompress_gt_indexes_io.size())
        return 0;

    vint_code::DecodeArray(decompress_gt_indexes_io.data() + index_start,
                           index_end - index_start,
                           sparse_matrix_cols);
    id_pos_size = sparse_matrix_cols.size();
    int prev = -1;
    while (whichByte_whereInRes[next_haplotype].first < last_byte && next_haplotype < whichByte_whereInRes.back().second)
    {
        if (cur_i < id_pos_size)
        {
            if (cur_i != prev_i)
            {
                cur_index = sparse_matrix_cols[cur_i]+prev;
                prev = cur_index;
                cur_pos = cur_index % 8;
                cur_id = cur_index / 8;
            }
            prev_i = cur_i;
            if (cur_id < whichByte_whereInRes[next_haplotype].first)
                cur_i += 1;
            else if (cur_id == whichByte_whereInRes[next_haplotype].first)
            {
                resAll[vec_start + whichByte_whereInRes[next_haplotype].second] += perm_lut8[cur_pos];
                cur_i += 1;
                count++;
            }
            else
            {
                next_haplotype++;
                if (count)
                {
                    if (next_haplotype < whichByte_whereInRes.back().second)
                        resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
                }
                else
                    resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
            }
        }
        else
        {
            next_haplotype++;
            if (count)
            {
                if (next_haplotype < whichByte_whereInRes.back().second)
                    resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
            }
            else
                resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
        }
    }
    // next_haplotype = 0;



    // if (parity)
    // {
    //     if (copy_bit_vector[parity][vector]) // If vector is a copy od other vector (certainly it is placed within the same block)
    //     {
    //         if (!full_decode)
    //         {
    //             is_uniqe_id = false;
    //             curr_copy++;
    //             return 0;
    //         }
    //         unsigned long long bit_pos = (curr_copy)*decompression_reader.used_bits_cp;
    //         decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);

    //         // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)

    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;

    //         curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
    //         // std::cerr<<"unique_pos_first_in_block:"<<curr_non_copy_vec_id<<" "<<unique_pos_first_in_block<<endl;
    //         is_uniqe_id = false;
    //         curr_copy++;
    //         memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id - unique_pos_first_in_block) * no_haplotypes, no_haplotypes);
            

    //         return 0;
    //     }
    //     else if (zeros_bit_vector[1][vector])
    //     {
    //         std::cerr<<zeros_bit_vector[1][vector]<<":"<<zeros_bit_vector[0][vector]<<":"<<vector<<":"<<vec_id<<endl;
    //         is_uniqe_id = false;
    //         curr_zeros++;
    //         if (!full_decode)
    //             return 0;
    //         fill_n(resAll + vec_start, no_haplotypes, 0);
    //         return 0;
    //     }
    //     else
    //     {
    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
    //         is_uniqe_id = true;
    //     }
    //     fill_n(resAll + no_haplotypes, no_haplotypes, 0);

    //     if (curr_non_copy_vec_id == 0)
    //         index_start = 0;
    //     else
    //         index_start = id_pos[curr_non_copy_vec_id - 1] + 1;
    //     curr_gt_indexes.assign(decompress_gt_indexes_io.begin()+ index_start, decompress_gt_indexes_io.begin() + id_pos[curr_non_copy_vec_id]);
    //     sparse_matrix_cols = vint_code::DecodeArray(curr_gt_indexes);
    //     id_pos_size = sparse_matrix_cols.size();
    //     int prev = -1;
    //     while (whichByte_whereInRes[next_haplotype].first < last_byte && next_haplotype < whichByte_whereInRes.back().second)
    //     {
    //         if (cur_i < id_pos_size)
    //         {
    //             if (cur_i != prev_i)
    //             {
    //                 cur_index = sparse_matrix_cols[cur_i]+prev;
    //                 prev = cur_index;
    //                 cur_pos = cur_index % 8;
    //                 cur_id = cur_index / 8;
    //             }
    //             prev_i = cur_i;
    //             if (cur_id < whichByte_whereInRes[next_haplotype].first)
    //                 cur_i += 1;
    //             else if (cur_id == whichByte_whereInRes[next_haplotype].first)
    //             {
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype].second] += perm_lut8[cur_pos];
    //                 cur_i += 1;
    //                 count++;
    //             }
    //             else
    //             {
    //                 next_haplotype++;
    //                 if (count)
    //                 {
    //                     if (next_haplotype < whichByte_whereInRes.back().second)
    //                         resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //                 }
    //                 else
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //             }
    //         }
    //         else
    //         {
    //             next_haplotype++;
    //             if (count)
    //             {
    //                 if (next_haplotype < whichByte_whereInRes.back().second)
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //             }
    //             else
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //         }
    //     }
    //     next_haplotype = 0;
    // }
    // else if (!zeros_bit_vector[0][vector])
    // {
    //     is_uniqe_id = false;
    //     curr_zeros++;
    //     if (!full_decode)
    //         return 0;
        
    //     fill_n(resAll + vec_start, no_haplotypes, 0);
    //     return 0;
    // }
    // else
    // {
    //     if (copy_bit_vector[parity][vector])
    //     {
    //         if (!full_decode)
    //         {
    //             is_uniqe_id = false;
    //             curr_copy++;
    //             return 0;
    //         }
    //         unsigned long long bit_pos = (curr_copy)*decompression_reader.used_bits_cp;
    //         decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);
    //         // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
    //         curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;

    //         is_uniqe_id = false;
    //         curr_copy++;
    //         memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id - unique_pos_first_in_block) * no_haplotypes, no_haplotypes);
            
    //         return 0;
    //     }
    //     else
    //     {
    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
    //         is_uniqe_id = true;
    //     }

    //     fill_n(resAll, no_haplotypes, 0);

    //     if (curr_non_copy_vec_id == 0)
    //         index_start = 0;
    //     else
    //         index_start = id_pos[curr_non_copy_vec_id - 1] + 1;
    //     curr_gt_indexes.assign(decompress_gt_indexes_io.begin()+ index_start, decompress_gt_indexes_io.begin() + id_pos[curr_non_copy_vec_id]);
    //     sparse_matrix_cols = vint_code::DecodeArray(curr_gt_indexes);
    //     int prev = -1;
    //     id_pos_size = sparse_matrix_cols.size();
    //     while (whichByte_whereInRes[next_haplotype].first < last_byte && next_haplotype < whichByte_whereInRes.back().second)
    //     {
    //         if (cur_i < id_pos_size)
    //         {
    //             if (cur_i != prev_i)
    //             {
    //                 cur_index = sparse_matrix_cols[cur_i]+prev;
    //                 prev = cur_index;
    //                 cur_pos = cur_index % 8;
    //                 cur_id = cur_index / 8;
    //             }
    //             prev_i = cur_i;
    //             if (cur_id < whichByte_whereInRes[next_haplotype].first)
    //                 cur_i += 1;
    //             else if (cur_id == whichByte_whereInRes[next_haplotype].first)
    //             {
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype].second] += perm_lut8[cur_pos];
    //                 cur_i += 1;
    //                 count++;
    //             }
    //             else
    //             {
    //                 next_haplotype++;
    //                 if (count)
    //                 {
    //                     if (next_haplotype < whichByte_whereInRes.back().second)
    //                         resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //                 }
    //                 else
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //             }
    //         }
    //         else
    //         {
    //             next_haplotype++;
    //             if (count)
    //             {
    //                 if (next_haplotype < whichByte_whereInRes.back().second)
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //             }
    //             else
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //         }
    //     }  
    // }
    return 0;
}

// // *****************************************************************************************************************
void Decompressor::calculate_start_position(int &_start_block,int &_start_position)
{
    pos_to_block.resize(fixed_variants_chunk_io.size(), 0);

    for (size_t i = 0; i < fixed_variants_chunk_io.size(); i++)
    {

        pos_to_block[i] = fixed_variants_chunk_io[i].data_compress[0].pos;
        
    }
    variant_desc_t cur_desc;
    _start_block = lower_bound(pos_to_block.begin(), pos_to_block.end(), range_1) - pos_to_block.begin() - 1;
    
    cur_desc.pos = range_1;

    if (_start_block < 0)
    {

        _start_block = 0;
        _start_position = 0;
    }
    else
    {
        if (cur_desc.pos > fixed_variants_chunk_io[_start_block].data_compress[fixed_variants_chunk_io[_start_block].data_compress.size() - 1].pos)
        {
            _start_position = 0;
            _start_block++;
        }
        else
        {
            _start_position = upper_bound(fixed_variants_chunk_io[_start_block].data_compress.begin(), 
                fixed_variants_chunk_io[_start_block].data_compress.end(), cur_desc, comp) - fixed_variants_chunk_io[_start_block].data_compress.begin();
            _start_position = _start_position > 0 ? _start_position : 0;
        }
    }    
    
}
void Decompressor::calculate_end_position(int &_end_block,int &_end_position){


    pos_to_block.resize(fixed_variants_chunk_io.size(), 0);

    for (size_t i = 0; i < fixed_variants_chunk_io.size(); i++)
    {

        pos_to_block[i] = fixed_variants_chunk_io[i].data_compress[0].pos;
    }
    variant_desc_t cur_desc;
    _end_block = upper_bound(pos_to_block.begin(), pos_to_block.end(), range_2) - pos_to_block.begin() - 1;
    
    cur_desc.pos = range_2;
    _end_position = upper_bound(fixed_variants_chunk_io[_end_block].data_compress.begin(), 
        fixed_variants_chunk_io[_end_block].data_compress.end(), cur_desc, comp1) - fixed_variants_chunk_io[_end_block].data_compress.begin();
    _end_position = _end_position > 0 ? _end_position : 0;

}
bool Decompressor::splitFileWriting(int file_num){

    auto logger = LogManager::Instance().Logger();
    split_files.resize(file_num);
    for(int i = 0; i < file_num; i++){
        char write_mode[5] = "wb-";
        if (out_type == file_type::BCF_File)
        {
            write_mode[3] = compression_level;
            write_mode[4] = '\0';
        }
        if (out_file_name != "-")
        {
            // char *gz_fname = (char *)malloc(strlen((out_file_name + decompression_reader.d_where_chrom[i].first).c_str()) + 5);

            if (out_type == file_type::VCF_File)
            {
                // snprintf(gz_fname, strlen((out_file_name + decompression_reader.d_where_chrom[i].first).c_str()) + 5, "%s.vcf", (out_file_name + decompression_reader.d_where_chrom[i].first).c_str());
                // split_files[i] = hts_open(gz_fname, "w");
                const std::string vcf_mode = vcf_write_mode(out_file_name, compression_level);
                split_files[i] = hts_open(out_file_name.c_str(), vcf_mode.c_str());
            }
            else
            {

                // snprintf(gz_fname, strlen((out_file_name + decompression_reader.d_where_chrom[i].first).c_str()) + 5, "%s.bcf", (out_file_name + decompression_reader.d_where_chrom[i].first).c_str());
                // split_files[i] = hts_open(gz_fname, write_mode);
                split_files[i] = hts_open(out_file_name.c_str(), write_mode);
            }

            // free(out_file_name);
            if (!split_files[i]){
                logger->error("Could not open split output file {} (index {})", out_file_name, i);
                return false;
            }
            else
            {
                hts_set_opt(split_files[i], HTS_OPT_CACHE_SIZE, 32000000);
                if (out_type == file_type::VCF_File) {
                    const std::string vcf_mode = vcf_write_mode(out_file_name, compression_level);
                    if (use_bgzf_threads(vcf_mode, params.no_threads)) {
                        hts_set_threads(split_files[i], params.no_threads);
                    }
                }
                rec = bcf_init();
            }
        }
    }
    return true;
}
// *****************************************************************************************************************
bool Decompressor::OpenForWriting()
{
    
    auto logger = LogManager::Instance().Logger();
    char write_mode[5] = "wb-";
    if (out_type == file_type::BCF_File)
    {
        write_mode[3] = compression_level;
        write_mode[4] = '\0';
    }


    if (out_file_name != "-")
    {
        if (out_type == file_type::VCF_File)
        {
            const std::string vcf_mode = vcf_write_mode(out_file_name, compression_level);
            out_file = hts_open(out_file_name.c_str(), vcf_mode.c_str());
            if (out_file && use_bgzf_threads(vcf_mode, params.no_threads))
                hts_set_threads(out_file, params.no_threads);
        }
        else if(out_type == file_type::BCF_File)
        {
            out_file = hts_open(out_file_name.c_str(), write_mode);
        }
        else if (out_type == file_type::BED_File)
        {
            if(!out_fam.Open(out_file_name + ".fam" , "w"))
	        {
                logger->error("Cannot open {}.fam", out_file_name);
		        return false;
	        }
            if(!out_bim.Open(out_file_name + ".bim" , "w"))
	        {
                logger->error("Cannot open {}.bim", out_file_name);
		        return false;
	        }
            if(!out_bed.Open(out_file_name + ".bed" , "wb"))
	        {
                logger->error("Cannot open {}.bed", out_file_name);
		        return false;
	        }
        }
        else if (out_type == file_type::PGEN_File)
        {
            if(!out_pgen.Open(out_file_name + ".pgen" , "wb"))
            {
                logger->error("Cannot open {}.pgen", out_file_name);
                return false;
            }
            if(!out_pvar.Open(out_file_name + ".pvar" , "w"))
            {
                logger->error("Cannot open {}.pvar", out_file_name);
                return false;
            }
            if(!out_psam.Open(out_file_name + ".psam" , "w"))
            {
                logger->error("Cannot open {}.psam", out_file_name);
                return false;
            }
        }
        else if (out_type == file_type::BGEN_File)
        {
            if(!out_bgen.Open(out_file_name + ".bgen" , "wb"))
            {
                logger->error("Cannot open {}.bgen", out_file_name);
                return false;
            }
            if(!out_sample.Open(out_file_name + ".sample" , "w"))
            {
                logger->error("Cannot open {}.sample", out_file_name);
                return false;
            }
        }
        else if (out_type == file_type::GDS_File)
        {
            logger->error("--make-gds is not implemented in the core binary; please convert via VCF + SeqArray.");
            return false;
        }
    }
    else
    {
        if (out_type == file_type::VCF_File)
        {
            const std::string vcf_mode = vcf_write_mode(out_file_name, compression_level);
            out_file = hts_open("-", vcf_mode.c_str());
            if (out_file && use_bgzf_threads(vcf_mode, params.no_threads))
                hts_set_threads(out_file, params.no_threads);
        }
        else if(out_type == file_type::BCF_File)
        {
            out_file = hts_open("-", write_mode);
        }
        else
        {
            logger->error("This output format requires --out <prefix> (cannot write to stdout).");
            return false;
        }
    }

    if(out_type == file_type::VCF_File || out_type == file_type::BCF_File){
        if (!out_file){
            logger->error("Could not open output file {}", out_file_name);
            return false;
        }
        else
        {
            hts_set_opt(out_file, HTS_OPT_CACHE_SIZE, 32000000);
            rec = bcf_init();
        }
    }
    return true;
}
// *****************************************************************************************************************
int Decompressor::analyzeInputSamples(vector<string> &v_samples)
{
    
    if (smpl.loadSamples(v_samples))
        return 1;
    
    return 0;
}
int Decompressor::initOutSplitFile(){

    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t";
    string str = "";
    smpl.get_all_samples(str);
    out_hdr = bcf_hdr_init("r");
    bcf_hdr_parse(out_hdr, (char *)header.c_str());
    bcf_hdr_sync(out_hdr);

    char delim = '\t';
    std::stringstream ss(str);
    std::string item;
    while (getline(ss, item, delim))
    {
        bcf_hdr_add_sample(out_hdr, item.c_str());
    }
    bcf_hdr_sync(out_hdr);
    
    for(size_t i = 0; i < split_files.size(); i++){
        if (bcf_hdr_write(split_files[i], out_hdr) < 0)
            return 1;
    }
    return 0;
}

void Decompressor::writeBgenHeaderAndSampleBlock()
{
    // BGEN v1.2, layout 2, zlib-compressed SNP blocks, with sample identifiers.
    const uint32_t sample_count = decompression_reader.n_samples;
    const uint32_t free_data_length = 0;
    const uint32_t header_length = 20 + free_data_length;
    const uint32_t flags = (1u << 0) | (2u << 2) | (1u << 31);

    // Sample identifier block size in bytes (including sample_count field).
    uint32_t sample_block_size = 8; // block_length (4) + sample_count (4)
    for (const auto& s : v_samples)
        sample_block_size += 2 + static_cast<uint32_t>(s.size());

    // First variant block absolute offset from file start.
    const uint32_t offset = 4 + 4 + header_length + sample_block_size;

    out_bgen.WriteUInt(offset, 4);
    out_bgen.WriteUInt(header_length, 4);
    out_bgen.WriteUInt(0, 4); // variant_count placeholder; patched on close
    out_bgen.WriteUInt(sample_count, 4);
    out_bgen.Write("bgen", 4);
    out_bgen.WriteUInt(free_data_length, 4);
    out_bgen.WriteUInt(flags, 4);

    // Sample identifier block
    uint32_t block_length = 4; // sample_count
    for (const auto& s : v_samples)
        block_length += 2 + static_cast<uint32_t>(s.size());
    out_bgen.WriteUInt(block_length, 4);
    out_bgen.WriteUInt(sample_count, 4);
    for (const auto& s : v_samples)
    {
        out_bgen.WriteUInt(static_cast<uint16_t>(s.size()), 2);
        out_bgen.Write(s.c_str(), s.size());
    }
}
// *****************************************************************************************************************
void Decompressor::WriteBEDMagicNumbers() {

    out_bed.PutByte(0x6C);  // 01101100
    out_bed.PutByte(0x1B);  // 00011011

    out_bed.PutByte(0x01);  // 00000001


    // PutByte(0x00);  // 00000000
}
int Decompressor::initOut()
{
    export_variant_count_ = 0;
    if(out_type == file_type::BED_File){
        string str = "";
        smpl.get_all_samples(str);
        char delim = '\t';
        std::stringstream ss(str);
        std::string item;
        string fam_line;
        while (getline(ss, item, delim))
        {
            fam_line = "0\t"+ item +"\t0\t0\t0\t-9\n";
            // std::cerr<< fam_line << endl;
            out_fam.Write(fam_line.c_str(), fam_line.size());
        }
        WriteBEDMagicNumbers();
    }
    else if (out_type == file_type::PGEN_File)
    {
        auto logger = LogManager::Instance().Logger();
        logger->warn("PGEN output is experimental; validate with plink2 before use.");

        // PGEN header (minimal placeholder; variant_count patched on close)
        out_pgen.PutByte(0x6C);
        out_pgen.PutByte(0x1B);
        out_pgen.PutByte(0x02);
        out_pgen.WriteUInt(0, 4); // variant_count placeholder
        out_pgen.WriteUInt(decompression_reader.n_samples, 4);

        out_pvar.Write("#CHROM\tPOS\tID\tREF\tALT\n", sizeof("#CHROM\tPOS\tID\tREF\tALT\n") - 1);
        out_psam.Write("#IID\n", 5);
        for (const auto& sample : v_samples)
        {
            std::string line = sample + "\n";
            out_psam.Write(line.c_str(), line.size());
        }
    }
    else if (out_type == file_type::BGEN_File)
    {
        auto logger = LogManager::Instance().Logger();
        logger->warn("BGEN output is experimental; validate with bgenix/qctool before use.");

        out_sample.Write("ID_1 ID_2 missing\n", sizeof("ID_1 ID_2 missing\n") - 1);
        out_sample.Write("0 0 0\n", 6);
        for (const auto& sample : v_samples)
        {
            std::string line = sample + " " + sample + " 0\n";
            out_sample.Write(line.c_str(), line.size());
        }

        writeBgenHeaderAndSampleBlock();
    }
    else{
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"; 
        string str = "";
        if (params.samples == "")
        {
            smpl.get_all_samples(str);
        }
        else
        {
            
            sampleIDs = smpl.setSamples(params.samples.c_str(), str);
            if (!sampleIDs)
            {
                auto logger = LogManager::Instance().Logger();
                logger->error("Failed to resolve sample subset: {}", params.samples);
                return 1;
            }
        }
           
        if (!out_file)
        {
            return 1;
        }
        out_hdr = bcf_hdr_init("r");
        if (out_genotypes)
            header += "FORMAT";
        bcf_hdr_parse(out_hdr, (char *)header.c_str());
        if(params.compress_mode == compress_mode_t::lossly_mode){
            bcf_hdr_remove(out_hdr,BCF_HL_INFO, NULL);
            bcf_hdr_remove(out_hdr,BCF_HL_FMT, NULL);
            bcf_hdr_append(out_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        }
        bcf_hdr_sync(out_hdr);
        if(params.out_AC_AN){
            bcf_hdr_append(out_hdr, "##INFO=<ID=AC,Number=A,Type=String,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
            bcf_hdr_append(out_hdr, "##INFO=<ID=AN,Number=A,Type=String,Description=\"Total number of alleles in called genotypes\">");
            bcf_hdr_sync(out_hdr);
        }

        if (out_genotypes)
        {
            char delim = '\t';
            std::stringstream ss(str);
            std::string item;
            while (getline(ss, item, delim))
            {
                if (out_genotypes)
                    bcf_hdr_add_sample(out_hdr, item.c_str());
            }
            bcf_hdr_sync(out_hdr);
        }

        const bool wants_no_header = !params.out_header_flag;
        if (wants_no_header && out_type == file_type::BCF_File)
        {
            LogManager::Instance().Logger()->warn("BCF output always includes header; ignoring --no-header");
        }

        const bool should_write_header =
            params.out_header_flag || params.header_only || (out_type == file_type::BCF_File);
        if (should_write_header)
        {
            if (bcf_hdr_write(out_file, out_hdr) < 0)
                return 1;
        }

        const char* low_copy_env = std::getenv("GSC_LOW_COPY_VCF");
        if (low_copy_env && *low_copy_env != '\0' && std::string(low_copy_env) != "0")
        {
            if (use_parallel_writing_)
                LogManager::Instance().Logger()->info("Low-copy VCF mode enabled; disabling parallel writer");
            use_parallel_writing_ = false;
        }

        // Initialize parallel writer if enabled
        if (use_parallel_writing_) {
            parallel_writer_ = new gsc::ParallelVCFWriter();

            // Use already-opened file handle (avoids duplicate header write and file conflicts)
            if (!parallel_writer_->InitializeWithHandle(out_file, out_hdr, true)) {
                LogManager::Instance().Logger()->error("Failed to initialize parallel VCF writer");
                delete parallel_writer_;
                parallel_writer_ = nullptr;
                use_parallel_writing_ = false;
                // Fall back to direct writing - out_file is still valid
            } else {
                LogManager::Instance().Logger()->info("Parallel VCF writer initialized successfully");
                // Transfer ownership: parallel writer will close the file in Finalize()
                out_file = nullptr;
            }
        }
    // }
    }
    return 0;
}
// *****************************************************************************************************************
// void Decompressor::setMemoryUsage()
// {
//     uint32_t no_haplotypes = smpl.no_samples * decompression_reader.ploidy;
//     // Set memory usage
//     if (params.MB_memory)
//     {
//         if (params.max_MB_memory)
//             max_stored_unique = (params.max_MB_memory * 1000000) / decompression_reader.vec_len;
//         else
//             max_stored_unique = decompression_reader.no_vec;
//         if (BLOCK_SIZE < max_stored_unique)
//             max_stored_unique = no_haplotypes*2;
//     }
//     else
//         max_stored_unique = 0;

// }

// // *****************************************************************************************************************
void Decompressor::decoded_vector_row(uint64_t vec_id, uint64_t offset, uint64_t _curr_non_copy_vec_id_offset, uint64_t length, int pos, uint8_t *decomp_data)
{
    auto logger = LogManager::Instance().Logger();

    uint64_t id, curr_non_copy_vec_id, toDelete;
    uint32_t tmp = 0;
    size_t id_pos_size;
    uint64_t index_start;
    int cur_id;
    int cur_pos;
    uint32_t cur_index;

    id = vec_id >> 1; // vec_id/2
    uint8_t parity = vec_id & 1;

    // Bounds check
    uint64_t bv_size = decompression_reader.rrr_zeros_bit_vector[parity].size();
    if (id >= bv_size) {
        logger->error("decoded_vector_row: vec_id={} id={} parity={} out of bounds! bv_size={}", vec_id, id, parity, bv_size);
        return;
    }

    if (decompression_reader.rrr_zeros_bit_vector[parity][id])
    {

        memcpy(decomp_data + pos, zeros_only_vector, length);

        // pos += decompression_reader.vec_len;
        
        return;
    }
    else if (decompression_reader.rrr_copy_bit_vector[parity][id])
    {

        unsigned long long bit_pos = (rrr_rank_copy_bit_vector[0](id + ((parity))) + rrr_rank_copy_bit_vector[1](id)) * decompression_reader.used_bits_cp;
        decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);      // /8
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);

        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
        // prev_chr_zeros_copy:Count the total number of all 1s and copied rows before the start position
        curr_non_copy_vec_id = vec_id - rrr_rank_zeros_bit_vector[0](id + ((parity))) - rrr_rank_zeros_bit_vector[1](id) -
                                rrr_rank_copy_bit_vector[0](id + ((parity))) - rrr_rank_copy_bit_vector[1](id) -_curr_non_copy_vec_id_offset;
        curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;

        got_it = done_unique.find(curr_non_copy_vec_id);
        
        if (got_it != done_unique.end())
        {
            memcpy(decomp_data + pos, got_it->second, length);
            // pos += decompression_reader.vec_len;
            return;
        }
    }
    else // unique and not a copy - no need to check if it was previously decompressed (got_it)
    {

        curr_non_copy_vec_id = vec_id - rrr_rank_zeros_bit_vector[0](id + ((parity))) - rrr_rank_zeros_bit_vector[1](id) -
                                rrr_rank_copy_bit_vector[0](id + ((parity))) - rrr_rank_copy_bit_vector[1](id) -_curr_non_copy_vec_id_offset;

        
    }
    if (vector_row_scratch_.size() < length)
        vector_row_scratch_.resize(length);
    uint8_t *cur_decomp_data = vector_row_scratch_.data();
    fill_n(cur_decomp_data, length, 0);

    // Bounds check for id_pos
    if (curr_non_copy_vec_id >= id_pos.size()) {
        logger->error("decoded_vector_row: curr_non_copy_vec_id={} out of bounds! id_pos.size()={}", curr_non_copy_vec_id, id_pos.size());
        return;
    }

    if (curr_non_copy_vec_id == 0)
        index_start = 0;
    else
        index_start = id_pos[curr_non_copy_vec_id - 1] + 1;

    // Bounds check for decompress_gt_indexes_io
    if (id_pos[curr_non_copy_vec_id] > decompress_gt_indexes_io.size()) {
        logger->error("decoded_vector_row: id_pos[{}]={} out of bounds! decompress_gt_indexes_io.size()={}",
                      curr_non_copy_vec_id, id_pos[curr_non_copy_vec_id], decompress_gt_indexes_io.size());
        return;
    }

    const size_t index_end = static_cast<size_t>(id_pos[curr_non_copy_vec_id]);
    if (index_end < index_start || index_end > decompress_gt_indexes_io.size())
    {
        logger->error("decoded_vector_row: invalid gt_index slice (start={}, end={}, size={})",
                      index_start, index_end, decompress_gt_indexes_io.size());
        return;
    }

    vint_code::DecodeArray(decompress_gt_indexes_io.data() + index_start,
                           index_end - index_start,
                           sparse_matrix_cols);
    id_pos_size = sparse_matrix_cols.size();
    int prev = -1;
    for(size_t i = 0; i < id_pos_size; ++i)
    {
        cur_index = sparse_matrix_cols[i]+prev;
        prev = cur_index;
        cur_pos = cur_index % 8;
        cur_id = cur_index / 8;
        cur_decomp_data[cur_id] += perm_lut8[cur_pos];
    }

    memcpy(decomp_data + pos, cur_decomp_data, length);
    // pos += decompression_reader.vec_len;

    if (max_stored_unique)
    {

        uint8_t *vector;
        if (done_unique.size() > max_stored_unique)
        {
            toDelete = stored_unique.back();
            stored_unique.pop_back();
            vector = done_unique[toDelete];
            done_unique.erase(toDelete);
        }
        else
        {
            vector = new uint8_t[length];
        }
        memcpy(vector, decomp_data + pos ,length);

        done_unique[curr_non_copy_vec_id] = vector;
        stored_unique.push_front(curr_non_copy_vec_id);

    }
}

bool Decompressor::decode_vector_row_partial_bytes(uint64_t vec_id,
                                                   uint64_t curr_non_copy_vec_id_offset,
                                                   const vector<ByteGroup> &groups,
                                                   const vector<pair<uint32_t, uint32_t>> &byte_where,
                                                   uint8_t *out_bytes,
                                                   uint32_t out_count)
{
    auto logger = LogManager::Instance().Logger();
    if (!out_bytes || !out_count)
        return true;

    fill_n(out_bytes, out_count, 0);
    if (groups.empty())
        return true;

    uint64_t id = vec_id >> 1;
    uint8_t parity = static_cast<uint8_t>(vec_id & 1);

    uint64_t bv_size = decompression_reader.rrr_zeros_bit_vector[parity].size();
    if (id >= bv_size)
    {
        logger->error("decode_vector_row_partial_bytes: vec_id={} id={} parity={} out of bounds (bv_size={})",
                      vec_id, id, parity, bv_size);
        return false;
    }

    if (decompression_reader.rrr_zeros_bit_vector[parity][id])
        return true;

    uint64_t non_copy_id = vec_id - rrr_rank_zeros_bit_vector[0](id + parity) - rrr_rank_zeros_bit_vector[1](id) -
                           rrr_rank_copy_bit_vector[0](id + parity) - rrr_rank_copy_bit_vector[1](id) - curr_non_copy_vec_id_offset;

    if (decompression_reader.rrr_copy_bit_vector[parity][id])
    {
        uint32_t tmp = 0;
        unsigned long long copy_rank = rrr_rank_copy_bit_vector[0](id + parity) + rrr_rank_copy_bit_vector[1](id);
        unsigned long long bit_pos = copy_rank * decompression_reader.used_bits_cp;
        decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7);
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);
        non_copy_id = non_copy_id - tmp - 1;
    }

    if (non_copy_id >= id_pos.size())
    {
        logger->error("decode_vector_row_partial_bytes: non_copy_id={} out of bounds (id_pos.size()={})",
                      non_copy_id, id_pos.size());
        return false;
    }

    size_t index_start = non_copy_id == 0 ? 0 : static_cast<size_t>(id_pos[non_copy_id - 1]) + 1;
    size_t index_end = static_cast<size_t>(id_pos[non_copy_id]);
    if (index_end > decompress_gt_indexes_io.size() || index_start > index_end)
    {
        logger->error("decode_vector_row_partial_bytes: invalid gt_index slice (start={}, end={}, size={})",
                      index_start, index_end, decompress_gt_indexes_io.size());
        return false;
    }

    if (index_end == index_start)
        return true;
    vint_code::DecodeArray(decompress_gt_indexes_io.data() + index_start,
                           index_end - index_start,
                           sparse_matrix_cols);

    const uint32_t last_byte = groups.back().byte_no;
    size_t group_idx = 0;

    int prev = -1;
    for (size_t i = 0; i < sparse_matrix_cols.size(); ++i)
    {
        uint32_t cur_index = sparse_matrix_cols[i] + prev;
        prev = static_cast<int>(cur_index);
        uint32_t cur_byte = cur_index >> 3;
        if (cur_byte > last_byte)
            break;

        while (group_idx < groups.size() && groups[group_idx].byte_no < cur_byte)
            ++group_idx;
        if (group_idx >= groups.size())
            break;
        if (groups[group_idx].byte_no != cur_byte)
            continue;

        const uint8_t bit = perm_lut8[cur_index & 7];
        for (uint32_t j = groups[group_idx].begin; j < groups[group_idx].end; ++j)
        {
            const uint32_t where = byte_where[j].second;
            if (where < out_count)
                out_bytes[where] += bit;
        }
    }

    return true;
}

// *****************************************************************************************************************
void Decompressor::decode_perm_rev(int vec2_start, const vector<uint32_t> &rev_perm, uint8_t *decomp_data_perm, uint8_t *decomp_data) 
{
    decode_perm_rev(vec2_start, decompression_reader.vec_len, standard_block_size, rev_perm, decomp_data_perm, decomp_data);
}

// *****************************************************************************************************************
void Decompressor::decode_perm_rev(int vec2_start, size_t vec_len, size_t block_haplotypes, const vector<uint32_t> &rev_perm,
                                   uint8_t *decomp_data_perm, uint8_t *decomp_data)
{

    for (size_t i = 0; i < vec_len * 2; i++)
        decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
    size_t x;
    for (x = 0; x + 8 < block_haplotypes;)
    {
        int x8 = x / 8;
        int d_x8 = decomp_data_perm[x8];
        int d2_x8 = decomp_data_perm[vec2_start + x8];
        if (!d_x8 && !d2_x8)
        {
            x += 8;
            continue;
        }

        int x_p = 1 << 7;
        for (int i = 0; i < 8; ++i, ++x, x_p >>= 1)
        {
            auto j = rev_perm[x];
            uint8_t j_p = perm_lut8[j % 8];
            int j8 = j / 8;
            if (d_x8 & x_p)
                decomp_data[j8] += j_p;
            if (d2_x8 & x_p)
                decomp_data[vec2_start + j8] += j_p;
        }
    }
    int x8 = x / 8;
    uint8_t d_x8 = decomp_data_perm[x8];
    uint8_t d2_x8 = decomp_data_perm[vec2_start + x8];

    for (; x < block_haplotypes; ++x)
    {
        auto j = rev_perm[x];
        uint8_t x_p = perm_lut8[x % 8];
        uint8_t j_p = perm_lut8[j % 8];
        int j8 = j / 8;
        if (d_x8 & x_p)
            decomp_data[j8] += j_p;
        if (d2_x8 & x_p)
            decomp_data[vec2_start + j8] += j_p;
    }

}
// *****************************************************************************************************************
void inline Decompressor::reverse_perm(const vector<uint32_t> &perm, vector<uint32_t> &rev_perm, int no_haplotypes)
{
    for (int i = 0; i < no_haplotypes; ++i)
        rev_perm[perm[i]] = i;
}
