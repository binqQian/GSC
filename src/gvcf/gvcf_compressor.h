/**
 * @file gvcf_compressor.h
 * @brief gVCF compression entry point
 *
 * Integrates gVCF compression module with the main GSC workflow.
 * Handles single-sample gVCF files using optimized compression strategies.
 */
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include "gvcf.h"
#include "gsc_params.h"

namespace gvcf {

// Magic number for gVCF compressed format (different from regular GSC)
constexpr uint32_t GVCF_FILE_MAGIC = 0x47564346; // "GVCF"
constexpr uint32_t GVCF_FILE_VERSION = 4;  // V4: block index for range query

// Backend selection functions
void SetGVCFBackend(compression_backend_t backend);
compression_backend_t GetGVCFBackend();

/**
 * Block index for range query support.
 * Stores the genomic range covered by each block.
 */
struct BlockIndex {
    std::string chrom;        // Chromosome name
    uint64_t start_pos;       // Start position (1-based, inclusive)
    uint64_t end_pos;         // End position (1-based, inclusive, using END field if available)
    uint64_t file_offset;     // File offset to the block
    uint32_t variant_count;   // Number of variants in the block

    BlockIndex() : start_pos(0), end_pos(0), file_offset(0), variant_count(0) {}
    BlockIndex(const std::string& c, uint64_t s, uint64_t e, uint64_t off, uint32_t cnt)
        : chrom(c), start_pos(s), end_pos(e), file_offset(off), variant_count(cnt) {}
};

/**
 * gVCF Compressor - handles single-sample gVCF file compression
 */
class GVCFCompressor {
public:
    explicit GVCFCompressor(const GSC_Params& params);
    ~GVCFCompressor();

    // Main compression entry point
    bool Compress();

    // Get compression statistics
    struct Statistics {
        uint64_t total_variants;
        uint64_t total_blocks;
        uint64_t original_size;
        uint64_t compressed_size;
        float compression_ratio;

        // Per-field statistics
        uint64_t chrom_compressed;
        uint64_t pos_compressed;
        uint64_t gt_compressed;
        uint64_t dp_compressed;
        uint64_t gq_compressed;
        uint64_t pl_compressed;
    };

    const Statistics& GetStatistics() const { return stats_; }

private:
    // VCF reading using htslib
    bool OpenInput();
    bool ReadHeader();
    bool ReadBlock(GVCFBlock& block);
    void CloseInput();

    // Compression
    bool CompressAndWriteBlock(const GVCFBlock& block);
    bool WriteFileHeader();
    bool WriteFileFooter();

    // gVCF detection
    bool DetectGVCFFormat();

    // Backend setup
    std::shared_ptr<CompressionBackend> CreateBackend();

    GSC_Params params_;
    Statistics stats_;

    // htslib handles
    void* vcf_file_;      // htsFile*
    void* vcf_header_;    // bcf_hdr_t*
    void* vcf_record_;    // bcf1_t*

    // Output file
    FILE* output_file_;
    std::string output_filename_;

    // Compression context
    std::shared_ptr<CompressionBackend> backend_;
    GVCFBlockConfig config_;

    // Block tracking
    uint32_t current_block_id_;
    std::vector<uint64_t> block_offsets_;
    std::vector<BlockIndex> block_indices_;  // V4: block index for range query

    // Sample info
    uint32_t num_samples_;
    std::vector<std::string> sample_names_;
    std::string header_text_;

    // gVCF detection flags
    bool is_gvcf_;
    bool has_end_field_;
    bool has_min_dp_;
    bool has_non_ref_;
};

/**
 * gVCF Decompressor - handles gVCF compressed file decompression
 */
class GVCFDecompressor {
public:
    explicit GVCFDecompressor(const GSC_Params& params);
    ~GVCFDecompressor();

    // Main decompression entry point
    bool Decompress();

    // Query support (deprecated, use GVCFQueryer instead)
    bool DecompressRange(const std::string& chrom, uint64_t start, uint64_t end);

private:
    bool OpenInput();
    bool ReadFileHeader();
    bool ReadBlock(uint32_t block_id, GVCFBlock& block);
    bool WriteVCFRecord(const GVCFBlock& block, uint32_t idx);
    void CloseInput();

    GSC_Params params_;

    // Input file
    FILE* input_file_;

    // Output file
    void* output_file_;   // htsFile*
    void* output_header_; // bcf_hdr_t*

    // Decompression context
    std::shared_ptr<CompressionBackend> backend_;

    // File metadata
    uint32_t num_blocks_;
    uint32_t file_version_;  // V4: to handle version differences
    std::vector<uint64_t> block_offsets_;
    std::vector<BlockIndex> block_indices_;  // V4: block index for range query
    std::string header_text_;
    std::vector<std::string> sample_names_;
};

/**
 * gVCF Range Query - supports querying variants by genomic range
 */
class GVCFQueryer {
public:
    explicit GVCFQueryer(const std::string& input_file);
    ~GVCFQueryer();

    // Open file and read index
    bool Open();

    // Query variants in a range, output to file
    bool QueryRange(const std::string& chrom, uint64_t start, uint64_t end,
                   const std::string& output_file);

    // Get block indices (for inspection)
    const std::vector<BlockIndex>& GetBlockIndices() const { return block_indices_; }

    // Get total variant count
    uint64_t GetTotalVariants() const { return total_variants_; }

private:
    // Find blocks that overlap with the query range
    std::vector<uint32_t> FindBlocksInRange(const std::string& chrom,
                                            uint64_t start, uint64_t end) const;

    // Read and decompress a specific block
    bool ReadBlock(uint32_t block_id, GVCFBlock& block);

    // Write VCF record
    bool WriteVCFRecord(void* output_file, void* header,
                       const GVCFBlock& block, uint32_t idx);

    std::string input_filename_;
    FILE* input_file_;

    // Decompression context
    std::shared_ptr<CompressionBackend> backend_;

    // File metadata
    uint32_t num_blocks_;
    uint32_t file_version_;
    uint64_t total_variants_;
    std::vector<uint64_t> block_offsets_;
    std::vector<BlockIndex> block_indices_;
    std::string header_text_;
    std::vector<std::string> sample_names_;
};

/**
 * Check if a VCF file appears to be gVCF format.
 * Returns true if the file has typical gVCF characteristics.
 */
bool IsGVCFFile(const std::string& filename);

/**
 * Check if a compressed file is in gVCF format.
 */
bool IsGVCFCompressed(const std::string& filename);

} // namespace gvcf
