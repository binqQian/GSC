/**
 * @file gvcf_block.h
 * @brief gVCF block structure definitions
 *
 * Defines data structures for gVCF-specific fields and blocks.
 * Reference: ref_code/gvcf_refCode/io/record.go
 */
#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <memory>

namespace gvcf {

// ============================================================================
// Genotype Data (GT field)
// ============================================================================

struct GenotypeData {
    int8_t allele1;   // First allele (-1 for missing)
    int8_t allele2;   // Second allele (-1 for missing)
    bool phased;      // True if phased (|), false if unphased (/)

    GenotypeData() : allele1(-1), allele2(-1), phased(false) {}
    GenotypeData(int8_t a1, int8_t a2, bool p = false)
        : allele1(a1), allele2(a2), phased(p) {}

    bool operator==(const GenotypeData& other) const {
        return allele1 == other.allele1 &&
               allele2 == other.allele2 &&
               phased == other.phased;
    }

    bool IsMissing() const {
        return allele1 < 0 && allele2 < 0;
    }

    bool IsHomRef() const {
        return allele1 == 0 && allele2 == 0;
    }

    bool IsHet() const {
        return allele1 != allele2 && allele1 >= 0 && allele2 >= 0;
    }

    // Convert to string representation (e.g., "0/0", "0|1", "./.")
    std::string ToString() const;

    // Parse from string
    static GenotypeData FromString(const std::string& str);
};

// ============================================================================
// gVCF Known Fields Block
// ============================================================================

struct GVCFPositionFields {
    std::vector<std::string> chrom;     // CHROM
    std::vector<uint64_t> pos;          // POS
    std::vector<std::string> id;        // ID (typically ".")
};

struct GVCFSequenceFields {
    std::vector<std::string> ref;                    // REF
    std::vector<std::vector<std::string>> alt;       // ALT (multi-allelic)
};

struct GVCFQualityFields {
    std::vector<float> qual;            // QUAL
    std::vector<std::string> filter;    // FILTER
};

struct GVCFInfoFields {
    std::vector<int64_t> end;           // INFO/END (gVCF block end position)
    // Other INFO fields can be added here
};

struct GVCFSampleFields {
    std::vector<GenotypeData> gt;       // GT
    std::vector<int32_t> dp;            // DP (read depth)
    std::vector<int32_t> gq;            // GQ (genotype quality)
    std::vector<int32_t> min_dp;        // MIN_DP (minimum depth in gVCF block)
    std::vector<std::vector<int32_t>> pl; // PL (phred-scaled likelihoods)
    std::vector<std::vector<int32_t>> ad; // AD (allelic depths)
};

// ============================================================================
// gVCF Block
// ============================================================================

struct GVCFBlock {
    // Known fields
    GVCFPositionFields position;
    GVCFSequenceFields sequence;
    GVCFQualityFields quality;
    GVCFInfoFields info;
    GVCFSampleFields sample;

    // Unknown/generic fields (INFO and FORMAT)
    std::unordered_map<std::string, std::vector<std::string>> unknown_info;
    std::unordered_map<std::string, std::vector<std::string>> unknown_format;

    // Block metadata
    uint32_t variant_count;
    uint32_t sample_count;

    GVCFBlock() : variant_count(0), sample_count(0) {}

    void Clear() {
        position.chrom.clear();
        position.pos.clear();
        position.id.clear();
        sequence.ref.clear();
        sequence.alt.clear();
        quality.qual.clear();
        quality.filter.clear();
        info.end.clear();
        sample.gt.clear();
        sample.dp.clear();
        sample.gq.clear();
        sample.min_dp.clear();
        sample.pl.clear();
        sample.ad.clear();
        unknown_info.clear();
        unknown_format.clear();
        variant_count = 0;
        sample_count = 0;
    }

    void Reserve(uint32_t n_variants) {
        position.chrom.reserve(n_variants);
        position.pos.reserve(n_variants);
        position.id.reserve(n_variants);
        sequence.ref.reserve(n_variants);
        sequence.alt.reserve(n_variants);
        quality.qual.reserve(n_variants);
        quality.filter.reserve(n_variants);
        info.end.reserve(n_variants);
    }
};

// ============================================================================
// Compressed Field Storage
// ============================================================================

enum class FieldCompressionMethod : uint8_t {
    NONE = 0,
    RLE = 1,
    DELTA = 2,
    MASK = 3,
    DICTIONARY = 4,
    DELTA_RLE = 5,      // Delta + RLE for indices
    MASK_RLE = 6,       // Mask with RLE-compressed bitmask
    ADAPTIVE = 7        // Auto-select best method
};

struct CompressedField {
    FieldCompressionMethod method;
    std::vector<uint8_t> data;
    uint32_t original_count;

    // Statistics (for debugging/analysis)
    float compression_ratio;
    float dominant_ratio;  // For mask encoding

    CompressedField()
        : method(FieldCompressionMethod::NONE)
        , original_count(0)
        , compression_ratio(0.0f)
        , dominant_ratio(0.0f) {}

    void Clear() {
        method = FieldCompressionMethod::NONE;
        data.clear();
        original_count = 0;
        compression_ratio = 0.0f;
        dominant_ratio = 0.0f;
    }
};

// ============================================================================
// Compressed gVCF Block
// ============================================================================

struct CompressedGVCFBlock {
    // Position fields
    CompressedField chrom;
    CompressedField pos;
    CompressedField id;

    // Sequence fields
    CompressedField ref;
    CompressedField alt;

    // Quality fields
    CompressedField qual;
    CompressedField filter;

    // INFO fields
    CompressedField info_end;

    // Sample fields (gVCF optimized)
    CompressedField gt_mask;    // Bitmask for dominant genotype
    CompressedField gt_patches; // Non-dominant genotypes
    CompressedField gt_phase;   // Phase information

    CompressedField dp;
    CompressedField gq;
    CompressedField min_dp;
    CompressedField dp_min_dp_diff; // DP - MIN_DP differences

    CompressedField pl;
    CompressedField ad;

    // Unknown fields
    std::unordered_map<std::string, CompressedField> unknown_info;
    std::unordered_map<std::string, CompressedField> unknown_format;

    // Metadata
    uint32_t variant_count;
    uint32_t sample_count;

    // gVCF specific flags
    bool has_end_field;
    bool has_min_dp;

    CompressedGVCFBlock()
        : variant_count(0)
        , sample_count(0)
        , has_end_field(false)
        , has_min_dp(false) {}

    void Clear() {
        chrom.Clear();
        pos.Clear();
        id.Clear();
        ref.Clear();
        alt.Clear();
        qual.Clear();
        filter.Clear();
        info_end.Clear();
        gt_mask.Clear();
        gt_patches.Clear();
        gt_phase.Clear();
        dp.Clear();
        gq.Clear();
        min_dp.Clear();
        dp_min_dp_diff.Clear();
        pl.Clear();
        ad.Clear();
        unknown_info.clear();
        unknown_format.clear();
        variant_count = 0;
        sample_count = 0;
        has_end_field = false;
        has_min_dp = false;
    }

    // Calculate total compressed size
    size_t TotalCompressedSize() const;

    // Serialize to byte buffer
    bool Serialize(std::vector<uint8_t>& buffer) const;

    // Deserialize from byte buffer
    bool Deserialize(const std::vector<uint8_t>& buffer);
    bool Deserialize(const uint8_t* buffer, size_t size);
};

// ============================================================================
// gVCF Block Configuration
// ============================================================================

struct GVCFBlockConfig {
    uint32_t block_size;           // Number of variants per block
    bool enable_adaptive;          // Auto-select compression method
    float mask_threshold;          // Min dominant ratio for mask encoding (default 0.7)

    // gVCF-specific options
    std::string default_gt;        // Expected dominant genotype (default "0/0")
    std::string default_alt;       // Expected dominant ALT (default "<NON_REF>")

    GVCFBlockConfig()
        : block_size(10000)
        , enable_adaptive(true)
        , mask_threshold(0.7f)
        , default_gt("0/0")
        , default_alt("<NON_REF>") {}
};

// ============================================================================
// gVCF Analysis Results
// ============================================================================

struct GVCFFieldAnalysis {
    std::string field_name;
    uint32_t total_count;
    uint32_t unique_count;
    float repetition_ratio;     // unique / total
    float dominant_ratio;       // Most frequent / total

    // For numeric fields
    int64_t min_value;
    int64_t max_value;
    double avg_value;

    // Recommended compression method
    FieldCompressionMethod recommended_method;
};

struct GVCFBlockAnalysis {
    uint32_t variant_count;
    uint32_t sample_count;

    // Per-field analysis
    GVCFFieldAnalysis chrom_analysis;
    GVCFFieldAnalysis pos_analysis;
    GVCFFieldAnalysis id_analysis;
    GVCFFieldAnalysis ref_analysis;
    GVCFFieldAnalysis alt_analysis;
    GVCFFieldAnalysis gt_analysis;
    GVCFFieldAnalysis dp_analysis;
    GVCFFieldAnalysis gq_analysis;
    GVCFFieldAnalysis min_dp_analysis;
    GVCFFieldAnalysis end_analysis;

    // Overall statistics
    float estimated_compression_ratio;
    bool is_typical_gvcf;  // Has expected gVCF characteristics
};

} // namespace gvcf
