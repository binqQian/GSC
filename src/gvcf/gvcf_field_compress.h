/**
 * @file gvcf_field_compress.h
 * @brief gVCF field-level compression
 *
 * Provides specialized compression for each gVCF field type.
 * Reference: ref_code/gvcf_refCode/compress/
 */
#pragma once

#include "gvcf_encoding.h"
#include "gvcf_block.h"
#include <functional>

namespace gvcf {

// Forward declaration for backend compression
class CompressionBackend;

// ============================================================================
// Compression Backend Interface
// ============================================================================

class CompressionBackend {
public:
    virtual ~CompressionBackend() = default;
    virtual bool Compress(const std::vector<uint8_t>& input,
                         std::vector<uint8_t>& output) = 0;
    virtual bool Decompress(const std::vector<uint8_t>& input,
                           std::vector<uint8_t>& output) = 0;
    virtual bool Decompress(const uint8_t* input, size_t input_size,
                           std::vector<uint8_t>& output) = 0;
};

// ============================================================================
// Field Compressor Base
// ============================================================================

class FieldCompressor {
protected:
    std::shared_ptr<CompressionBackend> backend_;
    GVCFBlockConfig config_;

public:
    FieldCompressor(std::shared_ptr<CompressionBackend> backend,
                   const GVCFBlockConfig& config)
        : backend_(backend), config_(config) {}

    virtual ~FieldCompressor() = default;

    // Apply backend compression to encoded data
    bool ApplyBackendCompression(const std::vector<uint8_t>& input,
                                 CompressedField& output);
};

// ============================================================================
// CHROM Field Compressor
// Strategy: RLE + Backend compression
// ============================================================================

class ChromCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::string>& data, CompressedField& output);

    // Analysis
    static GVCFFieldAnalysis Analyze(const std::vector<std::string>& data);
};

// ============================================================================
// POS Field Compressor
// Strategy: Delta encoding + Backend compression
// ============================================================================

class PosCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<uint64_t>& data, CompressedField& output);

    // Analysis
    static GVCFFieldAnalysis Analyze(const std::vector<uint64_t>& data);
};

// ============================================================================
// INFO/END Field Compressor
// Strategy: Inferred from next POS (99.97% can be inferred)
// Only stores: exception bitmap + exception values + last record's END
// ============================================================================

class EndCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    // Compress END using inference from next POS
    // For most records: END = next_POS - 1 (no storage needed)
    // Only stores exceptions where END != next_POS - 1
    bool Compress(const std::vector<int64_t>& end_values,
                 const std::vector<uint64_t>& pos_values,
                 CompressedField& output);

    // Compress END standalone (fallback, less efficient)
    bool CompressStandalone(const std::vector<int64_t>& data, CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<int64_t>& data);
};

// ============================================================================
// ID Field Compressor
// Strategy: Mask encoding (dominant ".") + Backend compression
// ============================================================================

class IdCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::string>& data, CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<std::string>& data);
};

// ============================================================================
// REF Field Compressor
// Strategy: Adaptive (Dictionary or RLE) + Backend compression
// ============================================================================

class RefCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::string>& data, CompressedField& output);

    // Select compression method based on data characteristics
    static FieldCompressionMethod SelectMethod(const std::vector<std::string>& data);
    static GVCFFieldAnalysis Analyze(const std::vector<std::string>& data);
};

// ============================================================================
// ALT Field Compressor
// Strategy: Mask encoding (dominant "<NON_REF>") + Dictionary + Backend
// ============================================================================

class AltCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::vector<std::string>>& data,
                 CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<std::vector<std::string>>& data);
};

// ============================================================================
// QUAL Field Compressor
// Strategy: Quantization + Backend compression
// ============================================================================

class QualCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<float>& data, CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<float>& data);
};

// ============================================================================
// FILTER Field Compressor
// Strategy: RLE or Mask + Backend compression
// ============================================================================

class FilterCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::string>& data, CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<std::string>& data);
};

// ============================================================================
// GT Field Compressor (gVCF optimized)
// Strategy: Mask (dominant 0/0) + RLE bitmask + Patches + Phase
// ============================================================================

class GTCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    // Compress GT data, outputs multiple compressed fields
    bool Compress(const std::vector<GenotypeData>& data,
                 CompressedField& mask_output,
                 CompressedField& patches_output,
                 CompressedField& phase_output);

    // Simplified compression to single field
    bool CompressSingle(const std::vector<GenotypeData>& data,
                       CompressedField& output);

    // Compress from string representation
    bool CompressStrings(const std::vector<std::string>& data,
                        CompressedField& mask_output,
                        CompressedField& patches_output,
                        CompressedField& phase_output);

    static GVCFFieldAnalysis Analyze(const std::vector<GenotypeData>& data);
    static GVCFFieldAnalysis AnalyzeStrings(const std::vector<std::string>& data);
};

// ============================================================================
// DP Field Compressor
// Strategy: Adaptive encoding + Backend compression
// ============================================================================

class DPCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<int32_t>& data, CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<int32_t>& data);
};

// ============================================================================
// MIN_DP Field Compressor (gVCF specific)
// Strategy: Encode differences from DP + Adaptive
// ============================================================================

class MinDPCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    // Compress MIN_DP relative to DP (stores DP - MIN_DP)
    bool Compress(const std::vector<int32_t>& min_dp,
                 const std::vector<int32_t>& dp,
                 CompressedField& diff_output,
                 CompressedField& dp_output);

    // Compress MIN_DP standalone
    bool CompressStandalone(const std::vector<int32_t>& data,
                           CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<int32_t>& min_dp,
                                     const std::vector<int32_t>& dp);
};

// ============================================================================
// GQ Field Compressor
// Strategy: Predict from PL + Store exceptions
// ============================================================================

class GQCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    // Compress GQ with PL prediction
    bool Compress(const std::vector<int32_t>& gq,
                 const std::vector<std::vector<int32_t>>& pl,
                 CompressedField& output);

    // Compress GQ standalone
    bool CompressStandalone(const std::vector<int32_t>& data,
                           CompressedField& output);

    // Predict GQ from PL values
    static int32_t PredictGQFromPL(const std::vector<int32_t>& pl);

    static GVCFFieldAnalysis Analyze(const std::vector<int32_t>& data);
};

// ============================================================================
// PL Field Compressor
// Strategy: Pattern dictionary + Adaptive encoding
// ============================================================================

class PLCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::vector<int32_t>>& data,
                 CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<std::vector<int32_t>>& data);

private:
    // Pattern encoding types
    enum class PLPatternType : uint8_t {
        PATTERN_REF = 0,     // Pattern index reference
        DELTA_ENCODED = 1,   // Delta from first value
        RLE_ENCODED = 2,     // RLE encoded
        RAW = 3              // Raw values
    };
};

// ============================================================================
// AD Field Compressor
// Strategy: Adaptive encoding based on allele count
// ============================================================================

class ADCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::vector<int32_t>>& data,
                 CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<std::vector<int32_t>>& data);
};

// ============================================================================
// Generic/Unknown Field Compressor
// Strategy: Dictionary + Backend compression
// ============================================================================

class GenericFieldCompressor : public FieldCompressor {
public:
    using FieldCompressor::FieldCompressor;

    bool Compress(const std::vector<std::string>& data, CompressedField& output);

    static GVCFFieldAnalysis Analyze(const std::vector<std::string>& data);
};

// ============================================================================
// gVCF Block Compressor
// Orchestrates compression of all fields
// ============================================================================

class GVCFBlockCompressor {
public:
    GVCFBlockCompressor(std::shared_ptr<CompressionBackend> backend,
                       const GVCFBlockConfig& config = GVCFBlockConfig());

    // Compress entire block
    bool Compress(const GVCFBlock& input, CompressedGVCFBlock& output);

    // Analyze block for compression characteristics
    GVCFBlockAnalysis Analyze(const GVCFBlock& block);

    // Get individual field compressors for custom use
    ChromCompressor& GetChromCompressor() { return chrom_comp_; }
    PosCompressor& GetPosCompressor() { return pos_comp_; }
    GTCompressor& GetGTCompressor() { return gt_comp_; }

private:
    std::shared_ptr<CompressionBackend> backend_;
    GVCFBlockConfig config_;

    // Field compressors
    ChromCompressor chrom_comp_;
    PosCompressor pos_comp_;
    EndCompressor end_comp_;
    IdCompressor id_comp_;
    RefCompressor ref_comp_;
    AltCompressor alt_comp_;
    QualCompressor qual_comp_;
    FilterCompressor filter_comp_;
    GTCompressor gt_comp_;
    DPCompressor dp_comp_;
    MinDPCompressor min_dp_comp_;
    GQCompressor gq_comp_;
    PLCompressor pl_comp_;
    ADCompressor ad_comp_;
    GenericFieldCompressor generic_comp_;
};

} // namespace gvcf
