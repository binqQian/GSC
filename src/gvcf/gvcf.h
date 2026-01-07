/**
 * @file gvcf.h
 * @brief Main header for gVCF compression module
 *
 * This file provides the public interface for gVCF compression.
 * Include this header to use gVCF compression features.
 */
#pragma once

// Core modules
#include "gvcf_encoding.h"
#include "gvcf_block.h"
#include "gvcf_field_compress.h"
#include "gvcf_field_decompress.h"

namespace gvcf {

// ============================================================================
// Version Information
// ============================================================================

constexpr uint32_t GVCF_MODULE_VERSION_MAJOR = 1;
constexpr uint32_t GVCF_MODULE_VERSION_MINOR = 0;
constexpr uint32_t GVCF_MODULE_VERSION_PATCH = 0;

inline std::string GetVersion() {
    return std::to_string(GVCF_MODULE_VERSION_MAJOR) + "." +
           std::to_string(GVCF_MODULE_VERSION_MINOR) + "." +
           std::to_string(GVCF_MODULE_VERSION_PATCH);
}

// ============================================================================
// Magic Numbers and File Format
// ============================================================================

constexpr uint32_t GVCF_MAGIC = 0x47564346; // "GVCF"
constexpr uint32_t GVCF_FORMAT_VERSION = 1;

// ============================================================================
// Backend Adapter for GSC Integration
// ============================================================================

/**
 * Adapter to connect gVCF compression to GSC's compression backends.
 * This class wraps GSC's CompressionStrategy interface.
 */
class GSCBackendAdapter : public CompressionBackend {
public:
    // Compression function signature
    using CompressFunc = std::function<bool(const std::vector<uint8_t>&,
                                            std::vector<uint8_t>&)>;
    using DecompressFunc = std::function<bool(const std::vector<uint8_t>&,
                                              std::vector<uint8_t>&)>;
    using DecompressPtrFunc = std::function<bool(const uint8_t*, size_t,
                                                 std::vector<uint8_t>&)>;

    GSCBackendAdapter(CompressFunc compress_fn,
                     DecompressFunc decompress_fn,
                     DecompressPtrFunc decompress_ptr_fn = nullptr)
        : compress_fn_(std::move(compress_fn))
        , decompress_fn_(std::move(decompress_fn))
        , decompress_ptr_fn_(std::move(decompress_ptr_fn)) {}

    bool Compress(const std::vector<uint8_t>& input,
                 std::vector<uint8_t>& output) override {
        if (compress_fn_) {
            return compress_fn_(input, output);
        }
        output = input;
        return true;
    }

    bool Decompress(const std::vector<uint8_t>& input,
                   std::vector<uint8_t>& output) override {
        if (decompress_fn_) {
            return decompress_fn_(input, output);
        }
        output = input;
        return true;
    }

    bool Decompress(const uint8_t* input, size_t input_size,
                   std::vector<uint8_t>& output) override {
        if (decompress_ptr_fn_) {
            return decompress_ptr_fn_(input, input_size, output);
        }
        // Fallback to vector-based decompress
        std::vector<uint8_t> tmp(input, input + input_size);
        return Decompress(tmp, output);
    }

private:
    CompressFunc compress_fn_;
    DecompressFunc decompress_fn_;
    DecompressPtrFunc decompress_ptr_fn_;
};

/**
 * No-op backend that doesn't apply any additional compression.
 * Useful for testing or when backend compression is handled elsewhere.
 */
class NoopBackend : public CompressionBackend {
public:
    bool Compress(const std::vector<uint8_t>& input,
                 std::vector<uint8_t>& output) override {
        output = input;
        return true;
    }

    bool Decompress(const std::vector<uint8_t>& input,
                   std::vector<uint8_t>& output) override {
        output = input;
        return true;
    }

    bool Decompress(const uint8_t* input, size_t input_size,
                   std::vector<uint8_t>& output) override {
        output.assign(input, input + input_size);
        return true;
    }
};

// ============================================================================
// High-Level API
// ============================================================================

/**
 * Create a compression backend from function pointers.
 */
inline std::shared_ptr<CompressionBackend> CreateBackend(
    GSCBackendAdapter::CompressFunc compress_fn,
    GSCBackendAdapter::DecompressFunc decompress_fn,
    GSCBackendAdapter::DecompressPtrFunc decompress_ptr_fn = nullptr) {
    return std::make_shared<GSCBackendAdapter>(
        std::move(compress_fn),
        std::move(decompress_fn),
        std::move(decompress_ptr_fn));
}

/**
 * Create a no-op backend.
 */
inline std::shared_ptr<CompressionBackend> CreateNoopBackend() {
    return std::make_shared<NoopBackend>();
}

/**
 * Quick compression of GT field (most common use case).
 */
inline bool CompressGT(const std::vector<std::string>& gt_strings,
                      std::shared_ptr<CompressionBackend> backend,
                      CompressedField& mask_output,
                      CompressedField& patches_output,
                      CompressedField& phase_output,
                      const GVCFBlockConfig& config = GVCFBlockConfig()) {
    GTCompressor compressor(backend, config);
    return compressor.CompressStrings(gt_strings, mask_output, patches_output, phase_output);
}

/**
 * Quick decompression of GT field.
 */
inline bool DecompressGT(const CompressedField& mask_input,
                        const CompressedField& patches_input,
                        const CompressedField& phase_input,
                        std::shared_ptr<CompressionBackend> backend,
                        std::vector<std::string>& output) {
    GTDecompressor decompressor(backend);
    return decompressor.DecompressStrings(mask_input, patches_input, phase_input, output);
}

/**
 * Analyze gVCF block to get compression statistics.
 */
inline GVCFBlockAnalysis AnalyzeBlock(const GVCFBlock& block,
                                      std::shared_ptr<CompressionBackend> backend = nullptr,
                                      const GVCFBlockConfig& config = GVCFBlockConfig()) {
    if (!backend) {
        backend = CreateNoopBackend();
    }
    GVCFBlockCompressor compressor(backend, config);
    return compressor.Analyze(block);
}

/**
 * Compress entire gVCF block.
 */
inline bool CompressBlock(const GVCFBlock& input,
                         CompressedGVCFBlock& output,
                         std::shared_ptr<CompressionBackend> backend = nullptr,
                         const GVCFBlockConfig& config = GVCFBlockConfig()) {
    if (!backend) {
        backend = CreateNoopBackend();
    }
    GVCFBlockCompressor compressor(backend, config);
    return compressor.Compress(input, output);
}

/**
 * Decompress entire gVCF block.
 */
inline bool DecompressBlock(const CompressedGVCFBlock& input,
                           GVCFBlock& output,
                           std::shared_ptr<CompressionBackend> backend = nullptr) {
    if (!backend) {
        backend = CreateNoopBackend();
    }
    GVCFBlockDecompressor decompressor(backend);
    return decompressor.Decompress(input, output);
}

/**
 * Validate compression by comparing original and decompressed blocks.
 */
inline bool ValidateCompression(const GVCFBlock& original,
                               std::shared_ptr<CompressionBackend> backend = nullptr,
                               const GVCFBlockConfig& config = GVCFBlockConfig()) {
    if (!backend) {
        backend = CreateNoopBackend();
    }

    // Compress
    CompressedGVCFBlock compressed;
    if (!CompressBlock(original, compressed, backend, config)) {
        return false;
    }

    // Decompress
    GVCFBlock decompressed;
    if (!DecompressBlock(compressed, decompressed, backend)) {
        return false;
    }

    // Validate
    GVCFBlockDecompressor decompressor(backend);
    return decompressor.Validate(original, decompressed);
}

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Check if data has typical gVCF characteristics.
 */
inline bool IsTypicalGVCF(const std::vector<std::string>& gt_data,
                         const std::vector<std::string>& alt_data) {
    if (gt_data.empty() || alt_data.empty()) {
        return false;
    }

    // Count 0/0 genotypes
    uint32_t hom_ref_count = 0;
    for (const auto& gt : gt_data) {
        if (gt == "0/0" || gt == "0|0") {
            hom_ref_count++;
        }
    }
    float gt_ratio = static_cast<float>(hom_ref_count) / gt_data.size();

    // Count <NON_REF> alleles
    uint32_t non_ref_count = 0;
    for (const auto& alt : alt_data) {
        if (alt == "<NON_REF>") {
            non_ref_count++;
        }
    }
    float alt_ratio = static_cast<float>(non_ref_count) / alt_data.size();

    return gt_ratio >= 0.7f && alt_ratio >= 0.7f;
}

/**
 * Estimate compression ratio for a block.
 */
inline float EstimateCompressionRatio(const GVCFBlock& block,
                                      const GVCFBlockConfig& config = GVCFBlockConfig()) {
    auto analysis = AnalyzeBlock(block, nullptr, config);

    // Estimate based on field characteristics
    float total_ratio = 0.0f;
    int field_count = 0;

    if (!block.position.chrom.empty()) {
        total_ratio += (1.0f - analysis.chrom_analysis.dominant_ratio) * 0.2f;
        field_count++;
    }

    if (!block.sample.gt.empty()) {
        total_ratio += (1.0f - analysis.gt_analysis.dominant_ratio) * 0.1f;
        field_count++;
    }

    if (!block.sequence.alt.empty()) {
        total_ratio += (1.0f - analysis.alt_analysis.dominant_ratio) * 0.15f;
        field_count++;
    }

    // Default estimate if no fields
    if (field_count == 0) {
        return 0.3f;
    }

    return total_ratio / field_count;
}

} // namespace gvcf
