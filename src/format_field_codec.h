#pragma once

#include <vector>
#include <cstdint>
#include <string>
#include <memory>
#include <functional>

namespace gsc {

// ============================================================================
// Field Features - Statistics collected during detection phase
// ============================================================================
struct FieldFeatures {
    // Basic type information
    bool is_array = false;          // Contains comma (array type)
    bool is_numeric = true;         // All values are numeric
    bool is_integer = true;         // All values are integers (no decimal point)

    // Missing rate
    uint32_t total_count = 0;       // Total samples observed
    uint32_t missing_count = 0;     // Count of "." values
    float missing_ratio = 0.0f;     // missing_count / total_count

    // Array features
    bool array_len_fixed = true;    // Array length is constant across samples
    uint32_t array_len = 0;         // If fixed, the array length
    uint32_t array_len_min = UINT32_MAX;
    uint32_t array_len_max = 0;

    // Value range (for integers)
    int64_t min_val = INT64_MAX;
    int64_t max_val = INT64_MIN;
    uint8_t max_bits = 0;           // Bits needed to store max_val

    // Sparsity and cardinality
    uint32_t unique_count = 0;      // Number of unique values
    float top1_freq = 0.0f;         // Frequency of most common value
    uint32_t zero_count = 0;        // Count of all-zero arrays

    // Cross-field correlation (set externally)
    float predictor_hit_ratio = 0.0f;  // Hit ratio for predictor (DP<-AD, GQ<-PL)

    // Pattern detection (for PL-like fields)
    float pattern_hit_ratio = 0.0f;    // Ratio matching known patterns

    // Computed flags for quick codec selection
    bool is_high_missing() const { return missing_ratio >= 0.8f; }
    bool is_predictable() const { return predictor_hit_ratio >= 0.95f; }
    bool is_pattern_match() const { return pattern_hit_ratio >= 0.9f; }
    bool fits_u8() const { return max_val >= 0 && max_val <= 0xFF; }
    bool fits_u16() const { return max_val >= 0 && max_val <= 0xFFFF; }
    bool fits_u32() const { return max_val >= 0 && max_val <= 0xFFFFFFFF; }
};

// ============================================================================
// Codec Type Enumeration
// ============================================================================
enum class CodecType : uint8_t {
    RawString       = 0,    // Fallback: store as-is
    PredictedScalar = 1,    // Predicted scalar (DP/GQ)
    BitTipArray     = 2,    // 2-bit tip array (AD/PL)
    SparseDictString= 3,    // Sparse dictionary string (PGT/PID)
    FixedWidthArray = 4,    // Fixed-width numeric array
    DeltaVarint     = 5,    // Delta + varint encoding
    DictArray       = 6,    // Dictionary-based array
};

// Convert CodecType to string for debugging
inline const char* codecTypeToString(CodecType type) {
    switch (type) {
        case CodecType::RawString:        return "RawString";
        case CodecType::PredictedScalar:  return "PredictedScalar";
        case CodecType::BitTipArray:      return "BitTipArray";
        case CodecType::SparseDictString: return "SparseDictString";
        case CodecType::FixedWidthArray:  return "FixedWidthArray";
        case CodecType::DeltaVarint:      return "DeltaVarint";
        case CodecType::DictArray:        return "DictArray";
        default:                          return "Unknown";
    }
}

// ============================================================================
// Codec Parameters - Serialized with encoded data for decoding
// ============================================================================
struct CodecParams {
    CodecType type = CodecType::RawString;
    uint32_t sample_count = 0;      // Number of samples in this row
    uint32_t array_len = 0;         // For array codecs
    uint8_t element_bytes = 1;      // 1/2/4 bytes per element
    uint32_t dict_size = 0;         // Dictionary size if applicable

    // Serialize to bytes
    void serialize(std::vector<uint8_t>& out) const;

    // Deserialize from bytes, returns bytes consumed
    size_t deserialize(const uint8_t* data, size_t len);
};

// ============================================================================
// Abstract Codec Interface
// ============================================================================
class FormatFieldCodec {
public:
    virtual ~FormatFieldCodec() = default;

    // -------------------------------------------------------------------------
    // Lifecycle: observe -> freeze -> encode* -> serialize
    // -------------------------------------------------------------------------

    // Observation phase: collect statistics for codec parameter selection
    // Called for first N samples before freeze()
    virtual void observe(const char* value, size_t len, uint32_t sample_pos) = 0;

    // Freeze: finalize codec parameters based on observations
    // After this, no more observe() calls; ready for encode()
    virtual void freeze() = 0;

    // Check if codec is frozen
    virtual bool isFrozen() const = 0;

    // -------------------------------------------------------------------------
    // Encoding phase (after freeze)
    // -------------------------------------------------------------------------

    // Encode a single sample's value
    virtual void encode(const char* value, size_t len, uint32_t sample_pos) = 0;

    // Serialize all encoded data to output buffer
    // Format: [CodecParams][Payload]
    virtual void serialize(std::vector<uint8_t>& out) const = 0;

    // -------------------------------------------------------------------------
    // Decoding phase
    // -------------------------------------------------------------------------

    // Deserialize from bytes, returns bytes consumed
    virtual size_t deserialize(const uint8_t* data, size_t len) = 0;

    // Decode a single sample's value
    virtual std::string decode(uint32_t sample_pos) const = 0;

    // -------------------------------------------------------------------------
    // Metadata
    // -------------------------------------------------------------------------

    // Get codec type
    virtual CodecType type() const = 0;

    // Get codec parameters (for serialization header)
    virtual CodecParams getParams() const = 0;

    // Reset state for next row
    virtual void reset() = 0;

    // -------------------------------------------------------------------------
    // Optional: Cross-field predictor support
    // -------------------------------------------------------------------------

    // Set a predictor function for predicted codecs
    // predictor(sample_pos) -> expected_value
    virtual void setPredictor(std::function<int64_t(uint32_t)> /*predictor*/) {}

    // Get computed value for cross-field prediction (e.g., sum(AD) for DP prediction)
    virtual int64_t getComputedValue(uint32_t /*sample_pos*/) const { return 0; }
};

// ============================================================================
// Factory function to create codec by type
// ============================================================================
std::unique_ptr<FormatFieldCodec> createCodec(CodecType type);

// Factory with field type hint (for field-specific optimizations)
std::unique_ptr<FormatFieldCodec> createCodec(CodecType type,
                                               const std::string& field_name,
                                               uint32_t array_len);

// ============================================================================
// Codec selection based on features
// ============================================================================
CodecType selectCodecType(const FieldFeatures& features,
                          bool has_predictor = false,
                          const std::string& field_name = "");

} // namespace gsc
