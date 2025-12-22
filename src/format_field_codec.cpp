#include "format_field_codec.h"
#include "vint_code.h"
#include "codecs/predicted_scalar_codec.h"
#include "codecs/bit_tip_array_codec.h"
#include "codecs/sparse_dict_codec.h"
#include <cstring>
#include <algorithm>

namespace gsc {

// ============================================================================
// CodecParams Implementation
// ============================================================================

void CodecParams::serialize(std::vector<uint8_t>& out) const {
    // Format: [type:1][sample_count:vint][array_len:vint][element_bytes:1][dict_size:vint]
    out.push_back(static_cast<uint8_t>(type));

    // Use vint encoding for variable-length integers
    vint_code::WriteVint(sample_count, out);
    vint_code::WriteVint(array_len, out);
    out.push_back(element_bytes);
    vint_code::WriteVint(dict_size, out);
}

size_t CodecParams::deserialize(const uint8_t* data, size_t len) {
    if (len < 3) return 0;  // Minimum: type + some vints

    size_t pos = 0;

    // Type
    type = static_cast<CodecType>(data[pos++]);

    // Sample count
    sample_count = vint_code::ReadVint(data, len, pos);

    // Array length
    array_len = vint_code::ReadVint(data, len, pos);

    // Element bytes
    if (pos >= len) return 0;
    element_bytes = data[pos++];

    // Dict size
    dict_size = vint_code::ReadVint(data, len, pos);

    return pos;
}

// ============================================================================
// RawStringCodec - Fallback codec that stores values as-is
// ============================================================================

class RawStringCodec : public FormatFieldCodec {
public:
    RawStringCodec() = default;

    void observe(const char* value, size_t len, uint32_t sample_pos) override {
        (void)value;
        (void)len;
        (void)sample_pos;
        ++observed_count_;
    }

    void freeze() override {
        frozen_ = true;
    }

    bool isFrozen() const override {
        return frozen_;
    }

    void encode(const char* value, size_t len, uint32_t sample_pos) override {
        (void)sample_pos;

        // Record offset for this sample
        offsets_.push_back(static_cast<uint32_t>(raw_data_.size()));

        // Append raw bytes
        raw_data_.insert(raw_data_.end(), value, value + len);

        ++sample_count_;
    }

    void serialize(std::vector<uint8_t>& out) const override {
        // Header
        CodecParams params = getParams();
        params.serialize(out);

        // Offsets (vint encoded)
        vint_code::WriteVint(static_cast<uint32_t>(offsets_.size()), out);
        for (uint32_t off : offsets_) {
            vint_code::WriteVint(off, out);
        }

        // Raw data size + data
        vint_code::WriteVint(static_cast<uint32_t>(raw_data_.size()), out);
        out.insert(out.end(), raw_data_.begin(), raw_data_.end());
    }

    size_t deserialize(const uint8_t* data, size_t len) override {
        size_t pos = 0;

        // Skip header (already parsed externally)
        // Read offsets
        uint32_t offset_count = vint_code::ReadVint(data, len, pos);

        offsets_.resize(offset_count);
        for (uint32_t i = 0; i < offset_count; ++i) {
            offsets_[i] = vint_code::ReadVint(data, len, pos);
        }

        // Read raw data
        uint32_t raw_size = vint_code::ReadVint(data, len, pos);

        raw_data_.assign(data + pos, data + pos + raw_size);
        pos += raw_size;

        sample_count_ = offset_count;
        frozen_ = true;

        return pos;
    }

    std::string decode(uint32_t sample_pos) const override {
        if (sample_pos >= offsets_.size()) {
            return ".";
        }

        uint32_t start = offsets_[sample_pos];
        uint32_t end = (sample_pos + 1 < offsets_.size())
                       ? offsets_[sample_pos + 1]
                       : static_cast<uint32_t>(raw_data_.size());

        return std::string(raw_data_.begin() + start, raw_data_.begin() + end);
    }

    CodecType type() const override {
        return CodecType::RawString;
    }

    CodecParams getParams() const override {
        CodecParams params;
        params.type = CodecType::RawString;
        params.sample_count = sample_count_;
        params.array_len = 0;
        params.element_bytes = 1;
        params.dict_size = 0;
        return params;
    }

    void reset() override {
        raw_data_.clear();
        offsets_.clear();
        sample_count_ = 0;
        observed_count_ = 0;
        frozen_ = false;
    }

private:
    std::vector<uint8_t> raw_data_;
    std::vector<uint32_t> offsets_;
    uint32_t sample_count_ = 0;
    uint32_t observed_count_ = 0;
    bool frozen_ = false;
};

// ============================================================================
// Factory function
// ============================================================================

std::unique_ptr<FormatFieldCodec> createCodec(CodecType type) {
    switch (type) {
        case CodecType::RawString:
            return std::make_unique<RawStringCodec>();

        case CodecType::PredictedScalar:
            return std::make_unique<PredictedScalarCodec>();

        case CodecType::BitTipArray:
            return std::make_unique<BitTipArrayCodec>();

        case CodecType::SparseDictString:
            return std::make_unique<SparseDictCodec>();

        // These still fall back to RawString for now
        case CodecType::FixedWidthArray:
        case CodecType::DeltaVarint:
        case CodecType::DictArray:
            // Use BitTipArray for arrays, RawString otherwise
            return std::make_unique<RawStringCodec>();

        default:
            return std::make_unique<RawStringCodec>();
    }
}

// Factory with field type hint for BitTipArrayCodec
std::unique_ptr<FormatFieldCodec> createCodec(CodecType type,
                                               const std::string& field_name,
                                               uint32_t array_len) {
    if (type == CodecType::BitTipArray) {
        BitTipArrayCodec::FieldType ft = BitTipArrayCodec::FieldType::GENERIC;

        // Detect field type from name
        if (field_name == "AD") {
            ft = BitTipArrayCodec::FieldType::AD_LIKE;
        } else if (field_name == "PL") {
            ft = BitTipArrayCodec::FieldType::PL_LIKE;
        }

        return std::make_unique<BitTipArrayCodec>(ft, array_len);
    }

    return createCodec(type);
}

// ============================================================================
// Codec selection based on features
// ============================================================================

CodecType selectCodecType(const FieldFeatures& features,
                          bool has_predictor,
                          const std::string& field_name) {
    (void)field_name;  // May be used for heuristics in future

    // Priority 1: High missing ratio with low cardinality -> SparseDictString
    if (features.is_high_missing() && features.unique_count <= 256) {
        return CodecType::SparseDictString;
    }

    // Priority 2: Predictable scalar (e.g., DP <- sum(AD), GQ <- f(PL))
    if (has_predictor && features.is_predictable() &&
        !features.is_array && features.is_integer) {
        return CodecType::PredictedScalar;
    }

    // Priority 3: Fixed-length array with pattern matching (AD/PL-like)
    if (features.is_array && features.array_len_fixed) {
        // Check for high zero ratio or pattern matching
        float zero_ratio = features.total_count > 0
                           ? static_cast<float>(features.zero_count) / features.total_count
                           : 0.0f;

        if (zero_ratio >= 0.5f || features.is_pattern_match()) {
            return CodecType::BitTipArray;
        }

        // Fixed-width array for moderate cardinality
        if (features.unique_count <= 65536) {
            if (features.fits_u8()) {
                return CodecType::FixedWidthArray;
            } else if (features.fits_u16()) {
                return CodecType::FixedWidthArray;
            }
        }

        // Dictionary-based array for high cardinality
        return CodecType::DictArray;
    }

    // Priority 4: Non-array integer with moderate cardinality
    if (!features.is_array && features.is_integer && features.unique_count <= 256) {
        return CodecType::DeltaVarint;
    }

    // Priority 5: String field with low cardinality
    if (!features.is_numeric && features.unique_count <= 256) {
        return CodecType::SparseDictString;
    }

    // Fallback: RawString
    return CodecType::RawString;
}

} // namespace gsc
