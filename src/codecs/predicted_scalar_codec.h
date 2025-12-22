#pragma once

#include "../format_field_codec.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <functional>

namespace gsc {

// ============================================================================
// PredictedScalarCodec - For fields like DP (predicted by sum(AD)) and GQ (predicted by f(PL))
//
// Encoding strategy:
// - Use a predictor function to compute expected value
// - Only store exceptions where actual != predicted
// - Missing values (.) are stored as special marker
//
// Format:
//   [CodecParams]
//   [exception_count: vint]
//   [exceptions: (pos, value) pairs as vint]
//   [missing_count: vint]
//   [missing_positions: vint array]
// ============================================================================

class PredictedScalarCodec : public FormatFieldCodec {
public:
    // Missing value marker
    static constexpr int64_t MISSING_MARKER = INT64_MIN;
    static constexpr int64_t NOT_PRESENT_MARKER = INT64_MIN + 1;

    PredictedScalarCodec();

    // -------------------------------------------------------------------------
    // FormatFieldCodec interface
    // -------------------------------------------------------------------------

    void observe(const char* value, size_t len, uint32_t sample_pos) override;
    void freeze() override;
    bool isFrozen() const override { return frozen_; }

    void encode(const char* value, size_t len, uint32_t sample_pos) override;
    void serialize(std::vector<uint8_t>& out) const override;

    size_t deserialize(const uint8_t* data, size_t len) override;
    std::string decode(uint32_t sample_pos) const override;

    CodecType type() const override { return CodecType::PredictedScalar; }
    CodecParams getParams() const override;
    void reset() override;

    // Set predictor function: predictor(sample_pos) -> expected_value
    void setPredictor(std::function<int64_t(uint32_t)> predictor) override;

    // Get the actual value for this sample (for use as predictor by other fields)
    int64_t getComputedValue(uint32_t sample_pos) const override;

private:
    bool frozen_ = false;
    uint32_t sample_count_ = 0;

    // Predictor function
    std::function<int64_t(uint32_t)> predictor_;

    // Default predicted value (used when no predictor, or predictor returns invalid)
    int64_t default_value_ = 0;

    // Exception storage: positions where actual != predicted
    std::vector<uint32_t> exception_positions_;
    std::vector<int64_t> exception_values_;

    // Missing value positions (value was ".")
    std::vector<uint32_t> missing_positions_;

    // Not present positions (field not in FORMAT for this sample)
    std::vector<uint32_t> not_present_positions_;

    // For decoding: exception lookup
    std::unordered_map<uint32_t, int64_t> exception_map_;
    std::unordered_set<uint32_t> missing_set_;
    std::unordered_set<uint32_t> not_present_set_;

    // Observation phase data
    std::vector<int64_t> observed_values_;
    uint32_t observed_count_ = 0;

    // Helper functions
    int64_t parseValue(const char* value, size_t len) const;
    int64_t getPredictedValue(uint32_t sample_pos) const;
};

} // namespace gsc
