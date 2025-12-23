#include "predicted_scalar_codec.h"
#include "../vint_code.h"
#include <cstdlib>
#include <cstring>
#include <limits>

namespace gsc {

PredictedScalarCodec::PredictedScalarCodec() {
    observed_values_.reserve(256);
    exception_positions_.reserve(64);
    exception_values_.reserve(64);
    missing_positions_.reserve(64);
}

void PredictedScalarCodec::observe(const char* value, size_t len, uint32_t sample_pos) {
    (void)sample_pos;
    if (frozen_) return;

    int64_t val = parseValue(value, len);
    observed_values_.push_back(val);
    ++observed_count_;
}

void PredictedScalarCodec::freeze() {
    if (frozen_) return;
    frozen_ = true;

    // Compute default value from observations if no predictor
    if (!predictor_ && !observed_values_.empty()) {
        // Use median or mode as default
        int64_t sum = 0;
        int64_t count = 0;
        for (int64_t v : observed_values_) {
            if (v != MISSING_MARKER && v != NOT_PRESENT_MARKER) {
                sum += v;
                ++count;
            }
        }
        if (count > 0) {
            default_value_ = sum / count;
        }
    }

    // Process observed values as if encoding
    for (uint32_t i = 0; i < observed_values_.size(); ++i) {
        int64_t actual = observed_values_[i];

        if (actual == MISSING_MARKER) {
            missing_positions_.push_back(i);
        } else if (actual == NOT_PRESENT_MARKER) {
            not_present_positions_.push_back(i);
        } else {
            int64_t predicted = getPredictedValue(i);
            if (actual != predicted) {
                exception_positions_.push_back(i);
                exception_values_.push_back(actual);
            }
        }
    }

    sample_count_ = observed_count_;
    observed_values_.clear();
}

void PredictedScalarCodec::encode(const char* value, size_t len, uint32_t sample_pos) {
    if (!frozen_) {
        observe(value, len, sample_pos);
        return;
    }

    int64_t actual = parseValue(value, len);

    if (actual == MISSING_MARKER) {
        missing_positions_.push_back(sample_pos);
    } else if (actual == NOT_PRESENT_MARKER) {
        not_present_positions_.push_back(sample_pos);
    } else {
        int64_t predicted = getPredictedValue(sample_pos);
        if (actual != predicted) {
            exception_positions_.push_back(sample_pos);
            exception_values_.push_back(actual);
        }
    }

    ++sample_count_;
}

void PredictedScalarCodec::serialize(std::vector<uint8_t>& out) const {
    // Header
    CodecParams params = getParams();
    params.serialize(out);

    // Default value (for decoding without predictor)
    // Encode as zigzag to handle negative values
    int64_t zigzag = (default_value_ << 1) ^ (default_value_ >> 63);
    vint_code::WriteVint(static_cast<uint32_t>(zigzag & 0xFFFFFFFF), out);
    vint_code::WriteVint(static_cast<uint32_t>((zigzag >> 32) & 0xFFFFFFFF), out);

    // Exceptions
    vint_code::WriteVint(static_cast<uint32_t>(exception_positions_.size()), out);
    uint32_t prev_pos = 0;
    for (size_t i = 0; i < exception_positions_.size(); ++i) {
        // Delta encode positions
        uint32_t delta = exception_positions_[i] - prev_pos;
        vint_code::WriteVint(delta, out);
        prev_pos = exception_positions_[i];

        // ZigZag encode values
        int64_t val = exception_values_[i];
        int64_t zz = (val << 1) ^ (val >> 63);
        vint_code::WriteVint(static_cast<uint32_t>(zz & 0xFFFFFFFF), out);
        vint_code::WriteVint(static_cast<uint32_t>((zz >> 32) & 0xFFFFFFFF), out);
    }

    // Missing positions
    vint_code::WriteVint(static_cast<uint32_t>(missing_positions_.size()), out);
    prev_pos = 0;
    for (uint32_t pos : missing_positions_) {
        uint32_t delta = pos - prev_pos;
        vint_code::WriteVint(delta, out);
        prev_pos = pos;
    }

    // Not present positions
    vint_code::WriteVint(static_cast<uint32_t>(not_present_positions_.size()), out);
    prev_pos = 0;
    for (uint32_t pos : not_present_positions_) {
        uint32_t delta = pos - prev_pos;
        vint_code::WriteVint(delta, out);
        prev_pos = pos;
    }
}

size_t PredictedScalarCodec::deserialize(const uint8_t* data, size_t len) {
    size_t pos = 0;

    // Default value
    uint32_t lo = vint_code::ReadVint(data, len, pos);
    uint32_t hi = vint_code::ReadVint(data, len, pos);
    int64_t zigzag = (static_cast<int64_t>(hi) << 32) | lo;
    default_value_ = (zigzag >> 1) ^ -(zigzag & 1);

    // Exceptions
    uint32_t exc_count = vint_code::ReadVint(data, len, pos);
    exception_map_.clear();
    exception_map_.reserve(exc_count);
    uint32_t cur_pos = 0;
    for (uint32_t i = 0; i < exc_count; ++i) {
        uint32_t delta = vint_code::ReadVint(data, len, pos);
        cur_pos += delta;

        lo = vint_code::ReadVint(data, len, pos);
        hi = vint_code::ReadVint(data, len, pos);
        zigzag = (static_cast<int64_t>(hi) << 32) | lo;
        int64_t val = (zigzag >> 1) ^ -(zigzag & 1);

        exception_map_[cur_pos] = val;
    }

    // Missing positions
    uint32_t miss_count = vint_code::ReadVint(data, len, pos);
    missing_set_.clear();
    missing_set_.reserve(miss_count);
    cur_pos = 0;
    for (uint32_t i = 0; i < miss_count; ++i) {
        uint32_t delta = vint_code::ReadVint(data, len, pos);
        cur_pos += delta;
        missing_set_.insert(cur_pos);
    }

    // Not present positions
    uint32_t np_count = vint_code::ReadVint(data, len, pos);
    not_present_set_.clear();
    not_present_set_.reserve(np_count);
    cur_pos = 0;
    for (uint32_t i = 0; i < np_count; ++i) {
        uint32_t delta = vint_code::ReadVint(data, len, pos);
        cur_pos += delta;
        not_present_set_.insert(cur_pos);
    }

    frozen_ = true;
    return pos;
}

std::string PredictedScalarCodec::decode(uint32_t sample_pos) const {
    // Check not present
    if (not_present_set_.count(sample_pos)) {
        return "";  // Field not present for this sample
    }

    // Check missing
    if (missing_set_.count(sample_pos)) {
        return ".";
    }

    // Check exception
    auto it = exception_map_.find(sample_pos);
    if (it != exception_map_.end()) {
        return std::to_string(it->second);
    }

    // Return predicted value
    int64_t predicted = getPredictedValue(sample_pos);
    return std::to_string(predicted);
}

bool PredictedScalarCodec::decodeToInt32(uint32_t sample_pos, int32_t& out) const {
    if (not_present_set_.count(sample_pos) || missing_set_.count(sample_pos)) {
        return false;
    }

    int64_t v = 0;
    auto it = exception_map_.find(sample_pos);
    if (it != exception_map_.end()) {
        v = it->second;
    } else {
        v = getPredictedValue(sample_pos);
    }

    if (v < std::numeric_limits<int32_t>::min() + 2 || v > std::numeric_limits<int32_t>::max()) {
        return false;
    }
    out = static_cast<int32_t>(v);
    return true;
}

CodecParams PredictedScalarCodec::getParams() const {
    CodecParams params;
    params.type = CodecType::PredictedScalar;
    params.sample_count = sample_count_;
    params.array_len = 0;
    params.element_bytes = 4;  // int32 typically
    params.dict_size = static_cast<uint32_t>(exception_positions_.size());
    return params;
}

void PredictedScalarCodec::reset() {
    frozen_ = false;
    sample_count_ = 0;
    exception_positions_.clear();
    exception_values_.clear();
    missing_positions_.clear();
    not_present_positions_.clear();
    exception_map_.clear();
    missing_set_.clear();
    not_present_set_.clear();
    observed_values_.clear();
    observed_count_ = 0;
}

void PredictedScalarCodec::setPredictor(std::function<int64_t(uint32_t)> predictor) {
    predictor_ = predictor;
}

int64_t PredictedScalarCodec::getComputedValue(uint32_t sample_pos) const {
    // Check not present
    if (not_present_set_.count(sample_pos)) {
        return NOT_PRESENT_MARKER;
    }

    // Check missing
    if (missing_set_.count(sample_pos)) {
        return MISSING_MARKER;
    }

    // Check exception
    auto it = exception_map_.find(sample_pos);
    if (it != exception_map_.end()) {
        return it->second;
    }

    // Return predicted value
    return getPredictedValue(sample_pos);
}

int64_t PredictedScalarCodec::parseValue(const char* value, size_t len) const {
    if (len == 0) {
        return NOT_PRESENT_MARKER;
    }

    if (len == 1 && value[0] == '.') {
        return MISSING_MARKER;
    }

    char* endptr;
    int64_t val = std::strtoll(value, &endptr, 10);

    if (endptr != value + len) {
        // Not a pure integer, try as float
        double fval = std::strtod(value, &endptr);
        if (endptr == value + len) {
            return static_cast<int64_t>(fval);
        }
        // Parsing failed, treat as missing
        return MISSING_MARKER;
    }

    return val;
}

int64_t PredictedScalarCodec::getPredictedValue(uint32_t sample_pos) const {
    if (predictor_) {
        int64_t predicted = predictor_(sample_pos);
        if (predicted != MISSING_MARKER && predicted != NOT_PRESENT_MARKER) {
            return predicted;
        }
    }
    return default_value_;
}

} // namespace gsc
