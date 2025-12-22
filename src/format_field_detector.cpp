#include "format_field_detector.h"
#include "vint_code.h"
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cmath>

namespace gsc {

// ============================================================================
// FormatFieldDetector Implementation
// ============================================================================

FormatFieldDetector::FormatFieldDetector(const std::string& field_name,
                                         const DetectorConfig& config)
    : field_name_(field_name)
    , config_(config)
{
    array_lengths_.reserve(config.observe_sample_limit);
    computed_values_.reserve(config.observe_sample_limit);
}

void FormatFieldDetector::feed(const char* value, size_t len) {
    if (finalized_ || sample_count_ >= config_.observe_sample_limit) {
        return;
    }

    ++sample_count_;

    // Check for missing value
    if (len == 1 && value[0] == '.') {
        ++missing_count_;
        computed_values_.push_back(INT64_MIN);  // Mark as missing
        return;
    }

    // Parse and analyze the value
    parseValue(value, len);
}

void FormatFieldDetector::setExpectedArrayLen(uint32_t len) {
    config_.expected_array_len = len;
}

bool FormatFieldDetector::isComplete() const {
    return sample_count_ >= config_.observe_sample_limit || finalized_;
}

void FormatFieldDetector::setPredictor(std::function<int64_t(uint32_t)> predictor) {
    predictor_ = predictor;
}

int64_t FormatFieldDetector::getComputedValue(uint32_t sample_pos) const {
    if (sample_pos < computed_values_.size()) {
        return computed_values_[sample_pos];
    }
    return INT64_MIN;
}

FieldFeatures FormatFieldDetector::finalize() {
    if (finalized_) {
        return features_;
    }

    finalized_ = true;

    // Basic counts
    features_.total_count = sample_count_;
    features_.missing_count = missing_count_;
    features_.missing_ratio = sample_count_ > 0
        ? static_cast<float>(missing_count_) / sample_count_
        : 0.0f;

    // Type information
    features_.is_array = has_comma_;
    features_.is_numeric = all_numeric_;
    features_.is_integer = all_integer_;

    // Value range
    features_.min_val = min_val_;
    features_.max_val = max_val_;
    features_.max_bits = computeBitsNeeded(max_val_);

    // Array length analysis
    if (!array_lengths_.empty()) {
        features_.array_len_min = *std::min_element(array_lengths_.begin(), array_lengths_.end());
        features_.array_len_max = *std::max_element(array_lengths_.begin(), array_lengths_.end());
        features_.array_len_fixed = (features_.array_len_min == features_.array_len_max);
        features_.array_len = features_.array_len_max;
    }

    // Cardinality
    features_.unique_count = static_cast<uint32_t>(value_freq_.size()) + overflow_count_;

    // Top frequency
    if (!value_freq_.empty()) {
        uint32_t max_freq = 0;
        for (const auto& kv : value_freq_) {
            max_freq = std::max(max_freq, kv.second);
        }
        features_.top1_freq = static_cast<float>(max_freq) / sample_count_;
    }

    // Zero array ratio
    features_.zero_count = zero_array_count_;

    // Predictor hit ratio
    if (predictor_ && (predictor_hit_count_ + predictor_miss_count_) > 0) {
        features_.predictor_hit_ratio = static_cast<float>(predictor_hit_count_)
            / (predictor_hit_count_ + predictor_miss_count_);
    }

    // Pattern hit ratio
    uint32_t non_missing = sample_count_ - missing_count_;
    if (non_missing > 0) {
        features_.pattern_hit_ratio = static_cast<float>(pattern_hit_count_) / non_missing;
    }

    return features_;
}

std::unique_ptr<FormatFieldCodec> FormatFieldDetector::createCodec() {
    if (!finalized_) {
        finalize();
    }

    CodecType type = recommendedCodecType();
    auto codec = gsc::createCodec(type);

    // Set predictor if applicable
    if (predictor_ && type == CodecType::PredictedScalar) {
        codec->setPredictor(predictor_);
    }

    return codec;
}

CodecType FormatFieldDetector::recommendedCodecType() const {
    return selectCodecType(features_, predictor_ != nullptr, field_name_);
}

void FormatFieldDetector::parseValue(const char* value, size_t len) {
    // Check if this is an array (contains comma)
    bool is_array = (std::memchr(value, ',', len) != nullptr);

    if (is_array) {
        has_comma_ = true;
        parseArrayValue(value, len);
    } else {
        parseScalarValue(value, len);
    }
}

void FormatFieldDetector::parseArrayValue(const char* value, size_t len) {
    std::vector<int64_t> elements;
    elements.reserve(16);

    const char* start = value;
    const char* end = value + len;
    const char* p = start;
    bool all_int = true;
    bool all_num = true;
    bool is_all_zero = true;

    while (p <= end) {
        if (p == end || *p == ',') {
            size_t elem_len = p - start;
            if (elem_len > 0) {
                // Parse element
                if (elem_len == 1 && *start == '.') {
                    // Missing element
                    elements.push_back(INT64_MIN);
                } else {
                    char* endptr;
                    int64_t val = std::strtoll(start, &endptr, 10);

                    if (endptr == start + elem_len) {
                        // Valid integer
                        elements.push_back(val);
                        if (val != 0) is_all_zero = false;
                        min_val_ = std::min(min_val_, val);
                        max_val_ = std::max(max_val_, val);
                    } else {
                        // Try as float
                        double fval = std::strtod(start, &endptr);
                        if (endptr == start + elem_len) {
                            all_int = false;
                            elements.push_back(static_cast<int64_t>(fval));
                            if (fval != 0.0) is_all_zero = false;
                        } else {
                            all_num = false;
                            all_int = false;
                        }
                    }
                }
            }
            start = p + 1;
        }
        ++p;
    }

    if (!all_num) all_numeric_ = false;
    if (!all_int) all_integer_ = false;

    // Track array length
    array_lengths_.push_back(static_cast<uint32_t>(elements.size()));

    // Track zero arrays
    if (is_all_zero && !elements.empty()) {
        ++zero_array_count_;
    }

    // Compute sum for cross-field prediction (e.g., sum(AD) for DP)
    int64_t sum = computeArraySum(elements);
    computed_values_.push_back(sum);

    // Check predictor if set
    if (predictor_) {
        int64_t expected = predictor_(sample_count_ - 1);
        if (expected != INT64_MIN && sum == expected) {
            ++predictor_hit_count_;
        } else if (expected != INT64_MIN) {
            ++predictor_miss_count_;
        }
    }

    // Track unique values (as string)
    if (value_freq_.size() < config_.unique_value_limit) {
        std::string key(value, len);
        ++value_freq_[key];
    } else {
        ++overflow_count_;
    }

    // Check AB pattern for PL-like fields
    if (all_int && elements.size() >= 3) {
        if (checkAbPattern(elements)) {
            ++pattern_hit_count_;
        }
    }
}

void FormatFieldDetector::parseScalarValue(const char* value, size_t len) {
    // Parse scalar value
    char* endptr;
    int64_t val = std::strtoll(value, &endptr, 10);

    if (endptr == value + len) {
        // Valid integer
        min_val_ = std::min(min_val_, val);
        max_val_ = std::max(max_val_, val);
        computed_values_.push_back(val);
    } else {
        // Try as float
        double fval = std::strtod(value, &endptr);
        if (endptr == value + len) {
            all_integer_ = false;
            computed_values_.push_back(static_cast<int64_t>(fval));
        } else {
            all_numeric_ = false;
            all_integer_ = false;
            computed_values_.push_back(INT64_MIN);
        }
    }

    // Check predictor
    if (predictor_ && !computed_values_.empty()) {
        int64_t expected = predictor_(sample_count_ - 1);
        int64_t actual = computed_values_.back();
        if (expected != INT64_MIN && actual != INT64_MIN && expected == actual) {
            ++predictor_hit_count_;
        } else if (expected != INT64_MIN && actual != INT64_MIN) {
            ++predictor_miss_count_;
        }
    }

    // Track unique values
    if (value_freq_.size() < config_.unique_value_limit) {
        std::string key(value, len);
        ++value_freq_[key];
    } else {
        ++overflow_count_;
    }
}

bool FormatFieldDetector::checkAbPattern(const std::vector<int64_t>& arr) {
    // PL pattern check from ref_code/fmt_comp
    // Pattern: [0, a, b, ...] where subsequent elements follow a formula
    // Simplified check: arr[0] == 0 and arr[1] > 0 and arr[2] > 0

    if (arr.size() < 3) return false;
    if (arr[0] != 0) return false;
    if (arr[1] <= 0 || arr[2] <= 0) return false;

    // Check if b == 15 * a (stronger pattern)
    int64_t a = arr[1];
    int64_t b = arr[2];
    if (b == 15 * a) return true;

    // General AB pattern: check if remaining elements follow expected formula
    // For now, accept if first 3 elements match basic pattern
    return true;
}

int64_t FormatFieldDetector::computeArraySum(const std::vector<int64_t>& arr) {
    int64_t sum = 0;
    for (int64_t v : arr) {
        if (v != INT64_MIN) {
            sum += v;
        }
    }
    return sum;
}

uint8_t FormatFieldDetector::computeBitsNeeded(int64_t max_val) {
    if (max_val <= 0) return 0;
    uint8_t bits = 0;
    uint64_t v = static_cast<uint64_t>(max_val);
    while (v > 0) {
        ++bits;
        v >>= 1;
    }
    return bits;
}

// ============================================================================
// FormatFieldManager Implementation
// ============================================================================

FormatFieldManager::FormatFieldManager() = default;

void FormatFieldManager::initRow(const std::vector<std::string>& format_keys,
                                  uint32_t allele_count) {
    current_format_keys_ = format_keys;
    allele_count_ = allele_count;
    samples_processed_ = 0;
    frozen_ = false;

    // Clear previous state
    detectors_.clear();
    codecs_.clear();
    sample_buffer_.clear();

    // Create detectors for each field (except GT if skip_gt_ is true)
    for (const auto& key : format_keys) {
        if (skip_gt_ && key == "GT") {
            continue;
        }

        DetectorConfig field_config = config_;
        field_config.expected_array_len = getExpectedArrayLen(key);

        detectors_[key] = std::make_unique<FormatFieldDetector>(key, field_config);
    }

    // Setup cross-field predictors
    setupCrossFieldPredictors();
}

void FormatFieldManager::processSample(const char* sample_str, size_t len,
                                        uint32_t sample_pos) {
    (void)sample_pos;

    // Parse sample fields
    std::unordered_map<std::string, std::string> fields;
    parseSampleFields(sample_str, len, fields);

    if (!frozen_) {
        // Observation phase: feed to detectors
        for (auto& kv : fields) {
            auto it = detectors_.find(kv.first);
            if (it != detectors_.end()) {
                it->second->feed(kv.second.c_str(), kv.second.size());
            }
        }

        // Buffer for later encoding
        BufferedSample bs;
        bs.fields = std::move(fields);
        sample_buffer_.push_back(std::move(bs));

        ++samples_processed_;

        // Check if we should freeze
        bool all_complete = true;
        for (const auto& det : detectors_) {
            if (!det.second->isComplete()) {
                all_complete = false;
                break;
            }
        }

        if (all_complete && samples_processed_ >= config_.observe_sample_limit) {
            freezeAndFlush();
        }
    } else {
        // Encoding phase: directly encode
        for (auto& kv : fields) {
            auto it = codecs_.find(kv.first);
            if (it != codecs_.end()) {
                it->second->encode(kv.second.c_str(), kv.second.size(), samples_processed_);
            }
        }
        ++samples_processed_;
    }
}

void FormatFieldManager::finalizeRow(std::vector<uint8_t>& out) {
    // Force freeze if not already done
    if (!frozen_) {
        freezeAndFlush();
    }

    // Serialize format: [field_count][field_name+codec_data]*
    uint32_t field_count = static_cast<uint32_t>(codecs_.size());
    vint_code::WriteVint(field_count, out);

    for (const auto& key : current_format_keys_) {
        if (skip_gt_ && key == "GT") continue;

        auto it = codecs_.find(key);
        if (it == codecs_.end()) continue;

        // Field name
        vint_code::WriteVint(static_cast<uint32_t>(key.size()), out);
        out.insert(out.end(), key.begin(), key.end());

        // Codec data
        std::vector<uint8_t> codec_data;
        it->second->serialize(codec_data);
        vint_code::WriteVint(static_cast<uint32_t>(codec_data.size()), out);
        out.insert(out.end(), codec_data.begin(), codec_data.end());
    }
}

std::string FormatFieldManager::decodeSample(const uint8_t* data, size_t len,
                                              uint32_t sample_pos,
                                              const std::vector<std::string>& format_keys) {
    // Decode format: [field_count][field_name+codec_data]*
    size_t pos = 0;
    uint32_t field_count = vint_code::ReadVint(data, len, pos);

    std::unordered_map<std::string, std::string> decoded_fields;

    for (uint32_t i = 0; i < field_count && pos < len; ++i) {
        // Field name
        uint32_t name_len = vint_code::ReadVint(data, len, pos);
        std::string field_name(reinterpret_cast<const char*>(data + pos), name_len);
        pos += name_len;

        // Codec data size
        uint32_t codec_size = vint_code::ReadVint(data, len, pos);

        // Parse codec params to determine type
        CodecParams params;
        size_t params_size = params.deserialize(data + pos, codec_size);
        (void)params_size;

        // Create codec and deserialize
        auto codec = gsc::createCodec(params.type);
        codec->deserialize(data + pos, codec_size);
        pos += codec_size;

        // Decode value
        decoded_fields[field_name] = codec->decode(sample_pos);
    }

    // Build output string in format key order
    std::string result;
    for (size_t i = 0; i < format_keys.size(); ++i) {
        if (i > 0) result += ':';

        const std::string& key = format_keys[i];
        auto it = decoded_fields.find(key);
        if (it != decoded_fields.end()) {
            result += it->second;
        } else {
            result += '.';  // Missing field
        }
    }

    return result;
}

void FormatFieldManager::parseSampleFields(const char* sample_str, size_t len,
                                            std::unordered_map<std::string, std::string>& fields) {
    // Parse "val1:val2:val3:..." into {key1: val1, key2: val2, ...}
    size_t key_idx = 0;
    const char* start = sample_str;
    const char* end = sample_str + len;
    const char* p = start;

    while (p <= end && key_idx < current_format_keys_.size()) {
        if (p == end || *p == ':') {
            const std::string& key = current_format_keys_[key_idx];
            fields[key] = std::string(start, p - start);

            ++key_idx;
            start = p + 1;
        }
        ++p;
    }
}

void FormatFieldManager::freezeAndFlush() {
    if (frozen_) return;

    frozen_ = true;

    // Finalize detectors and create codecs
    for (auto& kv : detectors_) {
        kv.second->finalize();
        codecs_[kv.first] = kv.second->createCodec();
        codecs_[kv.first]->freeze();
    }

    // Flush buffered samples
    uint32_t pos = 0;
    for (const auto& bs : sample_buffer_) {
        for (const auto& field : bs.fields) {
            auto it = codecs_.find(field.first);
            if (it != codecs_.end()) {
                it->second->encode(field.second.c_str(), field.second.size(), pos);
            }
        }
        ++pos;
    }

    sample_buffer_.clear();
}

void FormatFieldManager::setupCrossFieldPredictors() {
    // Setup DP <- sum(AD) predictor
    auto ad_it = detectors_.find("AD");
    auto dp_it = detectors_.find("DP");

    if (ad_it != detectors_.end() && dp_it != detectors_.end()) {
        FormatFieldDetector* ad_detector = ad_it->second.get();
        dp_it->second->setPredictor([ad_detector](uint32_t pos) -> int64_t {
            return ad_detector->getComputedValue(pos);
        });
    }

    // Setup GQ <- f(PL) predictor could be added here
    // For simplicity, we'll implement this in Phase 2
}

uint32_t FormatFieldManager::getExpectedArrayLen(const std::string& field_name) const {
    // AD: allele_count elements
    if (field_name == "AD") {
        return allele_count_;
    }

    // PL: n*(n+1)/2 elements where n = allele_count
    if (field_name == "PL") {
        return allele_count_ * (allele_count_ + 1) / 2;
    }

    // Default: unknown
    return 0;
}

} // namespace gsc
