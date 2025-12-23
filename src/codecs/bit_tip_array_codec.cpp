#include "bit_tip_array_codec.h"
#include "../vint_code.h"
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <algorithm>

namespace gsc {

BitTipArrayCodec::BitTipArrayCodec(FieldType field_type, uint32_t expected_len)
    : field_type_(field_type)
    , expected_len_(expected_len)
{
    observed_arrays_.reserve(256);
    sum_cache_.reserve(256);
}

void BitTipArrayCodec::observe(const char* value, size_t len, uint32_t sample_pos) {
    (void)sample_pos;
    if (frozen_) return;

    ObservedArray obs;

    // Check for missing
    if (len == 1 && value[0] == '.') {
        missing_positions_.push_back(static_cast<uint32_t>(observed_arrays_.size()));
        obs.sum = INT64_MIN;
        observed_arrays_.push_back(std::move(obs));
        return;
    }

    if (!parseArray(value, len, obs.values)) {
        // Parse failed, treat as missing
        missing_positions_.push_back(static_cast<uint32_t>(observed_arrays_.size()));
        obs.sum = INT64_MIN;
    } else {
        // Compute sum
        obs.sum = 0;
        for (int64_t v : obs.values) {
            if (v != INT64_MIN) {  // Skip missing elements
                obs.sum += v;
            }
        }
    }

    observed_arrays_.push_back(std::move(obs));
}

void BitTipArrayCodec::freeze() {
    if (frozen_) return;
    frozen_ = true;

    // Build missing set
    for (uint32_t pos : missing_positions_) {
        missing_set_.insert(pos);
    }

    // Determine array length from observations
    if (expected_len_ == 0) {
        for (const auto& obs : observed_arrays_) {
            if (!obs.values.empty()) {
                expected_len_ = static_cast<uint32_t>(obs.values.size());
                break;
            }
        }
    }

    // Process observed arrays
    sample_count_ = static_cast<uint32_t>(observed_arrays_.size());

    // Initialize tips storage
    size_t tips_bytes = (sample_count_ + 3) / 4;  // 2 bits per sample
    tips_.resize(tips_bytes, 0);

    // Classify each sample and encode
    uint32_t special1_idx = 0, special2_idx = 0;
    for (uint32_t i = 0; i < sample_count_; ++i) {
        if (missing_set_.count(i)) {
            // Missing: use TIP_GENERAL with special marker
            setTip(i, TIP_GENERAL);
            dict_ids_.push_back(UINT32_MAX);  // Special marker for missing
            sum_cache_.push_back(INT64_MIN);
            continue;
        }

        const auto& arr = observed_arrays_[i].values;
        int64_t sum = observed_arrays_[i].sum;
        sum_cache_.push_back(sum);

        TipType tip;
        int64_t a = 0, b = 0;

        switch (field_type_) {
            case FieldType::AD_LIKE:
                tip = classifyAD(arr, sum);
                break;
            case FieldType::PL_LIKE:
                tip = classifyPL(arr, a, b);
                break;
            default:
                tip = classifyGeneric(arr);
                break;
        }

        setTip(i, tip);

        switch (tip) {
            case TIP_ALL_ZERO:
                // No payload needed
                break;

            case TIP_SPECIAL_1:
                if (field_type_ == FieldType::PL_LIKE) {
                    // Store a and b
                    special1_payload_.push_back(static_cast<uint32_t>(a));
                    special2_payload_.push_back(static_cast<uint32_t>(b));
                }
                // For AD_LIKE, no payload (sum==2 is implicit)
                break;

            case TIP_SPECIAL_2:
                if (field_type_ == FieldType::AD_LIKE) {
                    // Store sum
                    special1_payload_.push_back(static_cast<uint32_t>(sum));
                } else if (field_type_ == FieldType::PL_LIKE) {
                    // Store a only (b = 15*a)
                    special1_payload_.push_back(static_cast<uint32_t>(a));
                }
                break;

            case TIP_GENERAL:
                // Store in dictionary
                dict_ids_.push_back(getOrCreateDictId(arr));
                break;
        }
    }

    observed_arrays_.clear();
}

void BitTipArrayCodec::encode(const char* value, size_t len, uint32_t sample_pos) {
    if (!frozen_) {
        observe(value, len, sample_pos);
        return;
    }

    // Extend tips if needed
    size_t tips_bytes = (sample_pos + 4) / 4;
    if (tips_.size() < tips_bytes) {
        tips_.resize(tips_bytes, 0);
    }

    // Check for missing
    if (len == 1 && value[0] == '.') {
        missing_positions_.push_back(sample_pos);
        missing_set_.insert(sample_pos);
        setTip(sample_pos, TIP_GENERAL);
        dict_ids_.push_back(UINT32_MAX);
        sum_cache_.push_back(INT64_MIN);
        ++sample_count_;
        return;
    }

    std::vector<int64_t> arr;
    if (!parseArray(value, len, arr)) {
        // Parse failed
        missing_positions_.push_back(sample_pos);
        missing_set_.insert(sample_pos);
        setTip(sample_pos, TIP_GENERAL);
        dict_ids_.push_back(UINT32_MAX);
        sum_cache_.push_back(INT64_MIN);
        ++sample_count_;
        return;
    }

    // Compute sum
    int64_t sum = 0;
    for (int64_t v : arr) {
        if (v != INT64_MIN) sum += v;
    }
    sum_cache_.push_back(sum);

    TipType tip;
    int64_t a = 0, b = 0;

    switch (field_type_) {
        case FieldType::AD_LIKE:
            tip = classifyAD(arr, sum);
            break;
        case FieldType::PL_LIKE:
            tip = classifyPL(arr, a, b);
            break;
        default:
            tip = classifyGeneric(arr);
            break;
    }

    setTip(sample_pos, tip);

    switch (tip) {
        case TIP_ALL_ZERO:
            break;
        case TIP_SPECIAL_1:
            if (field_type_ == FieldType::PL_LIKE) {
                special1_payload_.push_back(static_cast<uint32_t>(a));
                special2_payload_.push_back(static_cast<uint32_t>(b));
            }
            break;
        case TIP_SPECIAL_2:
            if (field_type_ == FieldType::AD_LIKE) {
                special1_payload_.push_back(static_cast<uint32_t>(sum));
            } else if (field_type_ == FieldType::PL_LIKE) {
                special1_payload_.push_back(static_cast<uint32_t>(a));
            }
            break;
        case TIP_GENERAL:
            dict_ids_.push_back(getOrCreateDictId(arr));
            break;
    }

    ++sample_count_;
}

void BitTipArrayCodec::serialize(std::vector<uint8_t>& out) const {
    // Header
    CodecParams params = getParams();
    params.serialize(out);

    // Field type and expected length
    out.push_back(static_cast<uint8_t>(field_type_));
    vint_code::WriteVint(expected_len_, out);

    // Tips
    vint_code::WriteVint(sample_count_, out);
    vint_code::WriteVint(static_cast<uint32_t>(tips_.size()), out);
    out.insert(out.end(), tips_.begin(), tips_.end());

    // Special1 payload
    vint_code::WriteVint(static_cast<uint32_t>(special1_payload_.size()), out);
    for (uint32_t v : special1_payload_) {
        vint_code::WriteVint(v, out);
    }

    // Special2 payload
    vint_code::WriteVint(static_cast<uint32_t>(special2_payload_.size()), out);
    for (uint32_t v : special2_payload_) {
        vint_code::WriteVint(v, out);
    }

    // Dictionary IDs
    vint_code::WriteVint(static_cast<uint32_t>(dict_ids_.size()), out);
    for (uint32_t id : dict_ids_) {
        vint_code::WriteVint(id, out);
    }

    // Dictionary content
    serializeDictionary(out);

    // Missing positions (for reconstruction)
    vint_code::WriteVint(static_cast<uint32_t>(missing_positions_.size()), out);
    uint32_t prev = 0;
    for (uint32_t pos : missing_positions_) {
        vint_code::WriteVint(pos - prev, out);
        prev = pos;
    }
}

size_t BitTipArrayCodec::deserialize(const uint8_t* data, size_t len) {
    size_t pos = 0;

    // Field type and expected length
    field_type_ = static_cast<FieldType>(data[pos++]);
    expected_len_ = vint_code::ReadVint(data, len, pos);

    // Tips
    sample_count_ = vint_code::ReadVint(data, len, pos);
    uint32_t tips_size = vint_code::ReadVint(data, len, pos);
    tips_.assign(data + pos, data + pos + tips_size);
    pos += tips_size;

    // Special1 payload
    uint32_t s1_size = vint_code::ReadVint(data, len, pos);
    special1_payload_.resize(s1_size);
    for (uint32_t i = 0; i < s1_size; ++i) {
        special1_payload_[i] = vint_code::ReadVint(data, len, pos);
    }

    // Special2 payload
    uint32_t s2_size = vint_code::ReadVint(data, len, pos);
    special2_payload_.resize(s2_size);
    for (uint32_t i = 0; i < s2_size; ++i) {
        special2_payload_[i] = vint_code::ReadVint(data, len, pos);
    }

    // Dictionary IDs
    uint32_t dict_id_count = vint_code::ReadVint(data, len, pos);
    dict_ids_.resize(dict_id_count);
    for (uint32_t i = 0; i < dict_id_count; ++i) {
        dict_ids_[i] = vint_code::ReadVint(data, len, pos);
    }

    // Dictionary content
    pos += deserializeDictionary(data + pos, len - pos);

    // Missing positions
    uint32_t miss_count = vint_code::ReadVint(data, len, pos);
    missing_set_.clear();
    uint32_t cur = 0;
    for (uint32_t i = 0; i < miss_count; ++i) {
        cur += vint_code::ReadVint(data, len, pos);
        missing_set_.insert(cur);
    }

    // Rebuild sum cache for cross-field predictors (e.g., DP <- sum(AD)).
    // This must be available after deserialize() to enable correct decoding of PredictedScalarCodec.
    sum_cache_.clear();
    sum_cache_.resize(sample_count_, INT64_MIN);

    uint32_t special1_idx = 0, special2_idx = 0, dict_idx = 0;
    std::vector<int64_t> arr;
    arr.reserve(expected_len_ > 0 ? expected_len_ : 8);

    for (uint32_t i = 0; i < sample_count_; ++i) {
        if (missing_set_.count(i)) {
            sum_cache_[i] = INT64_MIN;
            ++dict_idx;  // missing samples occupy a dict_ids_ slot
            continue;
        }

        TipType tip = getTip(i);
        reconstructArray(i, tip, special1_idx, special2_idx, dict_idx, arr);

        int64_t sum = 0;
        for (int64_t v : arr) {
            if (v != INT64_MIN) sum += v;
        }
        sum_cache_[i] = sum;

        switch (tip) {
            case TIP_SPECIAL_1:
                if (field_type_ == FieldType::PL_LIKE) {
                    ++special1_idx;
                    ++special2_idx;
                }
                break;
            case TIP_SPECIAL_2:
                ++special1_idx;
                break;
            case TIP_GENERAL:
                ++dict_idx;
                break;
            default:
                break;
        }
    }

    frozen_ = true;
    return pos;
}

std::string BitTipArrayCodec::decode(uint32_t sample_pos) const {
    if (missing_set_.count(sample_pos)) {
        return ".";
    }

    TipType tip = getTip(sample_pos);
    std::vector<int64_t> arr;

    // Count payloads up to this position
    uint32_t special1_idx = 0, special2_idx = 0, dict_idx = 0;

    for (uint32_t i = 0; i < sample_pos && i < sample_count_; ++i) {
        if (missing_set_.count(i)) {
            ++dict_idx;
            continue;
        }
        TipType t = getTip(i);
        switch (t) {
            case TIP_SPECIAL_1:
                if (field_type_ == FieldType::PL_LIKE) {
                    ++special1_idx;
                    ++special2_idx;
                }
                break;
            case TIP_SPECIAL_2:
                ++special1_idx;
                break;
            case TIP_GENERAL:
                ++dict_idx;
                break;
            default:
                break;
        }
    }

    reconstructArray(sample_pos, tip, special1_idx, special2_idx, dict_idx, arr);
    return arrayToString(arr);
}

CodecParams BitTipArrayCodec::getParams() const {
    CodecParams params;
    params.type = CodecType::BitTipArray;
    params.sample_count = sample_count_;
    params.array_len = expected_len_;
    params.element_bytes = 2;  // Default to uint16
    params.dict_size = static_cast<uint32_t>(dictionary_.size());
    return params;
}

void BitTipArrayCodec::reset() {
    frozen_ = false;
    sample_count_ = 0;
    tips_.clear();
    special1_payload_.clear();
    special2_payload_.clear();
    dictionary_.clear();
    dict_index_.clear();
    dict_ids_.clear();
    sum_cache_.clear();
    observed_arrays_.clear();
    missing_positions_.clear();
    missing_set_.clear();
}

int64_t BitTipArrayCodec::getComputedValue(uint32_t sample_pos) const {
    if (sample_pos < sum_cache_.size()) {
        return sum_cache_[sample_pos];
    }
    return INT64_MIN;
}

bool BitTipArrayCodec::parseArray(const char* value, size_t len,
                                   std::vector<int64_t>& out) const {
    out.clear();
    out.reserve(expected_len_ > 0 ? expected_len_ : 8);

    const char* start = value;
    const char* end = value + len;
    const char* p = start;

    while (p <= end) {
        if (p == end || *p == ',') {
            size_t elem_len = p - start;
            if (elem_len > 0) {
                if (elem_len == 1 && *start == '.') {
                    out.push_back(INT64_MIN);  // Missing element
                } else {
                    char* endptr;
                    int64_t val = std::strtoll(start, &endptr, 10);
                    if (endptr == start + elem_len) {
                        out.push_back(val);
                    } else {
                        // Try as float
                        double fval = std::strtod(start, &endptr);
                        if (endptr == start + elem_len) {
                            out.push_back(static_cast<int64_t>(fval));
                        } else {
                            return false;  // Parse error
                        }
                    }
                }
            }
            start = p + 1;
        }
        ++p;
    }

    return !out.empty();
}

BitTipArrayCodec::TipType BitTipArrayCodec::classifyAD(
    const std::vector<int64_t>& arr, int64_t sum) const {

    // Check all zeros
    bool all_zero = true;
    for (int64_t v : arr) {
        if (v != 0 && v != INT64_MIN) {
            all_zero = false;
            break;
        }
    }
    if (all_zero) return TIP_ALL_ZERO;

    // Check sum == arr[0]
    if (!arr.empty() && arr[0] != INT64_MIN && sum == arr[0]) {
        if (sum == 2) {
            return TIP_SPECIAL_1;  // sum==arr[0]==2, no payload
        } else {
            return TIP_SPECIAL_2;  // sum==arr[0]!=2, store sum
        }
    }

    return TIP_GENERAL;
}

BitTipArrayCodec::TipType BitTipArrayCodec::classifyPL(
    const std::vector<int64_t>& arr, int64_t& a, int64_t& b) const {

    // Check all zeros
    bool all_zero = true;
    for (int64_t v : arr) {
        if (v != 0 && v != INT64_MIN) {
            all_zero = false;
            break;
        }
    }
    if (all_zero) return TIP_ALL_ZERO;

    // Check AB pattern
    if (checkAbPattern(arr, a, b)) {
        if (b == 15 * a) {
            return TIP_SPECIAL_2;  // A*15 pattern
        }
        return TIP_SPECIAL_1;  // AB pattern
    }

    return TIP_GENERAL;
}

BitTipArrayCodec::TipType BitTipArrayCodec::classifyGeneric(
    const std::vector<int64_t>& arr) const {

    // Check all zeros
    bool all_zero = true;
    for (int64_t v : arr) {
        if (v != 0 && v != INT64_MIN) {
            all_zero = false;
            break;
        }
    }
    if (all_zero) return TIP_ALL_ZERO;

    return TIP_GENERAL;
}

bool BitTipArrayCodec::checkAbPattern(const std::vector<int64_t>& arr,
                                       int64_t& a, int64_t& b) const {
    // PL pattern: [0, a, b, ...] where:
    // - arr[0] == 0
    // - arr[1] = a (second minimum after 0)
    // - arr[2] = b
    // - Remaining elements follow specific formula

    // For correctness, only treat as special PL pattern when biallelic (len==3).
    if (expected_len_ != 3) return false;
    if (arr.size() != 3) return false;
    if (arr[0] != 0) return false;
    if (arr[1] == INT64_MIN || arr[2] == INT64_MIN) return false;
    if (arr[1] < 0 || arr[2] < 0) return false;

    a = arr[1];
    b = arr[2];

    // For now, accept basic pattern without full verification
    // Full verification would check that remaining elements follow formula
    return true;
}

void BitTipArrayCodec::setTip(uint32_t sample_pos, TipType tip) {
    size_t byte_idx = sample_pos / 4;
    size_t bit_offset = (sample_pos % 4) * 2;

    if (byte_idx >= tips_.size()) {
        tips_.resize(byte_idx + 1, 0);
    }

    // Clear existing bits and set new value
    tips_[byte_idx] &= ~(0x3 << bit_offset);
    tips_[byte_idx] |= (static_cast<uint8_t>(tip) << bit_offset);
}

BitTipArrayCodec::TipType BitTipArrayCodec::getTip(uint32_t sample_pos) const {
    size_t byte_idx = sample_pos / 4;
    size_t bit_offset = (sample_pos % 4) * 2;

    if (byte_idx >= tips_.size()) {
        return TIP_GENERAL;
    }

    return static_cast<TipType>((tips_[byte_idx] >> bit_offset) & 0x3);
}

uint32_t BitTipArrayCodec::getOrCreateDictId(const std::vector<int64_t>& arr) {
    std::string key = arrayToString(arr);

    auto it = dict_index_.find(key);
    if (it != dict_index_.end()) {
        return it->second;
    }

    uint32_t id = static_cast<uint32_t>(dictionary_.size());
    dict_index_[key] = id;

    DictEntry entry;
    entry.values = arr;

    // Determine element byte size
    int64_t max_val = 0;
    for (int64_t v : arr) {
        if (v != INT64_MIN) {
            max_val = std::max(max_val, std::abs(v));
        }
    }

    if (max_val < 0xFF) {
        entry.element_bytes = 1;
    } else if (max_val < 0xFFFF) {
        entry.element_bytes = 2;
    } else {
        entry.element_bytes = 4;
    }

    dictionary_.push_back(std::move(entry));
    return id;
}

std::string BitTipArrayCodec::arrayToString(const std::vector<int64_t>& arr) const {
    std::ostringstream oss;
    for (size_t i = 0; i < arr.size(); ++i) {
        if (i > 0) oss << ',';
        if (arr[i] == INT64_MIN) {
            oss << '.';
        } else {
            oss << arr[i];
        }
    }
    return oss.str();
}

void BitTipArrayCodec::reconstructArray(uint32_t sample_pos, TipType tip,
                                         uint32_t& special1_idx,
                                         uint32_t& special2_idx,
                                         uint32_t& dict_idx,
                                         std::vector<int64_t>& out) const {
    out.clear();
    out.resize(expected_len_, 0);

    switch (tip) {
        case TIP_ALL_ZERO:
            // All zeros, already initialized
            break;

        case TIP_SPECIAL_1:
            if (field_type_ == FieldType::AD_LIKE) {
                // sum == arr[0] == 2
                if (!out.empty()) {
                    out[0] = 2;
                }
            } else if (field_type_ == FieldType::PL_LIKE) {
                // AB pattern: arr[0]=0, arr[1]=a, arr[2]=b
                if (special1_idx < special1_payload_.size() &&
                    special2_idx < special2_payload_.size()) {
                    int64_t a = special1_payload_[special1_idx];
                    int64_t b = special2_payload_[special2_idx];
                    if (out.size() >= 3) {
                        out[0] = 0;
                        out[1] = a;
                        out[2] = b;
                        // Fill remaining based on pattern (simplified)
                    }
                }
            }
            break;

        case TIP_SPECIAL_2:
            if (field_type_ == FieldType::AD_LIKE) {
                // sum == arr[0], store sum
                if (special1_idx < special1_payload_.size() && !out.empty()) {
                    out[0] = special1_payload_[special1_idx];
                }
            } else if (field_type_ == FieldType::PL_LIKE) {
                // A*15 pattern
                if (special1_idx < special1_payload_.size()) {
                    int64_t a = special1_payload_[special1_idx];
                    int64_t b = 15 * a;
                    if (out.size() >= 3) {
                        out[0] = 0;
                        out[1] = a;
                        out[2] = b;
                    }
                }
            }
            break;

        case TIP_GENERAL:
            if (dict_idx < dict_ids_.size()) {
                uint32_t id = dict_ids_[dict_idx];
                if (id != UINT32_MAX && id < dictionary_.size()) {
                    out = dictionary_[id].values;
                }
            }
            break;
    }
}

void BitTipArrayCodec::serializeDictionary(std::vector<uint8_t>& out) const {
    vint_code::WriteVint(static_cast<uint32_t>(dictionary_.size()), out);

    for (const auto& entry : dictionary_) {
        out.push_back(entry.element_bytes);
        vint_code::WriteVint(static_cast<uint32_t>(entry.values.size()), out);

        for (int64_t v : entry.values) {
            // ZigZag encode for signed values
            int64_t zz = (v << 1) ^ (v >> 63);
            switch (entry.element_bytes) {
                case 1:
                    out.push_back(static_cast<uint8_t>(zz));
                    break;
                case 2:
                    out.push_back(static_cast<uint8_t>(zz & 0xFF));
                    out.push_back(static_cast<uint8_t>((zz >> 8) & 0xFF));
                    break;
                case 4:
                default:
                    out.push_back(static_cast<uint8_t>(zz & 0xFF));
                    out.push_back(static_cast<uint8_t>((zz >> 8) & 0xFF));
                    out.push_back(static_cast<uint8_t>((zz >> 16) & 0xFF));
                    out.push_back(static_cast<uint8_t>((zz >> 24) & 0xFF));
                    break;
            }
        }
    }
}

size_t BitTipArrayCodec::deserializeDictionary(const uint8_t* data, size_t len) {
    size_t pos = 0;

    uint32_t dict_count = vint_code::ReadVint(data, len, pos);
    dictionary_.resize(dict_count);

    for (uint32_t i = 0; i < dict_count; ++i) {
        DictEntry& entry = dictionary_[i];
        entry.element_bytes = data[pos++];

        uint32_t val_count = vint_code::ReadVint(data, len, pos);
        entry.values.resize(val_count);

        for (uint32_t j = 0; j < val_count; ++j) {
            int64_t zz = 0;
            switch (entry.element_bytes) {
                case 1:
                    zz = data[pos++];
                    break;
                case 2:
                    zz = data[pos] | (static_cast<int64_t>(data[pos + 1]) << 8);
                    pos += 2;
                    break;
                case 4:
                default:
                    zz = data[pos] |
                         (static_cast<int64_t>(data[pos + 1]) << 8) |
                         (static_cast<int64_t>(data[pos + 2]) << 16) |
                         (static_cast<int64_t>(data[pos + 3]) << 24);
                    pos += 4;
                    break;
            }
            // ZigZag decode
            entry.values[j] = (zz >> 1) ^ -(zz & 1);
        }

        // Rebuild dict_index for this entry
        dict_index_[arrayToString(entry.values)] = i;
    }

    return pos;
}

} // namespace gsc
