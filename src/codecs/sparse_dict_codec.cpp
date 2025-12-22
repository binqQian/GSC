#include "sparse_dict_codec.h"
#include "../vint_code.h"
#include <algorithm>

namespace gsc {

SparseDictCodec::SparseDictCodec() {
    observed_values_.reserve(256);
    non_missing_positions_.reserve(64);
    dict_ids_.reserve(64);
}

void SparseDictCodec::observe(const char* value, size_t len, uint32_t sample_pos) {
    if (frozen_) return;

    if (!isMissing(value, len)) {
        ObservedValue ov;
        ov.position = sample_pos;
        ov.value = std::string(value, len);
        observed_values_.push_back(std::move(ov));
    }

    sample_count_ = std::max(sample_count_, sample_pos + 1);
}

void SparseDictCodec::freeze() {
    if (frozen_) return;
    frozen_ = true;

    // Sort by position
    std::sort(observed_values_.begin(), observed_values_.end(),
              [](const ObservedValue& a, const ObservedValue& b) {
                  return a.position < b.position;
              });

    // Build dictionary and encode
    for (const auto& ov : observed_values_) {
        non_missing_positions_.push_back(ov.position);
        dict_ids_.push_back(getOrCreateDictId(ov.value));
    }

    // Build position lookup for decoding
    for (size_t i = 0; i < non_missing_positions_.size(); ++i) {
        position_to_idx_[non_missing_positions_[i]] = static_cast<uint32_t>(i);
    }

    observed_values_.clear();
}

void SparseDictCodec::encode(const char* value, size_t len, uint32_t sample_pos) {
    if (!frozen_) {
        observe(value, len, sample_pos);
        return;
    }

    if (!isMissing(value, len)) {
        non_missing_positions_.push_back(sample_pos);
        std::string val(value, len);
        dict_ids_.push_back(getOrCreateDictId(val));
        position_to_idx_[sample_pos] = static_cast<uint32_t>(dict_ids_.size() - 1);
    }

    sample_count_ = std::max(sample_count_, sample_pos + 1);
}

void SparseDictCodec::serialize(std::vector<uint8_t>& out) const {
    // Header
    CodecParams params = getParams();
    params.serialize(out);

    // Sample count
    vint_code::WriteVint(sample_count_, out);

    // Non-missing positions (delta encoded)
    vint_code::WriteVint(static_cast<uint32_t>(non_missing_positions_.size()), out);
    uint32_t prev_pos = 0;
    for (uint32_t pos : non_missing_positions_) {
        vint_code::WriteVint(pos - prev_pos, out);
        prev_pos = pos;
    }

    // Dictionary IDs
    for (uint32_t id : dict_ids_) {
        vint_code::WriteVint(id, out);
    }

    // Dictionary
    vint_code::WriteVint(static_cast<uint32_t>(dictionary_.size()), out);
    for (const auto& str : dictionary_) {
        vint_code::WriteVint(static_cast<uint32_t>(str.size()), out);
        out.insert(out.end(), str.begin(), str.end());
    }
}

size_t SparseDictCodec::deserialize(const uint8_t* data, size_t len) {
    size_t pos = 0;

    // Sample count
    sample_count_ = vint_code::ReadVint(data, len, pos);

    // Non-missing positions
    uint32_t nm_count = vint_code::ReadVint(data, len, pos);
    non_missing_positions_.resize(nm_count);
    uint32_t cur_pos = 0;
    for (uint32_t i = 0; i < nm_count; ++i) {
        cur_pos += vint_code::ReadVint(data, len, pos);
        non_missing_positions_[i] = cur_pos;
    }

    // Dictionary IDs
    dict_ids_.resize(nm_count);
    for (uint32_t i = 0; i < nm_count; ++i) {
        dict_ids_[i] = vint_code::ReadVint(data, len, pos);
    }

    // Dictionary
    uint32_t dict_count = vint_code::ReadVint(data, len, pos);
    dictionary_.resize(dict_count);
    for (uint32_t i = 0; i < dict_count; ++i) {
        uint32_t str_len = vint_code::ReadVint(data, len, pos);
        dictionary_[i] = std::string(reinterpret_cast<const char*>(data + pos), str_len);
        pos += str_len;
    }

    // Build position lookup
    position_to_idx_.clear();
    for (size_t i = 0; i < non_missing_positions_.size(); ++i) {
        position_to_idx_[non_missing_positions_[i]] = static_cast<uint32_t>(i);
    }

    frozen_ = true;
    return pos;
}

std::string SparseDictCodec::decode(uint32_t sample_pos) const {
    auto it = position_to_idx_.find(sample_pos);
    if (it == position_to_idx_.end()) {
        return ".";  // Missing
    }

    uint32_t idx = it->second;
    if (idx >= dict_ids_.size()) {
        return ".";
    }

    uint32_t dict_id = dict_ids_[idx];
    if (dict_id >= dictionary_.size()) {
        return ".";
    }

    return dictionary_[dict_id];
}

CodecParams SparseDictCodec::getParams() const {
    CodecParams params;
    params.type = CodecType::SparseDictString;
    params.sample_count = sample_count_;
    params.array_len = 0;
    params.element_bytes = 2;  // Dictionary ID size
    params.dict_size = static_cast<uint32_t>(dictionary_.size());
    return params;
}

void SparseDictCodec::reset() {
    frozen_ = false;
    sample_count_ = 0;
    non_missing_positions_.clear();
    dict_ids_.clear();
    dictionary_.clear();
    dict_index_.clear();
    position_to_idx_.clear();
    observed_values_.clear();
}

bool SparseDictCodec::isMissing(const char* value, size_t len) const {
    return (len == 0) || (len == 1 && value[0] == '.');
}

uint32_t SparseDictCodec::getOrCreateDictId(const std::string& value) {
    auto it = dict_index_.find(value);
    if (it != dict_index_.end()) {
        return it->second;
    }

    uint32_t id = static_cast<uint32_t>(dictionary_.size());
    dictionary_.push_back(value);
    dict_index_[value] = id;
    return id;
}

} // namespace gsc
