#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <memory>
#include <cstring>
#include "fmt_utils.h"
#include "fmt_dictionaries.h"
#include "vint_codec.h"

namespace fmt_compress {

// Special FORMAT field names that use optimized compression
constexpr const char* kSpecialFields[] = {"AD", "DP", "GQ", "PL", "PGT", "PID"};

// FMT field processor for optimized compression of AD/DP/PL/GQ/PGT/PID fields
class FmtFieldProcessor {
public:
    explicit FmtFieldProcessor(FmtDictionaries* dictionaries)
        : dictionaries_(dictionaries) {}

    ~FmtFieldProcessor() = default;

    // Check if a field name requires special processing
    static bool isSpecialField(const std::string& name) {
        for (const auto& f : kSpecialFields) {
            if (name == f) return true;
        }
        return false;
    }

    // Clear block-level data for next block
    void clear() {
        ad_tips_.clear();
        ad_ids_.clear();
        dp_exceptions_.clear();
        dp_normal_.clear();

        pl_tips_.clear();
        pl_ids_.clear();
        pl_a_values_.clear();
        pl_b_values_.clear();
        gq_exceptions_.clear();

        pgt_positions_.clear();
        pgt_data_.clear();
        pid_positions_.clear();
        pid_ids_.clear();

        sample_count_ = 0;
    }

    // Process AD and DP fields together (they are correlated)
    // AD: 2-bit tip encoding
    //   00: sum == 0 (all zeros)
    //   01: sum == AD[0] && sum == 2 (special single value)
    //   10: sum == AD[0] (single value, store sum)
    //   11: general case (store dict ID)
    // DP: predict from sum(AD), only store exceptions
    void processADDP(const char* ad_data, uint32_t ad_len,
                     uint32_t dp_value, uint32_t allele_count, uint32_t pos) {
        if (!ad_data || ad_len == 0) {
            // No AD, store DP directly if present
            if (dp_value != kMissingValue16) {
                dp_normal_.push_back(dp_value);
            }
            return;
        }

        // Parse AD values
        ad_vec_.clear();
        splitString2uint(ad_data, ad_len, ',', ad_vec_);

        // Calculate sum
        uint32_t sum = sumArray(static_cast<uint32_t>(ad_vec_.size()), ad_vec_.data());

        // Check DP exception (DP should equal sum(AD))
        if (dp_value != kMissingValue16 && dp_value != sum) {
            dp_exceptions_.emplace_back(dp_value, pos);
        }

        // Encode AD based on pattern
        if (sum == 0) {
            // Type 00: all zeros
            ad_tips_.push_back(0);
            ad_tips_.push_back(0);
        } else if (sum == ad_vec_[0]) {
            // Only first value is non-zero
            if (sum == 2) {
                // Type 01: common case sum==2
                ad_tips_.push_back(0);
                ad_tips_.push_back(1);
            } else {
                // Type 10: store single value
                ad_tips_.push_back(1);
                ad_tips_.push_back(0);
                ad_ids_.push_back(sum);
            }
        } else {
            // Type 11: general case, use dictionary
            ad_tips_.push_back(1);
            ad_tips_.push_back(1);

            ADItem item;
            buildADItem(sum, allele_count, item);
            uint32_t id = dictionaries_->getADItemId(item);
            ad_ids_.push_back(id);
        }
    }

    // Process PL and GQ fields together (they are correlated)
    // PL: 4 pattern types
    //   Type 0 (00): all zeros
    //   Type 1 (01): [0, a, b, a, b, b, ...] pattern, store a, b
    //   Type 2 (10): Type 1 with b == 15*a, store only a
    //   Type 3 (11): no pattern, use dictionary
    // GQ: predict from PL, only store exceptions
    void processPLGQ(const char* pl_data, uint32_t pl_len,
                     uint32_t gq_value, uint32_t pl_count, uint32_t pos) {
        if (!pl_data || pl_len == 0) {
            return;
        }

        // Parse PL values
        pl_vec_.clear();
        splitString2uint(pl_data, pl_len, ',', pl_vec_);

        uint32_t a_val, b_val;
        uint8_t type = checkPlPattern(pl_vec_.data(),
                                      static_cast<uint32_t>(pl_vec_.size()),
                                      a_val, b_val);

        // Encode based on type
        switch (type) {
            case 0:
                // All zeros
                pl_tips_.push_back(0);
                pl_tips_.push_back(0);
                // GQ should be missing for all-zero PL
                if (gq_value != kMissingValue16 && gq_value != 0) {
                    gq_exceptions_.emplace_back(gq_value, pos);
                }
                break;

            case 1:
                // Standard pattern [0, a, b, ...]
                pl_tips_.push_back(0);
                pl_tips_.push_back(1);
                pl_a_values_.push_back(a_val);
                pl_b_values_.push_back(b_val);
                // GQ should equal a_val
                if (gq_value != kMissingValue16 && gq_value != a_val) {
                    gq_exceptions_.emplace_back(gq_value, pos);
                }
                break;

            case 2:
                // Special pattern b == 15*a
                pl_tips_.push_back(1);
                pl_tips_.push_back(0);
                pl_a_values_.push_back(a_val);
                // GQ should equal a_val
                if (gq_value != kMissingValue16 && gq_value != a_val) {
                    gq_exceptions_.emplace_back(gq_value, pos);
                }
                break;

            case 3:
            default:
                // No pattern, use dictionary
                pl_tips_.push_back(1);
                pl_tips_.push_back(1);
                {
                    PLItem item;
                    buildPLItem(pl_vec_, item);
                    uint32_t id = dictionaries_->getPLItemId(item);
                    pl_ids_.push_back(id);

                    // GQ should equal second_min of PL
                    uint32_t expected_gq = getSecondMin(pl_vec_.data(),
                                                        static_cast<uint32_t>(pl_vec_.size()));
                    if (gq_value != kMissingValue16 && gq_value != expected_gq) {
                        gq_exceptions_.emplace_back(gq_value, pos);
                    }
                }
                break;
        }
    }

    // Process PGT field (sparse storage for non-"." values)
    void processPGT(const char* data, uint32_t len, uint32_t pos) {
        if (!data || len == 0 || (len == 1 && data[0] == '.')) {
            return;
        }

        // Store position and data
        pgt_positions_.push_back(pos);
        pgt_data_.push_back(static_cast<uint8_t>(len));
        pgt_data_.insert(pgt_data_.end(), data, data + len);
    }

    // Process PID field (sparse storage with dictionary)
    void processPID(const char* data, uint32_t len, uint32_t pos) {
        if (!data || len == 0 || (len == 1 && data[0] == '.')) {
            return;
        }

        // Store position
        pid_positions_.push_back(pos);

        // Build and store PID item
        PIDItem item;
        memset(item.data_, 0, sizeof(item.data_));
        memcpy(item.data_, &len, 2);
        memcpy(item.data_ + 2, data, std::min(static_cast<uint32_t>(len), 510u));

        uint32_t id = dictionaries_->getPIDItemId(item);
        pid_ids_.push_back(id);
    }

    // Serialize compressed data to output buffer
    void serialize(std::vector<uint8_t>& output) const {
        uint8_t buf[16];

        // AD section
        serializeTips(ad_tips_, output);
        serializeVintVector(ad_ids_, output);

        // DP section
        serializeAbnormalItems(dp_exceptions_, output);
        serializeVintVector(dp_normal_, output);

        // PL section
        serializeTips(pl_tips_, output);
        serializeVintVector(pl_ids_, output);
        serializeVintVector(pl_a_values_, output);
        serializeVintVector(pl_b_values_, output);

        // GQ section
        serializeAbnormalItems(gq_exceptions_, output);

        // PGT section
        serializeVintVector(pgt_positions_, output);
        uint8_t len = VintCodec::encode(pgt_data_.size(), buf);
        output.insert(output.end(), buf, buf + len);
        output.insert(output.end(), pgt_data_.begin(), pgt_data_.end());

        // PID section
        serializeVintVector(pid_positions_, output);
        serializeVintVector(pid_ids_, output);
    }

    // Deserialize compressed data from buffer
    size_t unserialize(const uint8_t* data, size_t data_len) {
        size_t offset = 0;

        // AD section
        offset += unserializeTips(data + offset, ad_tips_);
        offset += unserializeVintVector(data + offset, ad_ids_);

        // DP section
        offset += unserializeAbnormalItems(data + offset, dp_exceptions_);
        offset += unserializeVintVector(data + offset, dp_normal_);

        // PL section
        offset += unserializeTips(data + offset, pl_tips_);
        offset += unserializeVintVector(data + offset, pl_ids_);
        offset += unserializeVintVector(data + offset, pl_a_values_);
        offset += unserializeVintVector(data + offset, pl_b_values_);

        // GQ section
        offset += unserializeAbnormalItems(data + offset, gq_exceptions_);

        // PGT section
        offset += unserializeVintVector(data + offset, pgt_positions_);
        uint64_t pgt_size;
        offset += VintCodec::decode(data + offset, pgt_size);
        pgt_data_.assign(data + offset, data + offset + pgt_size);
        offset += pgt_size;

        // PID section
        offset += unserializeVintVector(data + offset, pid_positions_);
        offset += unserializeVintVector(data + offset, pid_ids_);

        // Build lookup sets for fast position checking
        buildPositionSets();

        return offset;
    }

    // Reconstruct AD field at given position
    bool reconstructAD(uint32_t tip_idx, uint32_t& id_idx,
                       uint8_t allele_count, std::string& output) const {
        if (tip_idx * 2 + 1 >= ad_tips_.size()) return false;

        uint8_t tip0 = ad_tips_[tip_idx * 2];
        uint8_t tip1 = ad_tips_[tip_idx * 2 + 1];

        if (tip0 == 0 && tip1 == 0) {
            // All zeros
            output = "0";
            for (uint8_t i = 1; i < allele_count; i++) {
                output += ",0";
            }
        } else if (tip0 == 0 && tip1 == 1) {
            // Single value = 2
            output = "2";
            for (uint8_t i = 1; i < allele_count; i++) {
                output += ",0";
            }
        } else if (tip0 == 1 && tip1 == 0) {
            // Single value stored
            uint32_t val = ad_ids_[id_idx++];
            output = std::to_string(val);
            for (uint8_t i = 1; i < allele_count; i++) {
                output += ",0";
            }
        } else {
            // Dictionary lookup
            uint32_t id = ad_ids_[id_idx++];
            reconstructADFromDict(id, allele_count, output);
        }

        return true;
    }

    // Reconstruct DP field at given position
    bool reconstructDP(uint32_t pos, uint32_t tip_idx,
                       std::string& output) const {
        // Check exceptions first
        auto it = dp_exception_set_.find(pos);
        if (it != dp_exception_set_.end()) {
            // Find exception value
            for (const auto& exc : dp_exceptions_) {
                if (exc.pos == pos) {
                    output = std::to_string(exc.val);
                    return true;
                }
            }
        }

        // Use predicted value (sum of AD)
        if (tip_idx * 2 + 1 >= ad_tips_.size()) {
            if (!dp_normal_.empty()) {
                // No AD, use normal DP
                // Note: need proper index tracking
            }
            return false;
        }

        // Reconstruct sum from AD tips
        uint8_t tip0 = ad_tips_[tip_idx * 2];
        uint8_t tip1 = ad_tips_[tip_idx * 2 + 1];

        if (tip0 == 0 && tip1 == 0) {
            output = "0";
        } else if (tip0 == 0 && tip1 == 1) {
            output = "2";
        } else {
            // Need to track AD id index for sum calculation
            // For now, return empty
            return false;
        }

        return true;
    }

    // Reconstruct PL field at given position
    bool reconstructPL(uint32_t tip_idx, uint32_t& id_idx,
                       uint32_t& a_idx, uint32_t& b_idx,
                       uint8_t allele_count, std::string& output) const {
        if (tip_idx * 2 + 1 >= pl_tips_.size()) return false;

        uint8_t tip0 = pl_tips_[tip_idx * 2];
        uint8_t tip1 = pl_tips_[tip_idx * 2 + 1];
        uint16_t pl_cnt = kPlCountTable[allele_count];

        if (tip0 == 0 && tip1 == 0) {
            // All zeros
            output = "0";
            for (uint16_t i = 1; i < pl_cnt; i++) {
                output += ",0";
            }
        } else if (tip0 == 0 && tip1 == 1) {
            // Pattern with a, b
            uint32_t a = pl_a_values_[a_idx++];
            uint32_t b = pl_b_values_[b_idx++];
            patternToPlString(a, b, allele_count, pl_cnt, output);
        } else if (tip0 == 1 && tip1 == 0) {
            // Pattern with a, b = 15*a
            uint32_t a = pl_a_values_[a_idx++];
            uint32_t b = a * 15;
            patternToPlString(a, b, allele_count, pl_cnt, output);
        } else {
            // Dictionary lookup
            uint32_t id = pl_ids_[id_idx++];
            reconstructPLFromDict(id, output);
        }

        return true;
    }

    // Reconstruct GQ field at given position
    bool reconstructGQ(uint32_t pos, uint32_t tip_idx,
                       uint32_t a_idx, std::string& output) const {
        // Check exceptions first
        auto it = gq_exception_set_.find(pos);
        if (it != gq_exception_set_.end()) {
            for (const auto& exc : gq_exceptions_) {
                if (exc.pos == pos) {
                    output = std::to_string(exc.val);
                    return true;
                }
            }
        }

        // Use predicted value from PL
        if (tip_idx * 2 + 1 >= pl_tips_.size()) return false;

        uint8_t tip0 = pl_tips_[tip_idx * 2];
        uint8_t tip1 = pl_tips_[tip_idx * 2 + 1];

        if (tip0 == 0 && tip1 == 0) {
            // PL all zeros -> GQ missing or 0
            output = ".";
        } else if ((tip0 == 0 && tip1 == 1) || (tip0 == 1 && tip1 == 0)) {
            // GQ = a_val
            if (a_idx < pl_a_values_.size()) {
                output = std::to_string(pl_a_values_[a_idx]);
            }
        } else {
            // Need to calculate second_min from PL dict
            // For now, return missing
            output = ".";
        }

        return true;
    }

    // Reconstruct PGT field at given position
    bool reconstructPGT(uint32_t pos, std::string& output) const {
        auto it = pgt_position_set_.find(pos);
        if (it == pgt_position_set_.end()) {
            output = ".";
            return true;
        }

        // Find data in pgt_data_
        size_t data_offset = 0;
        for (size_t i = 0; i < pgt_positions_.size(); i++) {
            if (pgt_positions_[i] == pos) {
                uint8_t len = pgt_data_[data_offset];
                output.assign(reinterpret_cast<const char*>(pgt_data_.data() + data_offset + 1), len);
                return true;
            }
            data_offset += 1 + pgt_data_[data_offset];
        }

        output = ".";
        return true;
    }

    // Reconstruct PID field at given position
    bool reconstructPID(uint32_t pos, std::string& output) const {
        auto it = pid_position_set_.find(pos);
        if (it == pid_position_set_.end()) {
            output = ".";
            return true;
        }

        // Find ID in pid_ids_
        for (size_t i = 0; i < pid_positions_.size(); i++) {
            if (pid_positions_[i] == pos) {
                uint32_t id = pid_ids_[i];
                const uint8_t* ptr = dictionaries_->getPIDItemPtr(id);
                if (ptr) {
                    uint16_t len;
                    memcpy(&len, ptr, 2);
                    output.assign(reinterpret_cast<const char*>(ptr + 2), len);
                }
                return true;
            }
        }

        output = ".";
        return true;
    }

    void setSampleCount(uint32_t count) { sample_count_ = count; }
    uint32_t getSampleCount() const { return sample_count_; }

private:
    // Build AD dictionary item
    void buildADItem(uint32_t sum, uint32_t allele_count, ADItem& item) const {
        memset(item.data_, 0, sizeof(item.data_));

        // Determine type based on value range
        uint8_t type;
        if (sum < 0xFF) {
            type = 0;  // 1 byte per value
        } else if (sum < 0xFFFF) {
            type = 1;  // 2 bytes per value
        } else {
            type = 2;  // 4 bytes per value
        }

        item.data_[1] = type;

        // Store values
        uint8_t* ptr = item.data_ + 2;
        for (size_t i = 0; i < ad_vec_.size() && ptr < item.data_ + 47; i++) {
            uint32_t val = ad_vec_[i];
            switch (type) {
                case 0:
                    *ptr++ = static_cast<uint8_t>(val);
                    break;
                case 1:
                    *reinterpret_cast<uint16_t*>(ptr) = static_cast<uint16_t>(val);
                    ptr += 2;
                    break;
                case 2:
                    *reinterpret_cast<uint32_t*>(ptr) = val;
                    ptr += 4;
                    break;
            }
        }

        item.data_[0] = static_cast<uint8_t>(ptr - item.data_ - 1);
    }

    // Build PL dictionary item
    void buildPLItem(const std::vector<uint32_t>& pl_values, PLItem& item) const {
        memset(item.data_, 0, sizeof(item.data_));

        uint32_t max_val = 0;
        for (auto v : pl_values) {
            if (v != kMissingValue && v > max_val) max_val = v;
        }

        uint8_t type;
        if (max_val < 0xFF) {
            type = 0;
        } else if (max_val < 0xFFFF) {
            type = 1;
        } else {
            type = 2;
        }

        item.data_[1] = type;

        uint8_t* ptr = item.data_ + 2;
        for (size_t i = 0; i < pl_values.size() && ptr < item.data_ + 127; i++) {
            uint32_t val = pl_values[i];
            switch (type) {
                case 0:
                    *ptr++ = static_cast<uint8_t>(val);
                    break;
                case 1:
                    *reinterpret_cast<uint16_t*>(ptr) = static_cast<uint16_t>(val);
                    ptr += 2;
                    break;
                case 2:
                    *reinterpret_cast<uint32_t*>(ptr) = val;
                    ptr += 4;
                    break;
            }
        }

        item.data_[0] = static_cast<uint8_t>(ptr - item.data_ - 1);
    }

    // Reconstruct AD from dictionary
    void reconstructADFromDict(uint32_t id, uint8_t allele_count,
                               std::string& output) const {
        const uint8_t* ptr = dictionaries_->getADItemPtr(id);
        if (!ptr) {
            output = ".";
            return;
        }

        // ptr[0] is data length (unused in reconstruction)
        uint8_t type = ptr[1];
        const uint8_t* data = ptr + 2;

        output.clear();
        for (uint8_t i = 0; i < allele_count; i++) {
            if (i > 0) output += ",";

            uint32_t val = 0;
            switch (type) {
                case 0:
                    val = data[i];
                    break;
                case 1:
                    val = reinterpret_cast<const uint16_t*>(data)[i];
                    break;
                case 2:
                    val = reinterpret_cast<const uint32_t*>(data)[i];
                    break;
            }
            output += std::to_string(val);
        }
    }

    // Reconstruct PL from dictionary
    void reconstructPLFromDict(uint32_t id, std::string& output) const {
        const uint8_t* ptr = dictionaries_->getPLItemPtr(id);
        if (!ptr) {
            output = ".";
            return;
        }

        uint8_t len = ptr[0];
        uint8_t type = ptr[1];
        const uint8_t* data = ptr + 2;

        output.clear();
        size_t count = 0;
        switch (type) {
            case 0: count = len - 1; break;
            case 1: count = (len - 1) / 2; break;
            case 2: count = (len - 1) / 4; break;
        }

        for (size_t i = 0; i < count; i++) {
            if (i > 0) output += ",";

            uint32_t val = 0;
            switch (type) {
                case 0:
                    val = data[i];
                    break;
                case 1:
                    val = reinterpret_cast<const uint16_t*>(data)[i];
                    break;
                case 2:
                    val = reinterpret_cast<const uint32_t*>(data)[i];
                    break;
            }
            output += std::to_string(val);
        }
    }

    // Serialize helper: tips (2-bit packed)
    void serializeTips(const std::vector<uint8_t>& tips,
                       std::vector<uint8_t>& output) const {
        uint8_t buf[16];
        uint8_t len = VintCodec::encode(tips.size(), buf);
        output.insert(output.end(), buf, buf + len);

        if (!tips.empty()) {
            size_t packed_size = (tips.size() + 7) / 8;
            size_t start = output.size();
            output.resize(start + packed_size, 0);
            packTips(tips, output.data() + start);
        }
    }

    // Serialize helper: vint vector
    void serializeVintVector(const std::vector<uint32_t>& vec,
                             std::vector<uint8_t>& output) const {
        uint8_t buf[16];
        uint8_t len = VintCodec::encode(vec.size(), buf);
        output.insert(output.end(), buf, buf + len);

        for (uint32_t val : vec) {
            len = VintCodec::encode(val, buf);
            output.insert(output.end(), buf, buf + len);
        }
    }

    // Serialize helper: abnormal items
    void serializeAbnormalItems(const std::vector<AbnormalItem>& items,
                                std::vector<uint8_t>& output) const {
        uint8_t buf[16];
        uint8_t len = VintCodec::encode(items.size(), buf);
        output.insert(output.end(), buf, buf + len);

        for (const auto& item : items) {
            len = VintCodec::encode(item.val, buf);
            output.insert(output.end(), buf, buf + len);
            len = VintCodec::encode(item.pos, buf);
            output.insert(output.end(), buf, buf + len);
        }
    }

    // Unserialize helper: tips
    size_t unserializeTips(const uint8_t* data, std::vector<uint8_t>& tips) {
        size_t offset = 0;
        uint64_t count;
        offset += VintCodec::decode(data, count);

        tips.clear();
        if (count > 0) {
            unpackTips(data + offset, static_cast<uint32_t>(count), tips);
            offset += (count + 7) / 8;
        }

        return offset;
    }

    // Unserialize helper: vint vector
    size_t unserializeVintVector(const uint8_t* data, std::vector<uint32_t>& vec) {
        size_t offset = 0;
        uint64_t count;
        offset += VintCodec::decode(data, count);

        vec.clear();
        vec.reserve(count);
        for (uint64_t i = 0; i < count; i++) {
            uint64_t val;
            offset += VintCodec::decode(data + offset, val);
            vec.push_back(static_cast<uint32_t>(val));
        }

        return offset;
    }

    // Unserialize helper: abnormal items
    size_t unserializeAbnormalItems(const uint8_t* data,
                                    std::vector<AbnormalItem>& items) {
        size_t offset = 0;
        uint64_t count;
        offset += VintCodec::decode(data, count);

        items.clear();
        items.reserve(count);
        for (uint64_t i = 0; i < count; i++) {
            uint64_t val, pos;
            offset += VintCodec::decode(data + offset, val);
            offset += VintCodec::decode(data + offset, pos);
            items.emplace_back(static_cast<uint32_t>(val), static_cast<uint32_t>(pos));
        }

        return offset;
    }

    // Build position sets for fast lookup
    void buildPositionSets() {
        dp_exception_set_.clear();
        for (const auto& item : dp_exceptions_) {
            dp_exception_set_.insert(item.pos);
        }

        gq_exception_set_.clear();
        for (const auto& item : gq_exceptions_) {
            gq_exception_set_.insert(item.pos);
        }

        pgt_position_set_.clear();
        for (uint32_t pos : pgt_positions_) {
            pgt_position_set_.insert(pos);
        }

        pid_position_set_.clear();
        for (uint32_t pos : pid_positions_) {
            pid_position_set_.insert(pos);
        }
    }

    // Global dictionaries (shared across all processors)
    FmtDictionaries* dictionaries_;

    // Block-level data
    uint32_t sample_count_ = 0;

    // AD data
    std::vector<uint8_t> ad_tips_;      // 2-bit per sample
    std::vector<uint32_t> ad_ids_;      // Dictionary IDs or single values
    std::vector<uint32_t> ad_vec_;      // Temp: parsed AD values

    // DP data
    std::vector<AbnormalItem> dp_exceptions_;
    std::vector<uint32_t> dp_normal_;   // DP values when no AD present

    // PL data
    std::vector<uint8_t> pl_tips_;      // 2-bit per sample
    std::vector<uint32_t> pl_ids_;      // Dictionary IDs
    std::vector<uint32_t> pl_a_values_; // Pattern a values
    std::vector<uint32_t> pl_b_values_; // Pattern b values
    std::vector<uint32_t> pl_vec_;      // Temp: parsed PL values

    // GQ data
    std::vector<AbnormalItem> gq_exceptions_;

    // PGT data
    std::vector<uint32_t> pgt_positions_;
    std::vector<uint8_t> pgt_data_;     // len + raw bytes

    // PID data
    std::vector<uint32_t> pid_positions_;
    std::vector<uint32_t> pid_ids_;

    // Position lookup sets (for decompression)
    std::unordered_set<uint32_t> dp_exception_set_;
    std::unordered_set<uint32_t> gq_exception_set_;
    std::unordered_set<uint32_t> pgt_position_set_;
    std::unordered_set<uint32_t> pid_position_set_;
};

}  // namespace fmt_compress
