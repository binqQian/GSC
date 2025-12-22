#pragma once

#include "../format_field_codec.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

namespace gsc {

// ============================================================================
// BitTipArrayCodec - For fixed-length numeric arrays like AD and PL
//
// Uses 2-bit "tips" per sample to indicate encoding type:
//
// For AD-like fields (sum-based patterns):
//   00 = all zeros
//   01 = sum == arr[0] && sum == 2 (common case, no payload)
//   10 = sum == arr[0] && sum != 2 (only store sum)
//   11 = general case (store in dictionary or inline)
//
// For PL-like fields (pattern-based):
//   00 = all zeros
//   01 = matches AB pattern (store a, b)
//   10 = matches A*15 pattern (store a only, b = 15*a)
//   11 = general case (dictionary)
//
// Format:
//   [CodecParams]
//   [tip_count: vint]
//   [tips: packed 2-bits per sample]
//   [payload for tip 01/10: vint array]
//   [payload for tip 11: dictionary + ids]
// ============================================================================

class BitTipArrayCodec : public FormatFieldCodec {
public:
    // Tip types
    enum TipType : uint8_t {
        TIP_ALL_ZERO = 0,      // 00: All elements are zero
        TIP_SPECIAL_1 = 1,     // 01: AD: sum==arr[0]&&sum==2; PL: AB pattern
        TIP_SPECIAL_2 = 2,     // 10: AD: sum==arr[0]&&sum!=2; PL: A*15 pattern
        TIP_GENERAL = 3        // 11: General case (dictionary)
    };

    // Field type affects encoding strategy
    enum class FieldType {
        AD_LIKE,    // Sum-based patterns
        PL_LIKE,    // AB patterns
        GENERIC     // No special patterns
    };

    BitTipArrayCodec(FieldType field_type = FieldType::GENERIC, uint32_t expected_len = 0);

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

    CodecType type() const override { return CodecType::BitTipArray; }
    CodecParams getParams() const override;
    void reset() override;

    // Get sum of array elements (for DP prediction)
    int64_t getComputedValue(uint32_t sample_pos) const override;

    // Set expected array length
    void setExpectedLen(uint32_t len) { expected_len_ = len; }

private:
    bool frozen_ = false;
    FieldType field_type_;
    uint32_t expected_len_;
    uint32_t sample_count_ = 0;

    // Tips: 2 bits per sample, packed into bytes
    std::vector<uint8_t> tips_;

    // Payloads for different tip types
    std::vector<uint32_t> special1_payload_;  // For tip 01
    std::vector<uint32_t> special2_payload_;  // For tip 10 (and second value for PL AB)

    // Dictionary for general case (tip 11)
    struct DictEntry {
        std::vector<int64_t> values;
        uint8_t element_bytes;  // 1, 2, or 4
    };
    std::vector<DictEntry> dictionary_;
    std::unordered_map<std::string, uint32_t> dict_index_;
    std::vector<uint32_t> dict_ids_;  // Dictionary ID for each tip-11 sample

    // Sum cache for cross-field prediction
    std::vector<int64_t> sum_cache_;

    // Observation data
    struct ObservedArray {
        std::vector<int64_t> values;
        int64_t sum;
    };
    std::vector<ObservedArray> observed_arrays_;

    // Missing sample positions
    std::vector<uint32_t> missing_positions_;
    std::unordered_set<uint32_t> missing_set_;

    // Helper functions
    bool parseArray(const char* value, size_t len, std::vector<int64_t>& out) const;
    TipType classifyAD(const std::vector<int64_t>& arr, int64_t sum) const;
    TipType classifyPL(const std::vector<int64_t>& arr, int64_t& a, int64_t& b) const;
    TipType classifyGeneric(const std::vector<int64_t>& arr) const;

    bool checkAbPattern(const std::vector<int64_t>& arr, int64_t& a, int64_t& b) const;

    void setTip(uint32_t sample_pos, TipType tip);
    TipType getTip(uint32_t sample_pos) const;

    uint32_t getOrCreateDictId(const std::vector<int64_t>& arr);
    std::string arrayToString(const std::vector<int64_t>& arr) const;

    void reconstructArray(uint32_t sample_pos, TipType tip,
                          uint32_t& special1_idx, uint32_t& special2_idx,
                          uint32_t& dict_idx,
                          std::vector<int64_t>& out) const;

    // Serialization helpers
    void serializeDictionary(std::vector<uint8_t>& out) const;
    size_t deserializeDictionary(const uint8_t* data, size_t len);
};

} // namespace gsc
