#pragma once

#include "../format_field_codec.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

namespace gsc {

// ============================================================================
// SparseDictCodec - For fields with high missing rate like PGT/PID
//
// Encoding strategy:
// - Assume most values are missing (".")
// - Only store non-missing positions and their dictionary IDs
// - String dictionary for unique values
//
// Format:
//   [CodecParams]
//   [sample_count: vint]
//   [non_missing_count: vint]
//   [non_missing_positions: delta-encoded vint array]
//   [dict_ids: vint array (one per non-missing)]
//   [dictionary: string array]
// ============================================================================

class SparseDictCodec : public FormatFieldCodec {
public:
    SparseDictCodec();

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

    CodecType type() const override { return CodecType::SparseDictString; }
    CodecParams getParams() const override;
    void reset() override;

private:
    bool frozen_ = false;
    uint32_t sample_count_ = 0;

    // Non-missing data
    std::vector<uint32_t> non_missing_positions_;
    std::vector<uint32_t> dict_ids_;

    // String dictionary
    std::vector<std::string> dictionary_;
    std::unordered_map<std::string, uint32_t> dict_index_;

    // For decoding: quick lookup
    std::unordered_map<uint32_t, uint32_t> position_to_idx_;  // pos -> index in dict_ids_

    // Observation phase
    struct ObservedValue {
        uint32_t position;
        std::string value;
    };
    std::vector<ObservedValue> observed_values_;

    // Helper functions
    bool isMissing(const char* value, size_t len) const;
    uint32_t getOrCreateDictId(const std::string& value);
};

} // namespace gsc
