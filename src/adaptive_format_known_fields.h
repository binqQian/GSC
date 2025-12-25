#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/vcf.h"
#include "variant.h"

namespace gsc {

// Adaptive FORMAT stream payload for known heavy fields:
// AD, DP, PL, GQ, PGT, PID.
//
// Row payload layout (V1):
//   "FMD1" (4 bytes)
//   allele_count (u16 LE)
//   presence_mask (u8):
//     bit0 AD, bit1 DP, bit2 PL, bit3 GQ, bit4 PGT, bit5 PID
//   Then, for each present field in fixed order AD,DP,PL,GQ,PGT,PID:
//     AD:
//       tips_size(vint) + tips_bytes
//       miss_count(vint) + miss_pos_deltas(vint...)
//       adid_count(vint) + adid_values(vint...)   // sums or AD dict ids for tips 10/11
//     DP:
//       ab_count(vint) + (pos_delta(vint), value(vint))...
//       normal_count(vint) + normal_values(vint)...
//     PL:
//       tips_size(vint) + tips_bytes
//       miss_count(vint) + miss_pos_deltas(vint...)
//       plid_count(vint) + plid_values(vint...)   // dict ids for tip 11
//       a_count(vint) + a_values(vint...)         // a for tip 01/10
//       b_count(vint) + b_values(vint...)         // b for tip 01 only
//     GQ:
//       ab_count(vint) + (pos_delta(vint), value(vint))...
//       normal_count(vint) + normal_values(vint)...
//     PGT:
//       pos_count(vint) + pos_deltas(vint...)
//       then pos_count strings: (len(vint) + bytes...)
//     PID:
//       pos_count(vint) + pos_deltas(vint...)
//       pid_id_count(vint) + pid_ids(vint...)     // dict ids for non-missing samples
//
// Dictionaries (optional, global across file):
// Streams:
//   adaptive_format_ad_dict, adaptive_format_pl_dict, adaptive_format_pid_dict
// Each stream data:
//   count(vint) + repeat count times: len(vint) + bytes[len]

struct AdaptiveKnownFieldsIndices
{
    int ad = -1;
    int dp = -1;
    int pl = -1;
    int gq = -1;
    int pgt = -1;
    int pid = -1;

    bool Any() const { return ad >= 0 || dp >= 0 || pl >= 0 || gq >= 0 || pgt >= 0 || pid >= 0; }
    bool IsKnownFieldIndex(int idx) const
    {
        return idx == ad || idx == dp || idx == pl || idx == gq || idx == pgt || idx == pid;
    }
};

struct AdaptiveKnownFieldsDicts
{
    std::unordered_map<std::string, uint32_t> ad_map;
    std::unordered_map<std::string, uint32_t> pl_map;
    std::unordered_map<std::string, uint32_t> pid_map;

    std::vector<std::string> ad_items;
    std::vector<std::string> pl_items;
    std::vector<std::string> pid_items;

    uint32_t GetOrAddAD(const std::vector<uint32_t>& values, uint32_t sum_hint);
    uint32_t GetOrAddPL(const std::vector<uint32_t>& values, uint32_t max_hint);
    uint32_t GetOrAddPID(const char* s, size_t len);

    void Serialize(std::vector<uint8_t>& out_ad,
                   std::vector<uint8_t>& out_pl,
                   std::vector<uint8_t>& out_pid) const;
    static bool Deserialize(const std::vector<uint8_t>& in,
                            std::vector<std::string>& out_items);
};

struct AdaptiveKnownFieldsRowV1
{
    static constexpr uint8_t kMagic[4] = {'F', 'M', 'D', '1'};
    static constexpr uint8_t kMaskAD = 1u << 0;
    static constexpr uint8_t kMaskDP = 1u << 1;
    static constexpr uint8_t kMaskPL = 1u << 2;
    static constexpr uint8_t kMaskGQ = 1u << 3;
    static constexpr uint8_t kMaskPGT = 1u << 4;
    static constexpr uint8_t kMaskPID = 1u << 5;

    uint16_t allele_count = 2;
    uint8_t presence_mask = 0;

    // AD
    std::vector<uint8_t> ad_tips;
    std::vector<uint32_t> ad_miss_pos;
    std::vector<uint32_t> ad_ids; // for tip 10/11

    // DP
    std::vector<uint32_t> dp_ab_pos;
    std::vector<uint32_t> dp_ab_val;
    std::vector<uint32_t> dp_normal;

    // PL
    std::vector<uint8_t> pl_tips;
    std::vector<uint32_t> pl_miss_pos;
    std::vector<uint32_t> pl_ids; // for tip 11
    std::vector<uint32_t> pl_a;   // for tip 01/10
    std::vector<uint32_t> pl_b;   // for tip 01

    // GQ
    std::vector<uint32_t> gq_ab_pos;
    std::vector<uint32_t> gq_ab_val;
    std::vector<uint32_t> gq_normal;

    // PGT
    std::vector<uint32_t> pgt_pos;
    std::vector<std::string> pgt_vals; // aligned with pgt_pos

    // PID
    std::vector<uint32_t> pid_pos;
    std::vector<uint32_t> pid_ids; // aligned with pid_pos

    bool Parse(const uint8_t* data, size_t len);

    bool Has(uint8_t mask) const { return (presence_mask & mask) != 0; }
    uint32_t ADCnt() const { return static_cast<uint32_t>(allele_count); }
    uint32_t PLCnt() const
    {
        const uint32_t a = static_cast<uint32_t>(allele_count);
        return a * (a + 1u) / 2u;
    }

    bool DecodeAD(uint32_t sample_count,
                  const std::vector<std::string>& ad_dict,
                  std::vector<int32_t>& out,
                  std::vector<uint32_t>* out_sum,
                  std::vector<uint8_t>* out_has_sum) const;

    bool DecodeDP(uint32_t sample_count,
                  const std::vector<uint32_t>& ad_sum,
                  const std::vector<uint8_t>& ad_has_sum,
                  std::vector<int32_t>& out) const;

    bool DecodePL(uint32_t sample_count,
                  const std::vector<std::string>& pl_dict,
                  std::vector<int32_t>& out,
                  std::vector<int32_t>* out_pred_gq,
                  std::vector<uint8_t>* out_has_pred) const;

    bool DecodeGQ(uint32_t sample_count,
                  const std::vector<int32_t>& pl_pred_gq,
                  const std::vector<uint8_t>& pl_has_pred,
                  std::vector<int32_t>& out) const;

    bool DecodePGT(uint32_t sample_count, std::vector<std::string>& out) const;
    bool DecodePID(uint32_t sample_count,
                   const std::vector<std::string>& pid_dict,
                   std::vector<std::string>& out) const;
};

// Encoder entry point.
bool EncodeAdaptiveKnownFieldsRowV1(const bcf1_t* rec,
                                   uint32_t sample_count,
                                   const std::vector<field_desc>& fields,
                                   const AdaptiveKnownFieldsIndices& idx,
                                   AdaptiveKnownFieldsDicts& dicts,
                                   std::vector<uint8_t>& out_row);

} // namespace gsc

