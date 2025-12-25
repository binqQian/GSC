#include "adaptive_format_known_fields.h"

#include <algorithm>
#include <cstring>
#include <limits>

#include "vint_code.h"

namespace gsc {
namespace {

static constexpr uint32_t kMissingU32 = 0xFFFFFFFFu;

static inline void writeU16LE(uint16_t v, std::vector<uint8_t>& out)
{
    out.push_back(static_cast<uint8_t>(v & 0xffu));
    out.push_back(static_cast<uint8_t>((v >> 8) & 0xffu));
}

static inline bool readU16LE(const uint8_t* data, size_t len, size_t& pos, uint16_t& out_v)
{
    if (pos + 2 > len) return false;
    out_v = static_cast<uint16_t>(data[pos] | (static_cast<uint16_t>(data[pos + 1]) << 8));
    pos += 2;
    return true;
}

static inline void setTip2(std::vector<uint8_t>& tips, uint32_t sample_idx, uint8_t tip)
{
    const uint32_t byte = sample_idx >> 2;           // /4
    const uint32_t shift = (sample_idx & 3u) << 1;   // *2
    tips[byte] &= static_cast<uint8_t>(~(0x3u << shift));
    tips[byte] |= static_cast<uint8_t>((tip & 0x3u) << shift);
}

static inline uint8_t getTip2(const std::vector<uint8_t>& tips, uint32_t sample_idx)
{
    const uint32_t byte = sample_idx >> 2;
    const uint32_t shift = (sample_idx & 3u) << 1;
    if (byte >= tips.size()) return 0;
    return static_cast<uint8_t>((tips[byte] >> shift) & 0x3u);
}

static inline bool isMissingScalarInt32(int32_t v)
{
    return v == bcf_int32_missing || v == bcf_int32_vector_end;
}

static inline uint32_t toU32OrMissing(int32_t v)
{
    if (v == bcf_int32_missing || v == bcf_int32_vector_end) return kMissingU32;
    return static_cast<uint32_t>(v);
}

static inline int32_t fromU32OrMissing(uint32_t v)
{
    if (v == kMissingU32) return bcf_int32_missing;
    return static_cast<int32_t>(v);
}

template <typename T>
static void getMaxAndSecondMin(const T* pdata, uint32_t cnt, T& max_val, T& second_min)
{
    max_val = pdata[0];
    T smallest = pdata[0];
    second_min = std::numeric_limits<T>::max();

    for (uint32_t i = 1; i < cnt; ++i)
    {
        if (pdata[i] > max_val) max_val = pdata[i];
        if (pdata[i] < smallest)
        {
            second_min = smallest;
            smallest = pdata[i];
        }
        else if (pdata[i] < second_min && pdata[i] != smallest)
        {
            second_min = pdata[i];
        }
    }
}

// checkAbPattern from ref_code/fmt_comp (operates on non-missing values).
static uint8_t checkAbPattern(const uint32_t* array, uint32_t size, uint32_t& a_val, uint32_t& b_val)
{
    if (size == 0) return 3;
    if (array[0] != 0) return 3;

    bool allZero = true;
    for (uint32_t i = 0; i < size; ++i)
    {
        if (array[i] != 0)
        {
            allZero = false;
            break;
        }
    }
    if (allZero)
    {
        a_val = 0;
        b_val = 0;
        return 0;
    }

    if (size < 3) return 3;
    a_val = array[1];
    b_val = array[2];

    if (size == 3)
    {
        if (a_val * 15u == b_val) return 2;
        return 1;
    }
    if (size == 6)
    {
        if (array[1] == array[3] && array[2] == array[4] && array[2] == array[5])
        {
            if (a_val * 15u == b_val) return 2;
            return 1;
        }
        return 3;
    }

    uint32_t pos = 3;
    uint32_t b_count = 2;
    while (pos < size)
    {
        if (array[pos] != a_val) return 3;
        ++pos;
        for (uint32_t i = 0; i < b_count && pos < size; ++i, ++pos)
        {
            if (array[pos] != b_val) return 3;
        }
        ++b_count;
    }

    if (a_val * 15u == b_val) return 2;
    return 1;
}

static void patternToPLArray(uint32_t a, uint32_t b, uint16_t allele, std::vector<int32_t>& out, size_t off)
{
    const uint32_t ac = static_cast<uint32_t>(allele);
    const uint32_t plcnt = ac * (ac + 1u) / 2u;
    if (off + plcnt > out.size()) return;
    out[off] = 0;
    size_t pos = off + 1;
    for (uint32_t i = 1; i < ac; ++i)
    {
        out[pos++] = static_cast<int32_t>(a);
        for (uint32_t j = 0; j < i; ++j)
            out[pos++] = static_cast<int32_t>(b);
    }
}

static std::string packFixedWidthArrayItem(const std::vector<uint32_t>& values, uint32_t width_bytes)
{
    const uint32_t cnt = static_cast<uint32_t>(values.size());
    const uint32_t payload_len = 1u + cnt * width_bytes; // [type]+values
    std::string out;
    out.resize(1u + payload_len);
    out[0] = static_cast<char>(payload_len);

    uint8_t type = 2;
    if (width_bytes == 1) type = 0;
    else if (width_bytes == 2) type = 1;
    else if (width_bytes == 4) type = 2;
    out[1] = static_cast<char>(type);

    size_t p = 2;
    for (uint32_t v : values)
    {
        uint32_t x = v;
        if (v == kMissingU32)
        {
            x = (width_bytes == 1) ? 0xffu : (width_bytes == 2) ? 0xffffu : 0xffffffffu;
        }

        for (uint32_t b = 0; b < width_bytes; ++b)
        {
            out[p++] = static_cast<char>(x & 0xffu);
            x >>= 8;
        }
    }
    return out;
}

static bool unpackFixedWidthArrayItem(const std::string& item,
                                      uint32_t expected_cnt,
                                      std::vector<uint32_t>& out_vals,
                                      uint32_t& width_bytes)
{
    out_vals.clear();
    width_bytes = 0;
    if (item.size() < 2) return false;
    const uint8_t payload_len = static_cast<uint8_t>(item[0]);
    if (item.size() != static_cast<size_t>(payload_len) + 1) return false;
    const uint8_t type = static_cast<uint8_t>(item[1]);
    if (type == 0) width_bytes = 1;
    else if (type == 1) width_bytes = 2;
    else if (type == 2) width_bytes = 4;
    else return false;

    const size_t need = 2 + static_cast<size_t>(expected_cnt) * width_bytes;
    if (need > item.size()) return false;

    out_vals.resize(expected_cnt);
    size_t p = 2;
    for (uint32_t i = 0; i < expected_cnt; ++i)
    {
        uint32_t x = 0;
        for (uint32_t b = 0; b < width_bytes; ++b)
        {
            x |= (static_cast<uint32_t>(static_cast<uint8_t>(item[p++]))) << (8u * b);
        }
        if ((width_bytes == 1 && x == 0xffu) ||
            (width_bytes == 2 && x == 0xffffu) ||
            (width_bytes == 4 && x == 0xffffffffu))
        {
            out_vals[i] = kMissingU32;
        }
        else
        {
            out_vals[i] = x;
        }
    }
    return true;
}

static bool parseFormatStrings(const field_desc& field, uint32_t sample_count, std::vector<std::string>& out)
{
    out.clear();
    out.resize(sample_count);
    if (!field.present || field.data == nullptr || field.data_size == 0) return false;

    const char* p = field.data;
    const char* end = field.data + field.data_size;
    for (uint32_t s = 0; s < sample_count; ++s)
    {
        if (p >= end) return false;
        size_t len = ::strnlen(p, static_cast<size_t>(end - p));
        out[s].assign(p, len);
        p += len + 1;
    }
    return true;
}

static inline void encodePosListDeltas(const std::vector<uint32_t>& pos, std::vector<uint8_t>& out)
{
    vint_code::WriteVint(static_cast<uint32_t>(pos.size()), out);
    uint32_t prev = 0;
    for (uint32_t p : pos)
    {
        vint_code::WriteVint(p - prev, out);
        prev = p;
    }
}

static inline bool decodePosListDeltas(const uint8_t* data, size_t len, size_t& pos, std::vector<uint32_t>& out_pos)
{
    out_pos.clear();
    if (pos >= len) return false;
    const uint32_t count = vint_code::ReadVint(data, len, pos);
    out_pos.resize(count);
    uint32_t cur = 0;
    for (uint32_t i = 0; i < count; ++i)
    {
        cur += vint_code::ReadVint(data, len, pos);
        out_pos[i] = cur;
    }
    return true;
}

} // namespace

uint32_t AdaptiveKnownFieldsDicts::GetOrAddAD(const std::vector<uint32_t>& values, uint32_t sum_hint)
{
    uint32_t max_v = 0;
    for (uint32_t v : values)
        if (v != kMissingU32) max_v = std::max(max_v, v);
    uint32_t hint = std::max(sum_hint, max_v);
    uint32_t width = (hint < 0xffu) ? 1 : (hint < 0xffffu) ? 2 : 4;
    std::string item = packFixedWidthArrayItem(values, width);
    auto it = ad_map.find(item);
    if (it != ad_map.end()) return it->second;
    uint32_t id = static_cast<uint32_t>(ad_items.size());
    ad_items.push_back(item);
    ad_map.emplace(ad_items.back(), id);
    return id;
}

uint32_t AdaptiveKnownFieldsDicts::GetOrAddPL(const std::vector<uint32_t>& values, uint32_t max_hint)
{
    uint32_t max_v = 0;
    for (uint32_t v : values)
        if (v != kMissingU32) max_v = std::max(max_v, v);
    uint32_t hint = std::max(max_hint, max_v);
    uint32_t width = (hint < 0xffu) ? 1 : (hint < 0xffffu) ? 2 : 4;
    std::string item = packFixedWidthArrayItem(values, width);
    auto it = pl_map.find(item);
    if (it != pl_map.end()) return it->second;
    uint32_t id = static_cast<uint32_t>(pl_items.size());
    pl_items.push_back(item);
    pl_map.emplace(pl_items.back(), id);
    return id;
}

uint32_t AdaptiveKnownFieldsDicts::GetOrAddPID(const char* s, size_t len)
{
    if (s == nullptr) return 0;
    std::string item;
    item.resize(2 + len);
    const uint16_t l16 = static_cast<uint16_t>(len);
    item[0] = static_cast<char>(l16 & 0xffu);
    item[1] = static_cast<char>((l16 >> 8) & 0xffu);
    if (len) std::memcpy(&item[2], s, len);

    auto it = pid_map.find(item);
    if (it != pid_map.end()) return it->second;
    uint32_t id = static_cast<uint32_t>(pid_items.size());
    pid_items.push_back(item);
    pid_map.emplace(pid_items.back(), id);
    return id;
}

void AdaptiveKnownFieldsDicts::Serialize(std::vector<uint8_t>& out_ad,
                                        std::vector<uint8_t>& out_pl,
                                        std::vector<uint8_t>& out_pid) const
{
    auto serializeOne = [](const std::vector<std::string>& items, std::vector<uint8_t>& out)
    {
        out.clear();
        vint_code::WriteVint(static_cast<uint32_t>(items.size()), out);
        for (const auto& s : items)
        {
            vint_code::WriteVint(static_cast<uint32_t>(s.size()), out);
            out.insert(out.end(), s.begin(), s.end());
        }
    };
    serializeOne(ad_items, out_ad);
    serializeOne(pl_items, out_pl);
    serializeOne(pid_items, out_pid);
}

bool AdaptiveKnownFieldsDicts::Deserialize(const std::vector<uint8_t>& in, std::vector<std::string>& out_items)
{
    out_items.clear();
    size_t pos = 0;
    if (in.empty()) return false;
    const uint32_t count = vint_code::ReadVint(in.data(), in.size(), pos);
    out_items.reserve(count);
    for (uint32_t i = 0; i < count && pos < in.size(); ++i)
    {
        uint32_t len = vint_code::ReadVint(in.data(), in.size(), pos);
        if (pos + len > in.size()) return false;
        out_items.emplace_back(reinterpret_cast<const char*>(in.data() + pos), len);
        pos += len;
    }
    return true;
}

bool EncodeAdaptiveKnownFieldsRowV1(const bcf1_t* rec,
                                   uint32_t sample_count,
                                   const std::vector<field_desc>& fields,
                                   const AdaptiveKnownFieldsIndices& idx,
                                   AdaptiveKnownFieldsDicts& dicts,
                                   std::vector<uint8_t>& out_row)
{
    out_row.clear();
    out_row.reserve(128);
    out_row.push_back('F');
    out_row.push_back('M');
    out_row.push_back('D');
    out_row.push_back('1');

    const uint16_t allele = rec ? static_cast<uint16_t>(rec->n_allele) : 2;
    writeU16LE(allele, out_row);

    uint8_t mask = 0;
    if (idx.ad >= 0 && static_cast<size_t>(idx.ad) < fields.size() && fields[idx.ad].present) mask |= AdaptiveKnownFieldsRowV1::kMaskAD;
    if (idx.dp >= 0 && static_cast<size_t>(idx.dp) < fields.size() && fields[idx.dp].present) mask |= AdaptiveKnownFieldsRowV1::kMaskDP;
    if (idx.pl >= 0 && static_cast<size_t>(idx.pl) < fields.size() && fields[idx.pl].present) mask |= AdaptiveKnownFieldsRowV1::kMaskPL;
    if (idx.gq >= 0 && static_cast<size_t>(idx.gq) < fields.size() && fields[idx.gq].present) mask |= AdaptiveKnownFieldsRowV1::kMaskGQ;
    if (idx.pgt >= 0 && static_cast<size_t>(idx.pgt) < fields.size() && fields[idx.pgt].present) mask |= AdaptiveKnownFieldsRowV1::kMaskPGT;
    if (idx.pid >= 0 && static_cast<size_t>(idx.pid) < fields.size() && fields[idx.pid].present) mask |= AdaptiveKnownFieldsRowV1::kMaskPID;
    out_row.push_back(mask);

    const uint32_t adcnt = static_cast<uint32_t>(allele);
    const uint32_t plcnt = static_cast<uint32_t>(allele) * (static_cast<uint32_t>(allele) + 1u) / 2u;

    std::vector<uint32_t> ad_sum(sample_count, 0);
    std::vector<uint8_t> ad_has_sum(sample_count, 0);
    std::vector<int32_t> pl_pred_gq(sample_count, bcf_int32_missing);
    std::vector<uint8_t> pl_has_pred(sample_count, 0);

    // AD
    if (mask & AdaptiveKnownFieldsRowV1::kMaskAD)
    {
        const field_desc& f = fields[idx.ad];
        const int32_t* p = reinterpret_cast<const int32_t*>(f.data);
        const uint32_t total = f.data_size / 4u;
        const uint32_t vps = (sample_count ? (total / sample_count) : 0);

        std::vector<uint8_t> tips((sample_count + 3u) / 4u, 0);
        std::vector<uint32_t> miss_pos;
        std::vector<uint32_t> adids;
        std::vector<uint32_t> vals;
        vals.resize(adcnt);

        for (uint32_t s = 0; s < sample_count; ++s)
        {
            bool missing = (p == nullptr || vps < adcnt);
            if (!missing)
            {
                const int32_t* sp = p + static_cast<size_t>(s) * vps;
                if (sp[0] == bcf_int32_vector_end) missing = true;
                else if (sp[0] == bcf_int32_missing && (vps == 1 || sp[1] == bcf_int32_vector_end)) missing = true;
            }
            if (missing)
            {
                miss_pos.push_back(s);
                setTip2(tips, s, 0);
                continue;
            }

            const int32_t* sp = p + static_cast<size_t>(s) * vps;
            bool has_missing_elem = false;
            uint64_t sum = 0;
            for (uint32_t j = 0; j < adcnt; ++j)
            {
                uint32_t u = toU32OrMissing(sp[j]);
                vals[j] = u;
                if (u == kMissingU32) has_missing_elem = true;
                else sum += u;
            }

            uint8_t tip = 3; // 11 general
            if (!has_missing_elem)
            {
                if (sum == 0)
                {
                    tip = 0;
                }
                else
                {
                    bool only_first = true;
                    for (uint32_t j = 1; j < adcnt; ++j)
                    {
                        if (vals[j] != 0) { only_first = false; break; }
                    }
                    if (only_first && vals[0] == sum)
                    {
                        if (sum == 2) tip = 1;
                        else tip = 2;
                    }
                }
            }

            setTip2(tips, s, tip);
            if (tip == 2)
            {
                adids.push_back(static_cast<uint32_t>(sum));
                ad_sum[s] = static_cast<uint32_t>(sum);
                ad_has_sum[s] = 1;
            }
            else if (tip == 1)
            {
                ad_sum[s] = 2;
                ad_has_sum[s] = 1;
            }
            else if (tip == 0)
            {
                ad_sum[s] = 0;
                ad_has_sum[s] = 1;
            }
            else
            {
                uint32_t id = dicts.GetOrAddAD(vals, static_cast<uint32_t>(sum));
                adids.push_back(id);
                // For prediction, only use fully-defined arrays.
                if (!has_missing_elem)
                {
                    ad_sum[s] = static_cast<uint32_t>(sum);
                    ad_has_sum[s] = 1;
                }
            }
        }

        vint_code::WriteVint(static_cast<uint32_t>(tips.size()), out_row);
        out_row.insert(out_row.end(), tips.begin(), tips.end());
        encodePosListDeltas(miss_pos, out_row);
        vint_code::WriteVint(static_cast<uint32_t>(adids.size()), out_row);
        for (uint32_t v : adids) vint_code::WriteVint(v, out_row);
    }

    // DP
    if (mask & AdaptiveKnownFieldsRowV1::kMaskDP)
    {
        const field_desc& f = fields[idx.dp];
        const int32_t* p = reinterpret_cast<const int32_t*>(f.data);
        const uint32_t total = f.data_size / 4u;
        const uint32_t vps = (sample_count ? (total / sample_count) : 0);

        std::vector<uint32_t> ab_pos;
        std::vector<uint32_t> ab_val;
        std::vector<uint32_t> normal;
        ab_pos.reserve(128);
        ab_val.reserve(128);

        for (uint32_t s = 0; s < sample_count; ++s)
        {
            uint32_t dpv = kMissingU32;
            if (p != nullptr && vps >= 1)
            {
                const int32_t* sp = p + static_cast<size_t>(s) * vps;
                dpv = toU32OrMissing(sp[0]);
            }

            if (ad_has_sum[s])
            {
                const uint32_t pred = ad_sum[s];
                if (dpv != pred)
                {
                    ab_pos.push_back(s);
                    ab_val.push_back(dpv);
                }
            }
            else
            {
                normal.push_back(dpv);
            }
        }

        vint_code::WriteVint(static_cast<uint32_t>(ab_pos.size()), out_row);
        uint32_t prev = 0;
        for (size_t i = 0; i < ab_pos.size(); ++i)
        {
            vint_code::WriteVint(ab_pos[i] - prev, out_row);
            prev = ab_pos[i];
            vint_code::WriteVint(ab_val[i], out_row);
        }
        vint_code::WriteVint(static_cast<uint32_t>(normal.size()), out_row);
        for (uint32_t v : normal) vint_code::WriteVint(v, out_row);
    }

    // PL
    if (mask & AdaptiveKnownFieldsRowV1::kMaskPL)
    {
        const field_desc& f = fields[idx.pl];
        const int32_t* p = reinterpret_cast<const int32_t*>(f.data);
        const uint32_t total = f.data_size / 4u;
        const uint32_t vps = (sample_count ? (total / sample_count) : 0);

        std::vector<uint8_t> tips((sample_count + 3u) / 4u, 0);
        std::vector<uint32_t> miss_pos;
        std::vector<uint32_t> plids;
        std::vector<uint32_t> avec;
        std::vector<uint32_t> bvec;
        std::vector<uint32_t> vals;
        vals.resize(plcnt);

        for (uint32_t s = 0; s < sample_count; ++s)
        {
            bool missing = (p == nullptr || vps < plcnt);
            if (!missing)
            {
                const int32_t* sp = p + static_cast<size_t>(s) * vps;
                if (sp[0] == bcf_int32_vector_end) missing = true;
                else if (sp[0] == bcf_int32_missing && (vps == 1 || sp[1] == bcf_int32_vector_end)) missing = true;
            }
            if (missing)
            {
                miss_pos.push_back(s);
                setTip2(tips, s, 0);
                continue;
            }

            const int32_t* sp = p + static_cast<size_t>(s) * vps;
            bool has_missing_elem = false;
            for (uint32_t j = 0; j < plcnt; ++j)
            {
                uint32_t u = toU32OrMissing(sp[j]);
                vals[j] = u;
                if (u == kMissingU32) has_missing_elem = true;
            }

            uint8_t tip = 3;
            uint32_t a = 0, b = 0;
            uint32_t max_val = 0, sec_min = 0;
            if (!has_missing_elem)
            {
                uint8_t t = checkAbPattern(vals.data(), plcnt, a, b);
                if (t == 0) tip = 0;
                else if (t == 1) tip = 1;
                else if (t == 2) tip = 2;
                else tip = 3;

                if (tip == 3)
                {
                    getMaxAndSecondMin(vals.data(), plcnt, max_val, sec_min);
                }
                else if (tip == 1 || tip == 2)
                {
                    sec_min = a;
                }
                else
                {
                    sec_min = kMissingU32;
                }
            }

            setTip2(tips, s, tip);
            if (tip == 1)
            {
                avec.push_back(a);
                bvec.push_back(b);
                pl_pred_gq[s] = static_cast<int32_t>(a);
                pl_has_pred[s] = 1;
            }
            else if (tip == 2)
            {
                avec.push_back(a);
                pl_pred_gq[s] = static_cast<int32_t>(a);
                pl_has_pred[s] = 1;
            }
            else if (tip == 3)
            {
                uint32_t id = dicts.GetOrAddPL(vals, max_val);
                plids.push_back(id);
                if (sec_min != kMissingU32)
                {
                    pl_pred_gq[s] = static_cast<int32_t>(sec_min);
                    pl_has_pred[s] = 1;
                }
            }
            else
            {
                // type0: predicted "." for GQ
                pl_has_pred[s] = 1;
                pl_pred_gq[s] = bcf_int32_missing;
            }
        }

        vint_code::WriteVint(static_cast<uint32_t>(tips.size()), out_row);
        out_row.insert(out_row.end(), tips.begin(), tips.end());
        encodePosListDeltas(miss_pos, out_row);
        vint_code::WriteVint(static_cast<uint32_t>(plids.size()), out_row);
        for (uint32_t v : plids) vint_code::WriteVint(v, out_row);
        vint_code::WriteVint(static_cast<uint32_t>(avec.size()), out_row);
        for (uint32_t v : avec) vint_code::WriteVint(v, out_row);
        vint_code::WriteVint(static_cast<uint32_t>(bvec.size()), out_row);
        for (uint32_t v : bvec) vint_code::WriteVint(v, out_row);
    }

    // GQ
    if (mask & AdaptiveKnownFieldsRowV1::kMaskGQ)
    {
        const field_desc& f = fields[idx.gq];
        const int32_t* p = reinterpret_cast<const int32_t*>(f.data);
        const uint32_t total = f.data_size / 4u;
        const uint32_t vps = (sample_count ? (total / sample_count) : 0);

        std::vector<uint32_t> ab_pos;
        std::vector<uint32_t> ab_val;
        std::vector<uint32_t> normal;

        for (uint32_t s = 0; s < sample_count; ++s)
        {
            uint32_t gqv = kMissingU32;
            if (p != nullptr && vps >= 1)
            {
                const int32_t* sp = p + static_cast<size_t>(s) * vps;
                gqv = toU32OrMissing(sp[0]);
            }

            if (pl_has_pred[s])
            {
                const int32_t pred = pl_pred_gq[s];
                const uint32_t pred_u = (pred == bcf_int32_missing) ? kMissingU32 : static_cast<uint32_t>(pred);
                if (gqv != pred_u)
                {
                    ab_pos.push_back(s);
                    ab_val.push_back(gqv);
                }
            }
            else
            {
                normal.push_back(gqv);
            }
        }

        vint_code::WriteVint(static_cast<uint32_t>(ab_pos.size()), out_row);
        uint32_t prev = 0;
        for (size_t i = 0; i < ab_pos.size(); ++i)
        {
            vint_code::WriteVint(ab_pos[i] - prev, out_row);
            prev = ab_pos[i];
            vint_code::WriteVint(ab_val[i], out_row);
        }
        vint_code::WriteVint(static_cast<uint32_t>(normal.size()), out_row);
        for (uint32_t v : normal) vint_code::WriteVint(v, out_row);
    }

    // PGT
    if (mask & AdaptiveKnownFieldsRowV1::kMaskPGT)
    {
        std::vector<std::string> vals;
        if (!parseFormatStrings(fields[idx.pgt], sample_count, vals))
        {
            encodePosListDeltas(std::vector<uint32_t>(), out_row);
        }
        else
        {
            std::vector<uint32_t> pos;
            std::vector<std::string> outv;
            for (uint32_t s = 0; s < sample_count; ++s)
            {
                if (vals[s] == "." || vals[s].empty()) continue;
                pos.push_back(s);
                outv.push_back(vals[s]);
            }
            encodePosListDeltas(pos, out_row);
            for (const auto& str : outv)
            {
                vint_code::WriteVint(static_cast<uint32_t>(str.size()), out_row);
                out_row.insert(out_row.end(), str.begin(), str.end());
            }
        }
    }

    // PID
    if (mask & AdaptiveKnownFieldsRowV1::kMaskPID)
    {
        std::vector<std::string> vals;
        if (!parseFormatStrings(fields[idx.pid], sample_count, vals))
        {
            encodePosListDeltas(std::vector<uint32_t>(), out_row);
            vint_code::WriteVint(0, out_row);
        }
        else
        {
            std::vector<uint32_t> pos;
            std::vector<uint32_t> ids;
            for (uint32_t s = 0; s < sample_count; ++s)
            {
                if (vals[s] == "." || vals[s].empty()) continue;
                pos.push_back(s);
                ids.push_back(dicts.GetOrAddPID(vals[s].data(), vals[s].size()));
            }
            encodePosListDeltas(pos, out_row);
            vint_code::WriteVint(static_cast<uint32_t>(ids.size()), out_row);
            for (uint32_t id : ids) vint_code::WriteVint(id, out_row);
        }
    }

    return true;
}

bool AdaptiveKnownFieldsRowV1::Parse(const uint8_t* data, size_t len)
{
    *this = AdaptiveKnownFieldsRowV1();
    if (data == nullptr || len < 7) return false;
    if (std::memcmp(data, kMagic, 4) != 0) return false;

    size_t pos = 4;
    if (!readU16LE(data, len, pos, allele_count)) return false;
    if (pos >= len) return false;
    presence_mask = data[pos++];

    auto readTips = [&](std::vector<uint8_t>& tips) -> bool
    {
        uint32_t sz = vint_code::ReadVint(data, len, pos);
        if (pos + sz > len) return false;
        tips.assign(data + pos, data + pos + sz);
        pos += sz;
        return true;
    };

    auto readVecU32 = [&](std::vector<uint32_t>& v) -> bool
    {
        uint32_t n = vint_code::ReadVint(data, len, pos);
        v.resize(n);
        for (uint32_t i = 0; i < n; ++i)
            v[i] = vint_code::ReadVint(data, len, pos);
        return true;
    };

    auto readAbPairs = [&](std::vector<uint32_t>& p, std::vector<uint32_t>& v) -> bool
    {
        uint32_t n = vint_code::ReadVint(data, len, pos);
        p.resize(n);
        v.resize(n);
        uint32_t cur = 0;
        for (uint32_t i = 0; i < n; ++i)
        {
            cur += vint_code::ReadVint(data, len, pos);
            p[i] = cur;
            v[i] = vint_code::ReadVint(data, len, pos);
        }
        return true;
    };

    if (Has(kMaskAD))
    {
        if (!readTips(ad_tips)) return false;
        if (!decodePosListDeltas(data, len, pos, ad_miss_pos)) return false;
        if (!readVecU32(ad_ids)) return false;
    }
    if (Has(kMaskDP))
    {
        if (!readAbPairs(dp_ab_pos, dp_ab_val)) return false;
        if (!readVecU32(dp_normal)) return false;
    }
    if (Has(kMaskPL))
    {
        if (!readTips(pl_tips)) return false;
        if (!decodePosListDeltas(data, len, pos, pl_miss_pos)) return false;
        if (!readVecU32(pl_ids)) return false;
        if (!readVecU32(pl_a)) return false;
        if (!readVecU32(pl_b)) return false;
    }
    if (Has(kMaskGQ))
    {
        if (!readAbPairs(gq_ab_pos, gq_ab_val)) return false;
        if (!readVecU32(gq_normal)) return false;
    }
    if (Has(kMaskPGT))
    {
        if (!decodePosListDeltas(data, len, pos, pgt_pos)) return false;
        pgt_vals.resize(pgt_pos.size());
        for (size_t i = 0; i < pgt_pos.size(); ++i)
        {
            uint32_t slen = vint_code::ReadVint(data, len, pos);
            if (pos + slen > len) return false;
            pgt_vals[i].assign(reinterpret_cast<const char*>(data + pos), slen);
            pos += slen;
        }
    }
    if (Has(kMaskPID))
    {
        if (!decodePosListDeltas(data, len, pos, pid_pos)) return false;
        uint32_t n = vint_code::ReadVint(data, len, pos);
        pid_ids.resize(n);
        for (uint32_t i = 0; i < n; ++i)
            pid_ids[i] = vint_code::ReadVint(data, len, pos);
    }

    return true;
}

bool AdaptiveKnownFieldsRowV1::DecodeAD(uint32_t sample_count,
                                       const std::vector<std::string>& ad_dict,
                                       std::vector<int32_t>& out,
                                       std::vector<uint32_t>* out_sum,
                                       std::vector<uint8_t>* out_has_sum) const
{
    if (!Has(kMaskAD)) return false;
    const uint32_t adcnt = ADCnt();
    out.assign(static_cast<size_t>(sample_count) * adcnt, bcf_int32_vector_end);
    if (out_sum) out_sum->assign(sample_count, 0);
    if (out_has_sum) out_has_sum->assign(sample_count, 0);

    size_t miss_i = 0;
    size_t id_i = 0;
    std::vector<uint32_t> dict_vals;
    uint32_t width = 0;

    for (uint32_t s = 0; s < sample_count; ++s)
    {
        const size_t off = static_cast<size_t>(s) * adcnt;
        const bool is_miss = (miss_i < ad_miss_pos.size() && ad_miss_pos[miss_i] == s);
        if (is_miss)
        {
            ++miss_i;
            out[off] = bcf_int32_missing;
            continue;
        }

        const uint8_t tip = getTip2(ad_tips, s);
        if (tip == 0)
        {
            for (uint32_t j = 0; j < adcnt; ++j) out[off + j] = 0;
            if (out_sum) (*out_sum)[s] = 0;
            if (out_has_sum) (*out_has_sum)[s] = 1;
        }
        else if (tip == 1)
        {
            out[off] = 2;
            for (uint32_t j = 1; j < adcnt; ++j) out[off + j] = 0;
            if (out_sum) (*out_sum)[s] = 2;
            if (out_has_sum) (*out_has_sum)[s] = 1;
        }
        else if (tip == 2)
        {
            if (id_i >= ad_ids.size()) { out[off] = bcf_int32_missing; continue; }
            const uint32_t sum = ad_ids[id_i++];
            out[off] = static_cast<int32_t>(sum);
            for (uint32_t j = 1; j < adcnt; ++j) out[off + j] = 0;
            if (out_sum) (*out_sum)[s] = sum;
            if (out_has_sum) (*out_has_sum)[s] = 1;
        }
        else
        {
            if (id_i >= ad_ids.size()) { out[off] = bcf_int32_missing; continue; }
            const uint32_t dict_id = ad_ids[id_i++];
            if (dict_id >= ad_dict.size()) { out[off] = bcf_int32_missing; continue; }
            if (!unpackFixedWidthArrayItem(ad_dict[dict_id], adcnt, dict_vals, width))
            {
                out[off] = bcf_int32_missing;
                continue;
            }
            uint64_t sum = 0;
            bool has_missing = false;
            for (uint32_t j = 0; j < adcnt; ++j)
            {
                uint32_t u = dict_vals[j];
                if (u == kMissingU32)
                {
                    out[off + j] = bcf_int32_missing;
                    has_missing = true;
                }
                else
                {
                    out[off + j] = static_cast<int32_t>(u);
                    sum += u;
                }
            }
            if (out_has_sum && !has_missing) (*out_has_sum)[s] = 1;
            if (out_sum && !has_missing) (*out_sum)[s] = static_cast<uint32_t>(sum);
        }
    }

    return true;
}

bool AdaptiveKnownFieldsRowV1::DecodeDP(uint32_t sample_count,
                                       const std::vector<uint32_t>& ad_sum,
                                       const std::vector<uint8_t>& ad_has_sum,
                                       std::vector<int32_t>& out) const
{
    if (!Has(kMaskDP)) return false;
    out.assign(sample_count, bcf_int32_missing);

    size_t ab_i = 0;
    size_t norm_i = 0;

    for (uint32_t s = 0; s < sample_count; ++s)
    {
        uint32_t v = kMissingU32;
        if (s < ad_has_sum.size() && ad_has_sum[s])
        {
            v = (s < ad_sum.size()) ? ad_sum[s] : kMissingU32;
        }
        else
        {
            if (norm_i < dp_normal.size()) v = dp_normal[norm_i++];
        }

        if (ab_i < dp_ab_pos.size() && dp_ab_pos[ab_i] == s)
        {
            v = dp_ab_val[ab_i];
            ++ab_i;
        }
        out[s] = fromU32OrMissing(v);
    }
    return true;
}

bool AdaptiveKnownFieldsRowV1::DecodePL(uint32_t sample_count,
                                       const std::vector<std::string>& pl_dict,
                                       std::vector<int32_t>& out,
                                       std::vector<int32_t>* out_pred_gq,
                                       std::vector<uint8_t>* out_has_pred) const
{
    if (!Has(kMaskPL)) return false;
    const uint32_t plcnt = PLCnt();
    out.assign(static_cast<size_t>(sample_count) * plcnt, bcf_int32_vector_end);
    if (out_pred_gq) out_pred_gq->assign(sample_count, bcf_int32_missing);
    if (out_has_pred) out_has_pred->assign(sample_count, 0);

    size_t miss_i = 0;
    size_t id_i = 0;
    size_t a_i = 0;
    size_t b_i = 0;

    std::vector<uint32_t> dict_vals;
    uint32_t width = 0;

    for (uint32_t s = 0; s < sample_count; ++s)
    {
        const size_t off = static_cast<size_t>(s) * plcnt;
        const bool is_miss = (miss_i < pl_miss_pos.size() && pl_miss_pos[miss_i] == s);
        if (is_miss)
        {
            ++miss_i;
            out[off] = bcf_int32_missing;
            continue;
        }

        const uint8_t tip = getTip2(pl_tips, s);
        if (tip == 0)
        {
            for (uint32_t j = 0; j < plcnt; ++j) out[off + j] = 0;
            if (out_has_pred) (*out_has_pred)[s] = 1;
            if (out_pred_gq) (*out_pred_gq)[s] = bcf_int32_missing;
        }
        else if (tip == 1)
        {
            if (a_i >= pl_a.size() || b_i >= pl_b.size()) { out[off] = bcf_int32_missing; continue; }
            uint32_t a = pl_a[a_i++];
            uint32_t b = pl_b[b_i++];
            patternToPLArray(a, b, allele_count, out, off);
            if (out_has_pred) (*out_has_pred)[s] = 1;
            if (out_pred_gq) (*out_pred_gq)[s] = static_cast<int32_t>(a);
        }
        else if (tip == 2)
        {
            if (a_i >= pl_a.size()) { out[off] = bcf_int32_missing; continue; }
            uint32_t a = pl_a[a_i++];
            uint32_t b = a * 15u;
            patternToPLArray(a, b, allele_count, out, off);
            if (out_has_pred) (*out_has_pred)[s] = 1;
            if (out_pred_gq) (*out_pred_gq)[s] = static_cast<int32_t>(a);
        }
        else
        {
            if (id_i >= pl_ids.size()) { out[off] = bcf_int32_missing; continue; }
            const uint32_t dict_id = pl_ids[id_i++];
            if (dict_id >= pl_dict.size()) { out[off] = bcf_int32_missing; continue; }
            if (!unpackFixedWidthArrayItem(pl_dict[dict_id], plcnt, dict_vals, width))
            {
                out[off] = bcf_int32_missing;
                continue;
            }
            std::vector<uint32_t> tmp = dict_vals;
            bool has_missing = false;
            for (uint32_t j = 0; j < plcnt; ++j)
            {
                uint32_t u = tmp[j];
                if (u == kMissingU32) { out[off + j] = bcf_int32_missing; has_missing = true; }
                else out[off + j] = static_cast<int32_t>(u);
            }
            if (!has_missing && out_has_pred && out_pred_gq)
            {
                uint32_t maxv = 0, secmin = 0;
                getMaxAndSecondMin(tmp.data(), plcnt, maxv, secmin);
                (*out_has_pred)[s] = 1;
                (*out_pred_gq)[s] = static_cast<int32_t>(secmin);
            }
        }
    }
    return true;
}

bool AdaptiveKnownFieldsRowV1::DecodeGQ(uint32_t sample_count,
                                       const std::vector<int32_t>& pl_pred_gq,
                                       const std::vector<uint8_t>& pl_has_pred,
                                       std::vector<int32_t>& out) const
{
    if (!Has(kMaskGQ)) return false;
    out.assign(sample_count, bcf_int32_missing);

    size_t ab_i = 0;
    size_t norm_i = 0;
    for (uint32_t s = 0; s < sample_count; ++s)
    {
        uint32_t v = kMissingU32;
        if (s < pl_has_pred.size() && pl_has_pred[s])
        {
            int32_t pred = (s < pl_pred_gq.size()) ? pl_pred_gq[s] : bcf_int32_missing;
            v = (pred == bcf_int32_missing) ? kMissingU32 : static_cast<uint32_t>(pred);
        }
        else
        {
            if (norm_i < gq_normal.size()) v = gq_normal[norm_i++];
        }

        if (ab_i < gq_ab_pos.size() && gq_ab_pos[ab_i] == s)
        {
            v = gq_ab_val[ab_i];
            ++ab_i;
        }
        out[s] = fromU32OrMissing(v);
    }
    return true;
}

bool AdaptiveKnownFieldsRowV1::DecodePGT(uint32_t sample_count, std::vector<std::string>& out) const
{
    if (!Has(kMaskPGT)) return false;
    out.assign(sample_count, ".");
    for (size_t i = 0; i < pgt_pos.size() && i < pgt_vals.size(); ++i)
    {
        if (pgt_pos[i] < sample_count) out[pgt_pos[i]] = pgt_vals[i];
    }
    return true;
}

bool AdaptiveKnownFieldsRowV1::DecodePID(uint32_t sample_count,
                                        const std::vector<std::string>& pid_dict,
                                        std::vector<std::string>& out) const
{
    if (!Has(kMaskPID)) return false;
    out.assign(sample_count, ".");
    for (size_t i = 0; i < pid_pos.size() && i < pid_ids.size(); ++i)
    {
        uint32_t p = pid_pos[i];
        uint32_t id = pid_ids[i];
        if (p >= sample_count || id >= pid_dict.size()) continue;
        const std::string& item = pid_dict[id];
        if (item.size() < 2) continue;
        uint16_t l = static_cast<uint8_t>(item[0]) | (static_cast<uint16_t>(static_cast<uint8_t>(item[1])) << 8);
        if (item.size() < static_cast<size_t>(l) + 2) continue;
        out[p].assign(item.data() + 2, l);
    }
    return true;
}

} // namespace gsc
