/**
 * @file gvcf_block.cpp
 * @brief Implementation of gVCF block structures
 */

#include "gvcf_block.h"
#include "gvcf_encoding.h"
#include <cstring>

namespace gvcf {

// ============================================================================
// CompressedGVCFBlock Implementation
// ============================================================================

size_t CompressedGVCFBlock::TotalCompressedSize() const {
    size_t total = 0;

    total += chrom.data.size();
    total += pos.data.size();
    total += id.data.size();
    total += ref.data.size();
    total += alt.data.size();
    total += qual.data.size();
    total += filter.data.size();
    total += info_end.data.size();

    total += gt_mask.data.size();
    total += gt_patches.data.size();
    total += gt_phase.data.size();

    total += dp.data.size();
    total += gq.data.size();
    total += min_dp.data.size();
    total += dp_min_dp_diff.data.size();
    total += pl.data.size();
    total += ad.data.size();

    for (const auto& kv : unknown_info) {
        const auto& name = kv.first;
        const auto& field = kv.second;
        total += name.size() + field.data.size();
    }

    for (const auto& kv : unknown_format) {
        const auto& name = kv.first;
        const auto& field = kv.second;
        total += name.size() + field.data.size();
    }

    return total;
}

namespace {

// Helper to serialize a CompressedField
void SerializeField(const CompressedField& field, std::vector<uint8_t>& buffer) {
    // Method
    buffer.push_back(static_cast<uint8_t>(field.method));

    // Original count
    VarIntUtil::WriteVarUint(field.original_count, buffer);

    // Data size and data
    VarIntUtil::WriteVarUint(field.data.size(), buffer);
    buffer.insert(buffer.end(), field.data.begin(), field.data.end());
}

// Helper to deserialize a CompressedField
bool DeserializeField(const uint8_t* buffer, size_t size, size_t& pos,
                     CompressedField& field) {
    if (pos >= size) return false;

    field.Clear();

    // Method
    field.method = static_cast<FieldCompressionMethod>(buffer[pos++]);

    // Original count
    field.original_count = static_cast<uint32_t>(
        VarIntUtil::ReadVarUint(buffer, size, pos));

    // Data size and data
    uint64_t data_size = VarIntUtil::ReadVarUint(buffer, size, pos);
    if (pos + data_size > size) return false;

    field.data.assign(buffer + pos, buffer + pos + data_size);
    pos += data_size;

    return true;
}

} // anonymous namespace

bool CompressedGVCFBlock::Serialize(std::vector<uint8_t>& buffer) const {
    buffer.clear();

    // Header
    // Magic number (4 bytes)
    buffer.push_back(0x47); // 'G'
    buffer.push_back(0x56); // 'V'
    buffer.push_back(0x43); // 'C'
    buffer.push_back(0x46); // 'F'

    // Version (1 byte)
    buffer.push_back(1);

    // Metadata
    VarIntUtil::WriteVarUint(variant_count, buffer);
    VarIntUtil::WriteVarUint(sample_count, buffer);

    // Flags
    uint8_t flags = 0;
    if (has_end_field) flags |= 0x01;
    if (has_min_dp) flags |= 0x02;
    buffer.push_back(flags);

    // Serialize position fields
    SerializeField(chrom, buffer);
    SerializeField(pos, buffer);
    SerializeField(id, buffer);

    // Serialize sequence fields
    SerializeField(ref, buffer);
    SerializeField(alt, buffer);

    // Serialize quality fields
    SerializeField(qual, buffer);
    SerializeField(filter, buffer);

    // Serialize INFO/END
    if (has_end_field) {
        SerializeField(info_end, buffer);
    }

    // Serialize sample fields
    SerializeField(gt_mask, buffer);
    SerializeField(gt_patches, buffer);
    SerializeField(gt_phase, buffer);

    SerializeField(dp, buffer);
    SerializeField(gq, buffer);

    if (has_min_dp) {
        SerializeField(min_dp, buffer);
        SerializeField(dp_min_dp_diff, buffer);
    }

    SerializeField(pl, buffer);
    SerializeField(ad, buffer);

    // Serialize unknown INFO fields
    VarIntUtil::WriteVarUint(unknown_info.size(), buffer);
    for (const auto& kv : unknown_info) {
        const auto& name = kv.first;
        const auto& field = kv.second;
        VarIntUtil::WriteVarUint(name.size(), buffer);
        buffer.insert(buffer.end(), name.begin(), name.end());
        SerializeField(field, buffer);
    }

    // Serialize unknown FORMAT fields
    VarIntUtil::WriteVarUint(unknown_format.size(), buffer);
    for (const auto& kv : unknown_format) {
        const auto& name = kv.first;
        const auto& field = kv.second;
        VarIntUtil::WriteVarUint(name.size(), buffer);
        buffer.insert(buffer.end(), name.begin(), name.end());
        SerializeField(field, buffer);
    }

    return true;
}

bool CompressedGVCFBlock::Deserialize(const std::vector<uint8_t>& buffer) {
    return Deserialize(buffer.data(), buffer.size());
}

bool CompressedGVCFBlock::Deserialize(const uint8_t* buffer, size_t size) {
    Clear();

    size_t pos = 0;

    // Check magic
    if (size < 5) return false;
    if (buffer[0] != 0x47 || buffer[1] != 0x56 ||
        buffer[2] != 0x43 || buffer[3] != 0x46) {
        return false;
    }
    pos = 4;

    // Version
    uint8_t version = buffer[pos++];
    if (version != 1) return false;

    // Metadata
    variant_count = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
    sample_count = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));

    // Flags
    uint8_t flags = buffer[pos++];
    has_end_field = (flags & 0x01) != 0;
    has_min_dp = (flags & 0x02) != 0;

    // Deserialize position fields
    if (!DeserializeField(buffer, size, pos, chrom)) return false;
    if (!DeserializeField(buffer, size, pos, this->pos)) return false;
    if (!DeserializeField(buffer, size, pos, id)) return false;

    // Deserialize sequence fields
    if (!DeserializeField(buffer, size, pos, ref)) return false;
    if (!DeserializeField(buffer, size, pos, alt)) return false;

    // Deserialize quality fields
    if (!DeserializeField(buffer, size, pos, qual)) return false;
    if (!DeserializeField(buffer, size, pos, filter)) return false;

    // Deserialize INFO/END
    if (has_end_field) {
        if (!DeserializeField(buffer, size, pos, info_end)) return false;
    }

    // Deserialize sample fields
    if (!DeserializeField(buffer, size, pos, gt_mask)) return false;
    if (!DeserializeField(buffer, size, pos, gt_patches)) return false;
    if (!DeserializeField(buffer, size, pos, gt_phase)) return false;

    if (!DeserializeField(buffer, size, pos, dp)) return false;
    if (!DeserializeField(buffer, size, pos, gq)) return false;

    if (has_min_dp) {
        if (!DeserializeField(buffer, size, pos, min_dp)) return false;
        if (!DeserializeField(buffer, size, pos, dp_min_dp_diff)) return false;
    }

    if (!DeserializeField(buffer, size, pos, pl)) return false;
    if (!DeserializeField(buffer, size, pos, ad)) return false;

    // Deserialize unknown INFO fields
    uint64_t info_count = VarIntUtil::ReadVarUint(buffer, size, pos);
    for (uint64_t i = 0; i < info_count; ++i) {
        uint64_t name_len = VarIntUtil::ReadVarUint(buffer, size, pos);
        std::string name(reinterpret_cast<const char*>(buffer + pos), name_len);
        pos += name_len;

        CompressedField field;
        if (!DeserializeField(buffer, size, pos, field)) return false;
        unknown_info[name] = std::move(field);
    }

    // Deserialize unknown FORMAT fields
    uint64_t format_count = VarIntUtil::ReadVarUint(buffer, size, pos);
    for (uint64_t i = 0; i < format_count; ++i) {
        uint64_t name_len = VarIntUtil::ReadVarUint(buffer, size, pos);
        std::string name(reinterpret_cast<const char*>(buffer + pos), name_len);
        pos += name_len;

        CompressedField field;
        if (!DeserializeField(buffer, size, pos, field)) return false;
        unknown_format[name] = std::move(field);
    }

    return true;
}

} // namespace gvcf
