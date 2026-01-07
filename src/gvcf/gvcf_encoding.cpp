/**
 * @file gvcf_encoding.cpp
 * @brief Implementation of gVCF encoding algorithms
 */

#include "gvcf_encoding.h"
#include <cstring>
#include <stdexcept>

namespace gvcf {

// ============================================================================
// VarIntUtil Implementation
// ============================================================================

size_t VarIntUtil::WriteVarUint(uint64_t value, std::vector<uint8_t>& buffer) {
    size_t bytes_written = 0;
    do {
        uint8_t byte = value & 0x7F;
        value >>= 7;
        if (value != 0) {
            byte |= 0x80;
        }
        buffer.push_back(byte);
        bytes_written++;
    } while (value != 0);
    return bytes_written;
}

size_t VarIntUtil::WriteVarUint(uint64_t value, uint8_t* buffer) {
    size_t bytes_written = 0;
    do {
        uint8_t byte = value & 0x7F;
        value >>= 7;
        if (value != 0) {
            byte |= 0x80;
        }
        buffer[bytes_written++] = byte;
    } while (value != 0);
    return bytes_written;
}

uint64_t VarIntUtil::ReadVarUint(const std::vector<uint8_t>& buffer, size_t& pos) {
    uint64_t result = 0;
    int shift = 0;
    while (pos < buffer.size()) {
        uint8_t byte = buffer[pos++];
        result |= static_cast<uint64_t>(byte & 0x7F) << shift;
        if ((byte & 0x80) == 0) {
            break;
        }
        shift += 7;
    }
    return result;
}

uint64_t VarIntUtil::ReadVarUint(const uint8_t* buffer, size_t size, size_t& pos) {
    uint64_t result = 0;
    int shift = 0;
    while (pos < size) {
        uint8_t byte = buffer[pos++];
        result |= static_cast<uint64_t>(byte & 0x7F) << shift;
        if ((byte & 0x80) == 0) {
            break;
        }
        shift += 7;
    }
    return result;
}

uint64_t VarIntUtil::ZigzagEncode(int64_t value) {
    return (static_cast<uint64_t>(value) << 1) ^ (value >> 63);
}

int64_t VarIntUtil::ZigzagDecode(uint64_t value) {
    return static_cast<int64_t>((value >> 1) ^ -(value & 1));
}

size_t VarIntUtil::WriteVarInt(int64_t value, std::vector<uint8_t>& buffer) {
    return WriteVarUint(ZigzagEncode(value), buffer);
}

int64_t VarIntUtil::ReadVarInt(const std::vector<uint8_t>& buffer, size_t& pos) {
    return ZigzagDecode(ReadVarUint(buffer, pos));
}

int64_t VarIntUtil::ReadVarInt(const uint8_t* buffer, size_t size, size_t& pos) {
    return ZigzagDecode(ReadVarUint(buffer, size, pos));
}

// ============================================================================
// BitmapUtil Implementation
// ============================================================================

std::vector<uint8_t> BitmapUtil::FromBools(const std::vector<bool>& data) {
    size_t byte_count = (data.size() + 7) / 8;
    std::vector<uint8_t> bitmap(byte_count, 0);

    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i]) {
            bitmap[i / 8] |= (1 << (i % 8));
        }
    }
    return bitmap;
}

std::vector<bool> BitmapUtil::ToBools(const std::vector<uint8_t>& bitmap, uint32_t len) {
    std::vector<bool> result(len, false);
    for (uint32_t i = 0; i < len; ++i) {
        result[i] = (bitmap[i / 8] & (1 << (i % 8))) != 0;
    }
    return result;
}

bool BitmapUtil::GetBit(const std::vector<uint8_t>& bitmap, uint32_t pos) {
    return (bitmap[pos / 8] & (1 << (pos % 8))) != 0;
}

void BitmapUtil::SetBit(std::vector<uint8_t>& bitmap, uint32_t pos, bool value) {
    if (value) {
        bitmap[pos / 8] |= (1 << (pos % 8));
    } else {
        bitmap[pos / 8] &= ~(1 << (pos % 8));
    }
}

uint32_t BitmapUtil::PopCount(const std::vector<uint8_t>& bitmap) {
    uint32_t count = 0;
    for (uint8_t byte : bitmap) {
        count += __builtin_popcount(byte);
    }
    return count;
}

uint32_t BitmapUtil::PopCount(const std::vector<uint8_t>& bitmap, uint32_t len) {
    uint32_t count = 0;
    for (uint32_t i = 0; i < len; ++i) {
        if (GetBit(bitmap, i)) {
            count++;
        }
    }
    return count;
}

bool BitmapUtil::CompressRLE(const std::vector<uint8_t>& bitmap, RLEByteResult& result) {
    return RLEByteEncoder::Compress(bitmap, result);
}

// ============================================================================
// RLEEncoder Implementation (String)
// ============================================================================

bool RLEEncoder::Compress(const std::vector<std::string>& data, RLEResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    std::string current = data[0];
    uint32_t count = 1;

    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] == current) {
            count++;
        } else {
            result.runs.push_back({current, count});
            current = data[i];
            count = 1;
        }
    }
    result.runs.push_back({current, count});

    return true;
}

bool RLEEncoder::Decompress(const RLEResult& result, std::vector<std::string>& data) {
    data.clear();
    data.reserve(result.original_len);

    for (const auto& run : result.runs) {
        for (uint32_t i = 0; i < run.count; ++i) {
            data.push_back(run.value);
        }
    }

    return data.size() == result.original_len;
}

bool RLEEncoder::Serialize(const RLEResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    // Write original length
    VarIntUtil::WriteVarUint(result.original_len, buffer);

    // Write number of runs
    VarIntUtil::WriteVarUint(result.runs.size(), buffer);

    // Write each run
    for (const auto& run : result.runs) {
        // Write string length and data
        VarIntUtil::WriteVarUint(run.value.size(), buffer);
        buffer.insert(buffer.end(), run.value.begin(), run.value.end());

        // Write count
        VarIntUtil::WriteVarUint(run.count, buffer);
    }

    return true;
}

bool RLEEncoder::Deserialize(const std::vector<uint8_t>& buffer, RLEResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool RLEEncoder::Deserialize(const uint8_t* buffer, size_t size, RLEResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
    uint64_t num_runs = VarIntUtil::ReadVarUint(buffer, size, pos);

    result.runs.reserve(num_runs);
    for (uint64_t i = 0; i < num_runs; ++i) {
        RLERun run;
        uint64_t str_len = VarIntUtil::ReadVarUint(buffer, size, pos);
        run.value.assign(reinterpret_cast<const char*>(buffer + pos), str_len);
        pos += str_len;
        run.count = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
        result.runs.push_back(std::move(run));
    }

    return true;
}

// ============================================================================
// RLEIntEncoder Implementation
// ============================================================================

bool RLEIntEncoder::Compress(const std::vector<int32_t>& data, RLEIntResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    int32_t current = data[0];
    uint32_t count = 1;

    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] == current) {
            count++;
        } else {
            result.runs.push_back({current, count});
            current = data[i];
            count = 1;
        }
    }
    result.runs.push_back({current, count});

    return true;
}

bool RLEIntEncoder::Decompress(const RLEIntResult& result, std::vector<int32_t>& data) {
    data.clear();
    data.reserve(result.original_len);

    for (const auto& run : result.runs) {
        for (uint32_t i = 0; i < run.count; ++i) {
            data.push_back(run.value);
        }
    }

    return data.size() == result.original_len;
}

bool RLEIntEncoder::Serialize(const RLEIntResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    VarIntUtil::WriteVarUint(result.original_len, buffer);
    VarIntUtil::WriteVarUint(result.runs.size(), buffer);

    for (const auto& run : result.runs) {
        VarIntUtil::WriteVarInt(run.value, buffer);
        VarIntUtil::WriteVarUint(run.count, buffer);
    }

    return true;
}

bool RLEIntEncoder::Deserialize(const std::vector<uint8_t>& buffer, RLEIntResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool RLEIntEncoder::Deserialize(const uint8_t* buffer, size_t size, RLEIntResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
    uint64_t num_runs = VarIntUtil::ReadVarUint(buffer, size, pos);

    result.runs.reserve(num_runs);
    for (uint64_t i = 0; i < num_runs; ++i) {
        RLEIntRun run;
        run.value = static_cast<int32_t>(VarIntUtil::ReadVarInt(buffer, size, pos));
        run.count = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
        result.runs.push_back(run);
    }

    return true;
}

// ============================================================================
// RLEByteEncoder Implementation
// ============================================================================

bool RLEByteEncoder::Compress(const std::vector<uint8_t>& data, RLEByteResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    uint8_t current = data[0];
    uint32_t count = 1;

    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] == current) {
            count++;
        } else {
            result.runs.push_back({current, count});
            current = data[i];
            count = 1;
        }
    }
    result.runs.push_back({current, count});

    return true;
}

bool RLEByteEncoder::Decompress(const RLEByteResult& result, std::vector<uint8_t>& data) {
    data.clear();
    data.reserve(result.original_len);

    for (const auto& run : result.runs) {
        for (uint32_t i = 0; i < run.count; ++i) {
            data.push_back(run.value);
        }
    }

    return data.size() == result.original_len;
}

bool RLEByteEncoder::Serialize(const RLEByteResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    VarIntUtil::WriteVarUint(result.original_len, buffer);
    VarIntUtil::WriteVarUint(result.runs.size(), buffer);

    for (const auto& run : result.runs) {
        buffer.push_back(run.value);
        VarIntUtil::WriteVarUint(run.count, buffer);
    }

    return true;
}

bool RLEByteEncoder::Deserialize(const std::vector<uint8_t>& buffer, RLEByteResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool RLEByteEncoder::Deserialize(const uint8_t* buffer, size_t size, RLEByteResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
    uint64_t num_runs = VarIntUtil::ReadVarUint(buffer, size, pos);

    result.runs.reserve(num_runs);
    for (uint64_t i = 0; i < num_runs; ++i) {
        RLEByteRun run;
        run.value = buffer[pos++];
        run.count = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
        result.runs.push_back(run);
    }

    return true;
}

// ============================================================================
// DeltaEncoder Implementation
// ============================================================================

bool DeltaEncoder::Compress(const std::vector<uint64_t>& data, DeltaResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());
    result.first_value = data[0];
    result.deltas.reserve(data.size() - 1);

    for (size_t i = 1; i < data.size(); ++i) {
        result.deltas.push_back(static_cast<int64_t>(data[i]) - static_cast<int64_t>(data[i - 1]));
    }

    return true;
}

bool DeltaEncoder::Decompress(const DeltaResult& result, std::vector<uint64_t>& data) {
    data.clear();
    if (result.original_len == 0) {
        return true;
    }

    data.reserve(result.original_len);
    data.push_back(result.first_value);

    uint64_t current = result.first_value;
    for (int64_t delta : result.deltas) {
        current = static_cast<uint64_t>(static_cast<int64_t>(current) + delta);
        data.push_back(current);
    }

    return data.size() == result.original_len;
}

bool DeltaEncoder::CompressSigned(const std::vector<int64_t>& data, DeltaResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());
    result.first_value = static_cast<uint64_t>(data[0]); // Interpret as unsigned for storage
    result.deltas.reserve(data.size() - 1);

    for (size_t i = 1; i < data.size(); ++i) {
        result.deltas.push_back(data[i] - data[i - 1]);
    }

    return true;
}

bool DeltaEncoder::DecompressSigned(const DeltaResult& result, std::vector<int64_t>& data) {
    data.clear();
    if (result.original_len == 0) {
        return true;
    }

    data.reserve(result.original_len);
    data.push_back(static_cast<int64_t>(result.first_value));

    int64_t current = static_cast<int64_t>(result.first_value);
    for (int64_t delta : result.deltas) {
        current += delta;
        data.push_back(current);
    }

    return data.size() == result.original_len;
}

bool DeltaEncoder::Serialize(const DeltaResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    VarIntUtil::WriteVarUint(result.original_len, buffer);
    VarIntUtil::WriteVarUint(result.first_value, buffer);

    // Write delta count
    VarIntUtil::WriteVarUint(result.deltas.size(), buffer);

    // Write deltas using zigzag encoding
    for (int64_t delta : result.deltas) {
        VarIntUtil::WriteVarInt(delta, buffer);
    }

    return true;
}

bool DeltaEncoder::Deserialize(const std::vector<uint8_t>& buffer, DeltaResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool DeltaEncoder::Deserialize(const uint8_t* buffer, size_t size, DeltaResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
    result.first_value = VarIntUtil::ReadVarUint(buffer, size, pos);

    uint64_t delta_count = VarIntUtil::ReadVarUint(buffer, size, pos);
    result.deltas.reserve(delta_count);

    for (uint64_t i = 0; i < delta_count; ++i) {
        result.deltas.push_back(VarIntUtil::ReadVarInt(buffer, size, pos));
    }

    return true;
}

// ============================================================================
// MaskEncoder Implementation
// ============================================================================

bool MaskEncoder::GetBit(const std::vector<uint8_t>& bitmask, uint32_t pos) {
    return BitmapUtil::GetBit(bitmask, pos);
}

void MaskEncoder::SetBit(std::vector<uint8_t>& bitmask, uint32_t pos, bool value) {
    BitmapUtil::SetBit(bitmask, pos, value);
}

bool MaskEncoder::Compress(const std::vector<std::string>& data, MaskResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    // Find dominant value (most frequent)
    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    std::string dominant;
    uint32_t max_count = 0;
    for (const auto& kv : freq) { const auto& val = kv.first; const auto& count = kv.second;
        if (count > max_count) {
            max_count = count;
            dominant = val;
        }
    }

    return CompressWithDominant(data, dominant, result);
}

bool MaskEncoder::CompressWithDominant(const std::vector<std::string>& data,
                                       const std::string& dominant,
                                       MaskResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());
    result.dominant_value = dominant;

    // Allocate bitmask
    size_t byte_count = (data.size() + 7) / 8;
    result.bitmask.resize(byte_count, 0);

    uint32_t dominant_count = 0;

    // Build bitmask and collect patches
    for (uint32_t i = 0; i < data.size(); ++i) {
        if (data[i] == dominant) {
            SetBit(result.bitmask, i, true);
            dominant_count++;
        } else {
            // SetBit(result.bitmask, i, false); // Already 0
            result.patches.push_back(data[i]);
            result.patch_indices.push_back(i);
        }
    }

    result.dominant_ratio = static_cast<float>(dominant_count) / data.size();

    return true;
}

bool MaskEncoder::Decompress(const MaskResult& result, std::vector<std::string>& data) {
    data.clear();
    if (result.original_len == 0) {
        return true;
    }

    data.resize(result.original_len, result.dominant_value);

    // Apply patches
    for (size_t i = 0; i < result.patches.size(); ++i) {
        data[result.patch_indices[i]] = result.patches[i];
    }

    return true;
}

bool MaskEncoder::Serialize(const MaskResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    // Header
    VarIntUtil::WriteVarUint(result.original_len, buffer);

    // Dominant value
    VarIntUtil::WriteVarUint(result.dominant_value.size(), buffer);
    buffer.insert(buffer.end(), result.dominant_value.begin(), result.dominant_value.end());

    // Bitmask (RLE compressed)
    RLEByteResult rle_mask;
    RLEByteEncoder::Compress(result.bitmask, rle_mask);
    std::vector<uint8_t> mask_buf;
    RLEByteEncoder::Serialize(rle_mask, mask_buf);

    VarIntUtil::WriteVarUint(mask_buf.size(), buffer);
    buffer.insert(buffer.end(), mask_buf.begin(), mask_buf.end());

    // Patches count
    VarIntUtil::WriteVarUint(result.patches.size(), buffer);

    // Patch indices (delta encoded)
    if (!result.patch_indices.empty()) {
        std::vector<uint64_t> indices_u64(result.patch_indices.begin(), result.patch_indices.end());
        DeltaResult delta_indices;
        DeltaEncoder::Compress(indices_u64, delta_indices);
        std::vector<uint8_t> idx_buf;
        DeltaEncoder::Serialize(delta_indices, idx_buf);

        VarIntUtil::WriteVarUint(idx_buf.size(), buffer);
        buffer.insert(buffer.end(), idx_buf.begin(), idx_buf.end());
    }

    // Patches (string values)
    for (const auto& patch : result.patches) {
        VarIntUtil::WriteVarUint(patch.size(), buffer);
        buffer.insert(buffer.end(), patch.begin(), patch.end());
    }

    return true;
}

bool MaskEncoder::Deserialize(const std::vector<uint8_t>& buffer, MaskResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool MaskEncoder::Deserialize(const uint8_t* buffer, size_t size, MaskResult& result) {
    result.clear();
    size_t pos = 0;

    // Header
    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));

    // Dominant value
    uint64_t dom_len = VarIntUtil::ReadVarUint(buffer, size, pos);
    result.dominant_value.assign(reinterpret_cast<const char*>(buffer + pos), dom_len);
    pos += dom_len;

    // Bitmask
    uint64_t mask_size = VarIntUtil::ReadVarUint(buffer, size, pos);
    RLEByteResult rle_mask;
    RLEByteEncoder::Deserialize(buffer + pos, mask_size, rle_mask);
    RLEByteEncoder::Decompress(rle_mask, result.bitmask);
    pos += mask_size;

    // Patches count
    uint64_t patch_count = VarIntUtil::ReadVarUint(buffer, size, pos);

    // Patch indices
    if (patch_count > 0) {
        uint64_t idx_size = VarIntUtil::ReadVarUint(buffer, size, pos);
        DeltaResult delta_indices;
        DeltaEncoder::Deserialize(buffer + pos, idx_size, delta_indices);
        std::vector<uint64_t> indices_u64;
        DeltaEncoder::Decompress(delta_indices, indices_u64);
        result.patch_indices.assign(indices_u64.begin(), indices_u64.end());
        pos += idx_size;
    }

    // Patches
    result.patches.reserve(patch_count);
    for (uint64_t i = 0; i < patch_count; ++i) {
        uint64_t str_len = VarIntUtil::ReadVarUint(buffer, size, pos);
        std::string patch(reinterpret_cast<const char*>(buffer + pos), str_len);
        result.patches.push_back(std::move(patch));
        pos += str_len;
    }

    // Calculate dominant ratio
    uint32_t dom_count = result.original_len - static_cast<uint32_t>(result.patches.size());
    result.dominant_ratio = result.original_len > 0 ?
        static_cast<float>(dom_count) / result.original_len : 0.0f;

    return true;
}

bool MaskEncoder::DeserializeWithPos(const std::vector<uint8_t>& buffer, size_t& pos, MaskResult& result) {
    result.clear();
    const uint8_t* buf = buffer.data();
    size_t size = buffer.size();

    // Header
    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buf, size, pos));

    // Dominant value
    uint64_t dom_len = VarIntUtil::ReadVarUint(buf, size, pos);
    result.dominant_value.assign(reinterpret_cast<const char*>(buf + pos), dom_len);
    pos += dom_len;

    // Bitmask
    uint64_t mask_size = VarIntUtil::ReadVarUint(buf, size, pos);
    RLEByteResult rle_mask;
    RLEByteEncoder::Deserialize(buf + pos, mask_size, rle_mask);
    RLEByteEncoder::Decompress(rle_mask, result.bitmask);
    pos += mask_size;

    // Patches count
    uint64_t patch_count = VarIntUtil::ReadVarUint(buf, size, pos);

    // Patch indices
    if (patch_count > 0) {
        uint64_t idx_size = VarIntUtil::ReadVarUint(buf, size, pos);
        DeltaResult delta_indices;
        DeltaEncoder::Deserialize(buf + pos, idx_size, delta_indices);
        std::vector<uint64_t> indices_u64;
        DeltaEncoder::Decompress(delta_indices, indices_u64);
        result.patch_indices.assign(indices_u64.begin(), indices_u64.end());
        pos += idx_size;
    }

    // Patches
    result.patches.reserve(patch_count);
    for (uint64_t i = 0; i < patch_count; ++i) {
        uint64_t str_len = VarIntUtil::ReadVarUint(buf, size, pos);
        std::string patch(reinterpret_cast<const char*>(buf + pos), str_len);
        result.patches.push_back(std::move(patch));
        pos += str_len;
    }

    // Calculate dominant ratio
    uint32_t dom_count = result.original_len - static_cast<uint32_t>(result.patches.size());
    result.dominant_ratio = result.original_len > 0 ?
        static_cast<float>(dom_count) / result.original_len : 0.0f;

    return true;
}

// ============================================================================
// MaskIntEncoder Implementation
// ============================================================================

bool MaskIntEncoder::Compress(const std::vector<int32_t>& data, MaskIntResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    // Find dominant value
    std::unordered_map<int32_t, uint32_t> freq;
    for (int32_t val : data) {
        freq[val]++;
    }

    int32_t dominant = 0;
    uint32_t max_count = 0;
    for (const auto& kv : freq) { const auto& val = kv.first; const auto& count = kv.second;
        if (count > max_count) {
            max_count = count;
            dominant = val;
        }
    }

    return CompressWithDominant(data, dominant, result);
}

bool MaskIntEncoder::CompressWithDominant(const std::vector<int32_t>& data,
                                          int32_t dominant,
                                          MaskIntResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());
    result.dominant_value = dominant;

    size_t byte_count = (data.size() + 7) / 8;
    result.bitmask.resize(byte_count, 0);

    uint32_t dominant_count = 0;

    for (uint32_t i = 0; i < data.size(); ++i) {
        if (data[i] == dominant) {
            BitmapUtil::SetBit(result.bitmask, i, true);
            dominant_count++;
        } else {
            result.patches.push_back(data[i]);
            result.patch_indices.push_back(i);
        }
    }

    result.dominant_ratio = static_cast<float>(dominant_count) / data.size();

    return true;
}

bool MaskIntEncoder::Decompress(const MaskIntResult& result, std::vector<int32_t>& data) {
    data.clear();
    if (result.original_len == 0) {
        return true;
    }

    data.resize(result.original_len, result.dominant_value);

    for (size_t i = 0; i < result.patches.size(); ++i) {
        data[result.patch_indices[i]] = result.patches[i];
    }

    return true;
}

bool MaskIntEncoder::Serialize(const MaskIntResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    VarIntUtil::WriteVarUint(result.original_len, buffer);
    VarIntUtil::WriteVarInt(result.dominant_value, buffer);

    // Bitmask (RLE compressed)
    RLEByteResult rle_mask;
    RLEByteEncoder::Compress(result.bitmask, rle_mask);
    std::vector<uint8_t> mask_buf;
    RLEByteEncoder::Serialize(rle_mask, mask_buf);

    VarIntUtil::WriteVarUint(mask_buf.size(), buffer);
    buffer.insert(buffer.end(), mask_buf.begin(), mask_buf.end());

    // Patches count
    VarIntUtil::WriteVarUint(result.patches.size(), buffer);

    // Patch indices (delta encoded)
    if (!result.patch_indices.empty()) {
        std::vector<uint64_t> indices_u64(result.patch_indices.begin(), result.patch_indices.end());
        DeltaResult delta_indices;
        DeltaEncoder::Compress(indices_u64, delta_indices);
        std::vector<uint8_t> idx_buf;
        DeltaEncoder::Serialize(delta_indices, idx_buf);

        VarIntUtil::WriteVarUint(idx_buf.size(), buffer);
        buffer.insert(buffer.end(), idx_buf.begin(), idx_buf.end());
    }

    // Patches (int values)
    for (int32_t patch : result.patches) {
        VarIntUtil::WriteVarInt(patch, buffer);
    }

    return true;
}

bool MaskIntEncoder::Deserialize(const std::vector<uint8_t>& buffer, MaskIntResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool MaskIntEncoder::Deserialize(const uint8_t* buffer, size_t size, MaskIntResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));
    result.dominant_value = static_cast<int32_t>(VarIntUtil::ReadVarInt(buffer, size, pos));

    // Bitmask
    uint64_t mask_size = VarIntUtil::ReadVarUint(buffer, size, pos);
    RLEByteResult rle_mask;
    RLEByteEncoder::Deserialize(buffer + pos, mask_size, rle_mask);
    RLEByteEncoder::Decompress(rle_mask, result.bitmask);
    pos += mask_size;

    // Patches count
    uint64_t patch_count = VarIntUtil::ReadVarUint(buffer, size, pos);

    // Patch indices
    if (patch_count > 0) {
        uint64_t idx_size = VarIntUtil::ReadVarUint(buffer, size, pos);
        DeltaResult delta_indices;
        DeltaEncoder::Deserialize(buffer + pos, idx_size, delta_indices);
        std::vector<uint64_t> indices_u64;
        DeltaEncoder::Decompress(delta_indices, indices_u64);
        result.patch_indices.assign(indices_u64.begin(), indices_u64.end());
        pos += idx_size;
    }

    // Patches
    result.patches.reserve(patch_count);
    for (uint64_t i = 0; i < patch_count; ++i) {
        result.patches.push_back(static_cast<int32_t>(VarIntUtil::ReadVarInt(buffer, size, pos)));
    }

    uint32_t dom_count = result.original_len - static_cast<uint32_t>(result.patches.size());
    result.dominant_ratio = result.original_len > 0 ?
        static_cast<float>(dom_count) / result.original_len : 0.0f;

    return true;
}

// ============================================================================
// DictEncoder Implementation
// ============================================================================

bool DictEncoder::Compress(const std::vector<std::string>& data, DictResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> value_to_index;

    for (const auto& val : data) {
        auto it = value_to_index.find(val);
        if (it == value_to_index.end()) {
            uint32_t idx = static_cast<uint32_t>(result.dictionary.size());
            result.dictionary.push_back(val);
            value_to_index[val] = idx;
            result.indices.push_back(idx);
        } else {
            result.indices.push_back(it->second);
        }
    }

    return true;
}

bool DictEncoder::Decompress(const DictResult& result, std::vector<std::string>& data) {
    data.clear();
    if (result.original_len == 0) {
        return true;
    }

    data.reserve(result.original_len);
    for (uint32_t idx : result.indices) {
        if (idx >= result.dictionary.size()) {
            return false; // Invalid index
        }
        data.push_back(result.dictionary[idx]);
    }

    return data.size() == result.original_len;
}

bool DictEncoder::Serialize(const DictResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    VarIntUtil::WriteVarUint(result.original_len, buffer);

    // Dictionary
    VarIntUtil::WriteVarUint(result.dictionary.size(), buffer);
    for (const auto& val : result.dictionary) {
        VarIntUtil::WriteVarUint(val.size(), buffer);
        buffer.insert(buffer.end(), val.begin(), val.end());
    }

    // Indices (RLE compressed if beneficial)
    RLEIntResult rle_indices;
    std::vector<int32_t> indices_signed(result.indices.begin(), result.indices.end());
    RLEIntEncoder::Compress(indices_signed, rle_indices);

    // Compare RLE vs raw size
    std::vector<uint8_t> rle_buf;
    RLEIntEncoder::Serialize(rle_indices, rle_buf);

    size_t raw_size = 0;
    for (uint32_t idx : result.indices) {
        raw_size += (idx < 128) ? 1 : ((idx < 16384) ? 2 : 3);
    }

    if (rle_buf.size() < raw_size) {
        // Use RLE
        buffer.push_back(1); // RLE flag
        VarIntUtil::WriteVarUint(rle_buf.size(), buffer);
        buffer.insert(buffer.end(), rle_buf.begin(), rle_buf.end());
    } else {
        // Use raw
        buffer.push_back(0); // Raw flag
        VarIntUtil::WriteVarUint(result.indices.size(), buffer);
        for (uint32_t idx : result.indices) {
            VarIntUtil::WriteVarUint(idx, buffer);
        }
    }

    return true;
}

bool DictEncoder::Deserialize(const std::vector<uint8_t>& buffer, DictResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool DictEncoder::Deserialize(const uint8_t* buffer, size_t size, DictResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));

    // Dictionary
    uint64_t dict_size = VarIntUtil::ReadVarUint(buffer, size, pos);
    result.dictionary.reserve(dict_size);
    for (uint64_t i = 0; i < dict_size; ++i) {
        uint64_t str_len = VarIntUtil::ReadVarUint(buffer, size, pos);
        std::string val(reinterpret_cast<const char*>(buffer + pos), str_len);
        result.dictionary.push_back(std::move(val));
        pos += str_len;
    }

    // Indices
    uint8_t use_rle = buffer[pos++];
    if (use_rle) {
        uint64_t rle_size = VarIntUtil::ReadVarUint(buffer, size, pos);
        RLEIntResult rle_indices;
        RLEIntEncoder::Deserialize(buffer + pos, rle_size, rle_indices);
        std::vector<int32_t> indices_signed;
        RLEIntEncoder::Decompress(rle_indices, indices_signed);
        result.indices.assign(indices_signed.begin(), indices_signed.end());
        pos += rle_size;
    } else {
        uint64_t idx_count = VarIntUtil::ReadVarUint(buffer, size, pos);
        result.indices.reserve(idx_count);
        for (uint64_t i = 0; i < idx_count; ++i) {
            result.indices.push_back(static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos)));
        }
    }

    return true;
}

// ============================================================================
// DictIntEncoder Implementation
// ============================================================================

bool DictIntEncoder::Compress(const std::vector<int32_t>& data, DictIntResult& result) {
    result.clear();
    if (data.empty()) {
        return true;
    }

    result.original_len = static_cast<uint32_t>(data.size());

    std::unordered_map<int32_t, uint32_t> value_to_index;

    for (int32_t val : data) {
        auto it = value_to_index.find(val);
        if (it == value_to_index.end()) {
            uint32_t idx = static_cast<uint32_t>(result.dictionary.size());
            result.dictionary.push_back(val);
            value_to_index[val] = idx;
            result.indices.push_back(idx);
        } else {
            result.indices.push_back(it->second);
        }
    }

    return true;
}

bool DictIntEncoder::Decompress(const DictIntResult& result, std::vector<int32_t>& data) {
    data.clear();
    if (result.original_len == 0) {
        return true;
    }

    data.reserve(result.original_len);
    for (uint32_t idx : result.indices) {
        if (idx >= result.dictionary.size()) {
            return false;
        }
        data.push_back(result.dictionary[idx]);
    }

    return data.size() == result.original_len;
}

bool DictIntEncoder::Serialize(const DictIntResult& result, std::vector<uint8_t>& buffer) {
    buffer.clear();

    VarIntUtil::WriteVarUint(result.original_len, buffer);

    // Dictionary
    VarIntUtil::WriteVarUint(result.dictionary.size(), buffer);
    for (int32_t val : result.dictionary) {
        VarIntUtil::WriteVarInt(val, buffer);
    }

    // Indices (RLE if beneficial)
    RLEIntResult rle_indices;
    std::vector<int32_t> indices_signed(result.indices.begin(), result.indices.end());
    RLEIntEncoder::Compress(indices_signed, rle_indices);

    std::vector<uint8_t> rle_buf;
    RLEIntEncoder::Serialize(rle_indices, rle_buf);

    size_t raw_size = 0;
    for (uint32_t idx : result.indices) {
        raw_size += (idx < 128) ? 1 : ((idx < 16384) ? 2 : 3);
    }

    if (rle_buf.size() < raw_size) {
        buffer.push_back(1);
        VarIntUtil::WriteVarUint(rle_buf.size(), buffer);
        buffer.insert(buffer.end(), rle_buf.begin(), rle_buf.end());
    } else {
        buffer.push_back(0);
        VarIntUtil::WriteVarUint(result.indices.size(), buffer);
        for (uint32_t idx : result.indices) {
            VarIntUtil::WriteVarUint(idx, buffer);
        }
    }

    return true;
}

bool DictIntEncoder::Deserialize(const std::vector<uint8_t>& buffer, DictIntResult& result) {
    return Deserialize(buffer.data(), buffer.size(), result);
}

bool DictIntEncoder::Deserialize(const uint8_t* buffer, size_t size, DictIntResult& result) {
    result.clear();
    size_t pos = 0;

    result.original_len = static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos));

    // Dictionary
    uint64_t dict_size = VarIntUtil::ReadVarUint(buffer, size, pos);
    result.dictionary.reserve(dict_size);
    for (uint64_t i = 0; i < dict_size; ++i) {
        result.dictionary.push_back(static_cast<int32_t>(VarIntUtil::ReadVarInt(buffer, size, pos)));
    }

    // Indices
    uint8_t use_rle = buffer[pos++];
    if (use_rle) {
        uint64_t rle_size = VarIntUtil::ReadVarUint(buffer, size, pos);
        RLEIntResult rle_indices;
        RLEIntEncoder::Deserialize(buffer + pos, rle_size, rle_indices);
        std::vector<int32_t> indices_signed;
        RLEIntEncoder::Decompress(rle_indices, indices_signed);
        result.indices.assign(indices_signed.begin(), indices_signed.end());
        pos += rle_size;
    } else {
        uint64_t idx_count = VarIntUtil::ReadVarUint(buffer, size, pos);
        result.indices.reserve(idx_count);
        for (uint64_t i = 0; i < idx_count; ++i) {
            result.indices.push_back(static_cast<uint32_t>(VarIntUtil::ReadVarUint(buffer, size, pos)));
        }
    }

    return true;
}

} // namespace gvcf
