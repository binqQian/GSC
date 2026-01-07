/**
 * @file gvcf_encoding.h
 * @brief gVCF encoding algorithms: RLE, Delta, Mask, Dictionary
 *
 * Reference: ref_code/gvcf_refCode/encoding/
 * These algorithms are optimized for gVCF data characteristics.
 */
#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <algorithm>

namespace gvcf {

// ============================================================================
// RLE (Run-Length Encoding)
// Used for: CHROM, FILTER, bitmasks - fields with high repetition
// ============================================================================

struct RLERun {
    std::string value;
    uint32_t count;
};

struct RLEResult {
    std::vector<RLERun> runs;
    uint32_t original_len;

    void clear() {
        runs.clear();
        original_len = 0;
    }
};

// RLE for string arrays
class RLEEncoder {
public:
    static bool Compress(const std::vector<std::string>& data, RLEResult& result);
    static bool Decompress(const RLEResult& result, std::vector<std::string>& data);
    static bool Serialize(const RLEResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, RLEResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, RLEResult& result);
};

// RLE for integer arrays (used for bitmasks, indices)
struct RLEIntRun {
    int32_t value;
    uint32_t count;
};

struct RLEIntResult {
    std::vector<RLEIntRun> runs;
    uint32_t original_len;

    void clear() {
        runs.clear();
        original_len = 0;
    }
};

class RLEIntEncoder {
public:
    static bool Compress(const std::vector<int32_t>& data, RLEIntResult& result);
    static bool Decompress(const RLEIntResult& result, std::vector<int32_t>& data);
    static bool Serialize(const RLEIntResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, RLEIntResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, RLEIntResult& result);
};

// RLE for byte arrays (used for bitmasks)
struct RLEByteRun {
    uint8_t value;
    uint32_t count;
};

struct RLEByteResult {
    std::vector<RLEByteRun> runs;
    uint32_t original_len;

    void clear() {
        runs.clear();
        original_len = 0;
    }
};

class RLEByteEncoder {
public:
    static bool Compress(const std::vector<uint8_t>& data, RLEByteResult& result);
    static bool Decompress(const RLEByteResult& result, std::vector<uint8_t>& data);
    static bool Serialize(const RLEByteResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, RLEByteResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, RLEByteResult& result);
};

// ============================================================================
// Delta Encoding
// Used for: POS, END - monotonically increasing sequences
// ============================================================================

struct DeltaResult {
    uint64_t first_value;        // First value in sequence
    std::vector<int64_t> deltas; // Delta values (can be negative for END field)
    uint32_t original_len;

    void clear() {
        first_value = 0;
        deltas.clear();
        original_len = 0;
    }
};

class DeltaEncoder {
public:
    // For unsigned sequences (POS)
    static bool Compress(const std::vector<uint64_t>& data, DeltaResult& result);
    static bool Decompress(const DeltaResult& result, std::vector<uint64_t>& data);

    // For signed sequences (END - POS differences)
    static bool CompressSigned(const std::vector<int64_t>& data, DeltaResult& result);
    static bool DecompressSigned(const DeltaResult& result, std::vector<int64_t>& data);

    static bool Serialize(const DeltaResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, DeltaResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, DeltaResult& result);
};

// ============================================================================
// Mask Encoding (Dominant Value Encoding)
// Used for: GT (0/0 dominant), ALT (<NON_REF> dominant), ID ('.' dominant)
// ============================================================================

struct MaskResult {
    std::string dominant_value;       // Most frequent value
    std::vector<uint8_t> bitmask;     // 1 = dominant, 0 = patch
    std::vector<std::string> patches; // Non-dominant values
    std::vector<uint32_t> patch_indices; // Positions of patches
    uint32_t original_len;
    float dominant_ratio;             // Ratio of dominant values

    void clear() {
        dominant_value.clear();
        bitmask.clear();
        patches.clear();
        patch_indices.clear();
        original_len = 0;
        dominant_ratio = 0.0f;
    }
};

class MaskEncoder {
public:
    // Compress string array with dominant value detection
    static bool Compress(const std::vector<std::string>& data, MaskResult& result);

    // Compress with specified dominant value (for known patterns like "0/0", "<NON_REF>")
    static bool CompressWithDominant(const std::vector<std::string>& data,
                                     const std::string& dominant,
                                     MaskResult& result);

    static bool Decompress(const MaskResult& result, std::vector<std::string>& data);

    static bool Serialize(const MaskResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, MaskResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, MaskResult& result);
    // Deserialize and return end position (for reading trailing data)
    static bool DeserializeWithPos(const std::vector<uint8_t>& buffer, size_t& pos, MaskResult& result);

    // Get bit at position
    static bool GetBit(const std::vector<uint8_t>& bitmask, uint32_t pos);
    // Set bit at position
    static void SetBit(std::vector<uint8_t>& bitmask, uint32_t pos, bool value);
};

// Mask encoding for integers (GT allele indices)
struct MaskIntResult {
    int32_t dominant_value;
    std::vector<uint8_t> bitmask;
    std::vector<int32_t> patches;
    std::vector<uint32_t> patch_indices;
    uint32_t original_len;
    float dominant_ratio;

    void clear() {
        dominant_value = 0;
        bitmask.clear();
        patches.clear();
        patch_indices.clear();
        original_len = 0;
        dominant_ratio = 0.0f;
    }
};

class MaskIntEncoder {
public:
    static bool Compress(const std::vector<int32_t>& data, MaskIntResult& result);
    static bool CompressWithDominant(const std::vector<int32_t>& data,
                                     int32_t dominant,
                                     MaskIntResult& result);
    static bool Decompress(const MaskIntResult& result, std::vector<int32_t>& data);
    static bool Serialize(const MaskIntResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, MaskIntResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, MaskIntResult& result);
};

// ============================================================================
// Dictionary Encoding
// Used for: REF, ALT variants, unknown fields with moderate cardinality
// ============================================================================

struct DictResult {
    std::vector<std::string> dictionary;   // Unique values
    std::vector<uint32_t> indices;         // Index into dictionary for each position
    uint32_t original_len;

    void clear() {
        dictionary.clear();
        indices.clear();
        original_len = 0;
    }
};

class DictEncoder {
public:
    static bool Compress(const std::vector<std::string>& data, DictResult& result);
    static bool Decompress(const DictResult& result, std::vector<std::string>& data);
    static bool Serialize(const DictResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, DictResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, DictResult& result);
};

// Dictionary for integers
struct DictIntResult {
    std::vector<int32_t> dictionary;
    std::vector<uint32_t> indices;
    uint32_t original_len;

    void clear() {
        dictionary.clear();
        indices.clear();
        original_len = 0;
    }
};

class DictIntEncoder {
public:
    static bool Compress(const std::vector<int32_t>& data, DictIntResult& result);
    static bool Decompress(const DictIntResult& result, std::vector<int32_t>& data);
    static bool Serialize(const DictIntResult& result, std::vector<uint8_t>& buffer);
    static bool Deserialize(const std::vector<uint8_t>& buffer, DictIntResult& result);
    static bool Deserialize(const uint8_t* buffer, size_t size, DictIntResult& result);
};

// ============================================================================
// Bitmap utilities
// Used for: presence/absence markers, phase information
// ============================================================================

class BitmapUtil {
public:
    // Create bitmap from boolean vector
    static std::vector<uint8_t> FromBools(const std::vector<bool>& data);

    // Extract boolean vector from bitmap
    static std::vector<bool> ToBools(const std::vector<uint8_t>& bitmap, uint32_t len);

    // Compress bitmap with RLE
    static bool CompressRLE(const std::vector<uint8_t>& bitmap, RLEByteResult& result);

    // Get/Set bit operations
    static bool GetBit(const std::vector<uint8_t>& bitmap, uint32_t pos);
    static void SetBit(std::vector<uint8_t>& bitmap, uint32_t pos, bool value);

    // Count set bits
    static uint32_t PopCount(const std::vector<uint8_t>& bitmap);
    static uint32_t PopCount(const std::vector<uint8_t>& bitmap, uint32_t len);
};

// ============================================================================
// Variable-length integer encoding utilities
// ============================================================================

class VarIntUtil {
public:
    // Write variable-length unsigned integer
    static size_t WriteVarUint(uint64_t value, std::vector<uint8_t>& buffer);
    static size_t WriteVarUint(uint64_t value, uint8_t* buffer);

    // Read variable-length unsigned integer
    static uint64_t ReadVarUint(const std::vector<uint8_t>& buffer, size_t& pos);
    static uint64_t ReadVarUint(const uint8_t* buffer, size_t size, size_t& pos);

    // Write variable-length signed integer (zigzag encoding)
    static size_t WriteVarInt(int64_t value, std::vector<uint8_t>& buffer);

    // Read variable-length signed integer
    static int64_t ReadVarInt(const std::vector<uint8_t>& buffer, size_t& pos);
    static int64_t ReadVarInt(const uint8_t* buffer, size_t size, size_t& pos);

    // Zigzag encoding for signed integers
    static uint64_t ZigzagEncode(int64_t value);
    static int64_t ZigzagDecode(uint64_t value);
};

} // namespace gvcf
