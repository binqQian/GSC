#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <iomanip>
#include <cstring>

// Include XXHash from zstd library
#define XXH_STATIC_LINKING_ONLY
#include "zstd-1.5.2/lib/common/xxhash.h"

namespace gsc {

// Hash algorithm type enumeration
enum class HashAlgorithm : uint8_t {
    NONE = 0,
    XXHASH64 = 1
};

// File format version constants for backward compatibility
static constexpr uint32_t GSC_FILE_FORMAT_VERSION_V1 = 1u;  // Original format (no hash)
static constexpr uint32_t GSC_FILE_FORMAT_VERSION_V2 = 2u;  // With hash support
static constexpr uint32_t GSC_FILE_FORMAT_VERSION_LATEST = GSC_FILE_FORMAT_VERSION_V2;

// Integrity flag bit positions
static constexpr uint8_t GSC_INTEGRITY_FLAG_HAS_HASH = 0x01;
static constexpr uint8_t GSC_INTEGRITY_ALGO_MASK = 0x06;     // bits 1-2
static constexpr uint8_t GSC_INTEGRITY_ALGO_NONE = 0x00;
static constexpr uint8_t GSC_INTEGRITY_ALGO_XXHASH64 = 0x02;

// Hash result structure (fixed size for serialization)
struct HashResult {
    HashAlgorithm algorithm = HashAlgorithm::NONE;
    uint64_t xxh64_value = 0;
    bool valid = false;

    // Serialize to bytes for file storage
    std::vector<uint8_t> Serialize() const {
        std::vector<uint8_t> result;
        result.reserve(SerializedSize());

        result.push_back(static_cast<uint8_t>(algorithm));

        // XXHash64 value (little-endian)
        for (int i = 0; i < 8; ++i) {
            result.push_back(static_cast<uint8_t>((xxh64_value >> (i * 8)) & 0xFF));
        }

        result.push_back(valid ? 1 : 0);

        return result;
    }

    // Deserialize from bytes
    bool Deserialize(const uint8_t* data, size_t size) {
        if (size < SerializedSize()) return false;

        size_t pos = 0;
        algorithm = static_cast<HashAlgorithm>(data[pos++]);

        xxh64_value = 0;
        for (int i = 0; i < 8; ++i) {
            xxh64_value |= static_cast<uint64_t>(data[pos++]) << (i * 8);
        }

        valid = (data[pos++] != 0);

        return true;
    }

    // Fixed serialized size: 1 (algo) + 8 (hash) + 1 (valid) = 10 bytes
    static constexpr size_t SerializedSize() { return 10; }

    // Comparison operator
    bool operator==(const HashResult& other) const {
        if (algorithm != other.algorithm) return false;
        if (!valid || !other.valid) return false;

        switch (algorithm) {
            case HashAlgorithm::XXHASH64:
                return xxh64_value == other.xxh64_value;
            default:
                return false;
        }
    }

    bool operator!=(const HashResult& other) const { return !(*this == other); }

    // Convert to hex string for logging
    std::string ToHexString() const {
        std::ostringstream oss;

        switch (algorithm) {
            case HashAlgorithm::XXHASH64:
                oss << "xxh64:" << std::hex << std::setfill('0')
                    << std::setw(16) << xxh64_value;
                break;
            default:
                oss << "none";
        }

        return oss.str();
    }
};

// Incremental hash calculator (supports streaming for large files)
class IncrementalHasher {
public:
    explicit IncrementalHasher(HashAlgorithm algo = HashAlgorithm::XXHASH64);
    ~IncrementalHasher();

    // Disable copy
    IncrementalHasher(const IncrementalHasher&) = delete;
    IncrementalHasher& operator=(const IncrementalHasher&) = delete;

    // Enable move
    IncrementalHasher(IncrementalHasher&& other) noexcept;
    IncrementalHasher& operator=(IncrementalHasher&& other) noexcept;

    // Core interface
    void Reset();                                    // Reset state
    void Update(const void* data, size_t length);   // Incremental update
    void Update(const std::vector<uint8_t>& data);  // Convenience interface
    void Update(const std::string& data);           // String interface
    HashResult Finalize();                          // Get final result

    // One-shot computation (static method)
    static HashResult Compute(const void* data, size_t length,
                              HashAlgorithm algo = HashAlgorithm::XXHASH64);

private:
    HashAlgorithm algorithm_;
    XXH64_state_t* xxh64_state_ = nullptr;
};

// Integrity verification manager (unified interface)
class IntegrityManager {
public:
    IntegrityManager() = default;
    ~IntegrityManager() = default;

    // For compression
    void EnableCompression(HashAlgorithm algo);
    void UpdateRawData(const void* data, size_t length);
    void UpdateRawData(const std::vector<uint8_t>& data);
    void UpdateRawData(const std::string& data);
    HashResult FinalizeCompression();

    // For decompression
    void EnableDecompression(const HashResult& expected);
    void UpdateDecompressedData(const void* data, size_t length);
    void UpdateDecompressedData(const std::vector<uint8_t>& data);
    void UpdateDecompressedData(const std::string& data);
    bool VerifyDecompression();  // Returns verification result
    HashResult GetComputedHash();  // Get computed hash for logging

    // State query
    bool IsEnabled() const { return enabled_; }
    HashAlgorithm GetAlgorithm() const { return algorithm_; }
    const HashResult& GetExpectedHash() const { return expected_hash_; }

private:
    bool enabled_ = false;
    HashAlgorithm algorithm_ = HashAlgorithm::NONE;
    std::unique_ptr<IncrementalHasher> hasher_;
    HashResult expected_hash_;
    bool is_compression_mode_ = true;
};

} // namespace gsc
