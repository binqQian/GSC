#include "integrity_checker.h"
#include "logger.h"

namespace gsc {

//==============================================================================
// IncrementalHasher Implementation
//==============================================================================

IncrementalHasher::IncrementalHasher(HashAlgorithm algo) : algorithm_(algo) {
    if (algorithm_ == HashAlgorithm::XXHASH64) {
        xxh64_state_ = XXH64_createState();
        if (xxh64_state_) {
            XXH64_reset(xxh64_state_, 0);
        }
    }
}

IncrementalHasher::~IncrementalHasher() {
    if (xxh64_state_) {
        XXH64_freeState(xxh64_state_);
        xxh64_state_ = nullptr;
    }
}

IncrementalHasher::IncrementalHasher(IncrementalHasher&& other) noexcept
    : algorithm_(other.algorithm_), xxh64_state_(other.xxh64_state_) {
    other.xxh64_state_ = nullptr;
}

IncrementalHasher& IncrementalHasher::operator=(IncrementalHasher&& other) noexcept {
    if (this != &other) {
        if (xxh64_state_) {
            XXH64_freeState(xxh64_state_);
        }
        algorithm_ = other.algorithm_;
        xxh64_state_ = other.xxh64_state_;
        other.xxh64_state_ = nullptr;
    }
    return *this;
}

void IncrementalHasher::Reset() {
    if (algorithm_ == HashAlgorithm::XXHASH64 && xxh64_state_) {
        XXH64_reset(xxh64_state_, 0);
    }
}

void IncrementalHasher::Update(const void* data, size_t length) {
    if (length == 0 || !data) return;

    if (algorithm_ == HashAlgorithm::XXHASH64 && xxh64_state_) {
        XXH64_update(xxh64_state_, data, length);
    }
}

void IncrementalHasher::Update(const std::vector<uint8_t>& data) {
    Update(data.data(), data.size());
}

void IncrementalHasher::Update(const std::string& data) {
    Update(data.data(), data.size());
}

HashResult IncrementalHasher::Finalize() {
    HashResult result;
    result.algorithm = algorithm_;

    if (algorithm_ == HashAlgorithm::XXHASH64 && xxh64_state_) {
        result.xxh64_value = XXH64_digest(xxh64_state_);
        result.valid = true;
    }

    return result;
}

HashResult IncrementalHasher::Compute(const void* data, size_t length, HashAlgorithm algo) {
    HashResult result;
    result.algorithm = algo;

    if (algo == HashAlgorithm::XXHASH64) {
        result.xxh64_value = XXH64(data, length, 0);
        result.valid = true;
    }

    return result;
}

//==============================================================================
// IntegrityManager Implementation
//==============================================================================

void IntegrityManager::EnableCompression(HashAlgorithm algo) {
    enabled_ = true;
    algorithm_ = algo;
    is_compression_mode_ = true;
    hasher_ = std::make_unique<IncrementalHasher>(algo);
}

void IntegrityManager::UpdateRawData(const void* data, size_t length) {
    if (enabled_ && hasher_ && is_compression_mode_) {
        hasher_->Update(data, length);
    }
}

void IntegrityManager::UpdateRawData(const std::vector<uint8_t>& data) {
    UpdateRawData(data.data(), data.size());
}

void IntegrityManager::UpdateRawData(const std::string& data) {
    UpdateRawData(data.data(), data.size());
}

HashResult IntegrityManager::FinalizeCompression() {
    if (!enabled_ || !hasher_ || !is_compression_mode_) {
        return HashResult{};
    }
    return hasher_->Finalize();
}

void IntegrityManager::EnableDecompression(const HashResult& expected) {
    enabled_ = true;
    algorithm_ = expected.algorithm;
    is_compression_mode_ = false;
    expected_hash_ = expected;
    hasher_ = std::make_unique<IncrementalHasher>(expected.algorithm);
}

void IntegrityManager::UpdateDecompressedData(const void* data, size_t length) {
    if (enabled_ && hasher_ && !is_compression_mode_) {
        hasher_->Update(data, length);
    }
}

void IntegrityManager::UpdateDecompressedData(const std::vector<uint8_t>& data) {
    UpdateDecompressedData(data.data(), data.size());
}

void IntegrityManager::UpdateDecompressedData(const std::string& data) {
    UpdateDecompressedData(data.data(), data.size());
}

HashResult IntegrityManager::GetComputedHash() {
    if (!enabled_ || !hasher_) {
        return HashResult{};
    }
    return hasher_->Finalize();
}

bool IntegrityManager::VerifyDecompression() {
    if (!enabled_ || !hasher_ || is_compression_mode_) {
        return false;
    }

    auto logger = LogManager::Instance().Logger();
    HashResult computed = hasher_->Finalize();

    bool match = (computed == expected_hash_);
    if (match) {
        logger->info("Integrity verification PASSED: {}", computed.ToHexString());
    } else {
        logger->error("Integrity verification FAILED!");
        logger->error("  Expected: {}", expected_hash_.ToHexString());
        logger->error("  Computed: {}", computed.ToHexString());
    }

    return match;
}

} // namespace gsc
