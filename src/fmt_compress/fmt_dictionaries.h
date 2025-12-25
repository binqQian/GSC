#pragma once

#include <vector>
#include <unordered_map>
#include <mutex>
#include <cstring>
#include "fmt_utils.h"
#include "vint_codec.h"

namespace fmt_compress {

// Global dictionaries for AD, PL, and PID items
// Thread-safe with mutex protection for writes
class FmtDictionaries {
public:
    FmtDictionaries() = default;
    ~FmtDictionaries() = default;

    // Get or create AD item ID (thread-safe)
    uint32_t getADItemId(const ADItem& item) {
        std::lock_guard<std::mutex> lock(ad_mutex_);
        auto it = ad_map_.find(item);
        if (it != ad_map_.end()) {
            return it->second;
        }

        uint32_t id = static_cast<uint32_t>(ad_ptrs_.size());
        ad_map_[item] = id;

        // Store item data
        size_t offset = ad_items_.size();
        uint8_t len = item.data_[0] + 1;
        ad_items_.insert(ad_items_.end(), item.data_, item.data_ + len);
        ad_ptrs_.push_back(ad_items_.data() + offset);

        return id;
    }

    // Get or create PL item ID (thread-safe)
    uint32_t getPLItemId(const PLItem& item) {
        std::lock_guard<std::mutex> lock(pl_mutex_);
        auto it = pl_map_.find(item);
        if (it != pl_map_.end()) {
            return it->second;
        }

        uint32_t id = static_cast<uint32_t>(pl_ptrs_.size());
        pl_map_[item] = id;

        // Store item data
        size_t offset = pl_items_.size();
        uint8_t len = item.data_[0] + 1;
        pl_items_.insert(pl_items_.end(), item.data_, item.data_ + len);
        pl_ptrs_.push_back(pl_items_.data() + offset);

        return id;
    }

    // Get or create PID item ID (thread-safe)
    uint32_t getPIDItemId(const PIDItem& item) {
        std::lock_guard<std::mutex> lock(pid_mutex_);
        auto it = pid_map_.find(item);
        if (it != pid_map_.end()) {
            return it->second;
        }

        uint32_t id = static_cast<uint32_t>(pid_ptrs_.size());
        pid_map_[item] = id;

        // Store item data
        size_t offset = pid_items_.size();
        uint16_t len;
        memcpy(&len, item.data_, 2);
        len += 2;
        pid_items_.insert(pid_items_.end(), item.data_, item.data_ + len);
        pid_ptrs_.push_back(pid_items_.data() + offset);

        return id;
    }

    // Get AD item by ID
    const uint8_t* getADItemPtr(uint32_t id) const {
        return (id < ad_ptrs_.size()) ? ad_ptrs_[id] : nullptr;
    }

    // Get PL item by ID
    const uint8_t* getPLItemPtr(uint32_t id) const {
        return (id < pl_ptrs_.size()) ? pl_ptrs_[id] : nullptr;
    }

    // Get PID item by ID
    const uint8_t* getPIDItemPtr(uint32_t id) const {
        return (id < pid_ptrs_.size()) ? pid_ptrs_[id] : nullptr;
    }

    // Serialize all dictionaries to output buffer
    void serialize(std::vector<uint8_t>& output) const {
        uint8_t buf[16];

        // AD dictionary
        uint8_t len = VintCodec::encode(ad_items_.size(), buf);
        output.insert(output.end(), buf, buf + len);
        output.insert(output.end(), ad_items_.begin(), ad_items_.end());

        // PL dictionary
        len = VintCodec::encode(pl_items_.size(), buf);
        output.insert(output.end(), buf, buf + len);
        output.insert(output.end(), pl_items_.begin(), pl_items_.end());

        // PID dictionary
        len = VintCodec::encode(pid_items_.size(), buf);
        output.insert(output.end(), buf, buf + len);
        output.insert(output.end(), pid_items_.begin(), pid_items_.end());
    }

    // Deserialize dictionaries from input buffer
    size_t unserialize(const uint8_t* data, size_t data_len) {
        size_t offset = 0;
        uint64_t size;

        // AD dictionary
        offset += VintCodec::decode(data + offset, size);
        ad_items_.assign(data + offset, data + offset + size);
        rebuildADPtrs();
        offset += size;

        // PL dictionary
        offset += VintCodec::decode(data + offset, size);
        pl_items_.assign(data + offset, data + offset + size);
        rebuildPLPtrs();
        offset += size;

        // PID dictionary
        offset += VintCodec::decode(data + offset, size);
        pid_items_.assign(data + offset, data + offset + size);
        rebuildPIDPtrs();
        offset += size;

        return offset;
    }

    // Clear all dictionaries
    void clear() {
        std::lock_guard<std::mutex> lock1(ad_mutex_);
        std::lock_guard<std::mutex> lock2(pl_mutex_);
        std::lock_guard<std::mutex> lock3(pid_mutex_);

        ad_map_.clear();
        ad_items_.clear();
        ad_ptrs_.clear();

        pl_map_.clear();
        pl_items_.clear();
        pl_ptrs_.clear();

        pid_map_.clear();
        pid_items_.clear();
        pid_ptrs_.clear();
    }

    // Get dictionary sizes for statistics
    size_t getADCount() const { return ad_ptrs_.size(); }
    size_t getPLCount() const { return pl_ptrs_.size(); }
    size_t getPIDCount() const { return pid_ptrs_.size(); }

private:
    // Rebuild AD pointers after deserialization
    void rebuildADPtrs() {
        ad_ptrs_.clear();
        size_t offset = 0;
        while (offset < ad_items_.size()) {
            ad_ptrs_.push_back(ad_items_.data() + offset);
            uint8_t len = ad_items_[offset] + 1;
            offset += len;
        }
    }

    // Rebuild PL pointers after deserialization
    void rebuildPLPtrs() {
        pl_ptrs_.clear();
        size_t offset = 0;
        while (offset < pl_items_.size()) {
            pl_ptrs_.push_back(pl_items_.data() + offset);
            uint8_t len = pl_items_[offset] + 1;
            offset += len;
        }
    }

    // Rebuild PID pointers after deserialization
    void rebuildPIDPtrs() {
        pid_ptrs_.clear();
        size_t offset = 0;
        while (offset < pid_items_.size()) {
            pid_ptrs_.push_back(pid_items_.data() + offset);
            uint16_t len;
            memcpy(&len, pid_items_.data() + offset, 2);
            offset += len + 2;
        }
    }

    // Mutexes for thread safety
    std::mutex ad_mutex_;
    std::mutex pl_mutex_;
    std::mutex pid_mutex_;

    // AD dictionary
    std::unordered_map<ADItem, uint32_t> ad_map_;
    std::vector<uint8_t> ad_items_;
    std::vector<const uint8_t*> ad_ptrs_;

    // PL dictionary
    std::unordered_map<PLItem, uint32_t> pl_map_;
    std::vector<uint8_t> pl_items_;
    std::vector<const uint8_t*> pl_ptrs_;

    // PID dictionary
    std::unordered_map<PIDItem, uint32_t> pid_map_;
    std::vector<uint8_t> pid_items_;
    std::vector<const uint8_t*> pid_ptrs_;
};

}  // namespace fmt_compress
