#pragma once

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <fstream>

#ifdef __linux__
#include <unistd.h>
#include <sys/resource.h>
#endif

namespace gsc {

/**
 * Memory monitoring utilities for Linux systems.
 * Provides functions to query current memory usage and system limits.
 */
class MemoryMonitor {
public:
    struct MemoryInfo {
        uint64_t virtual_size_kb = 0;    // Virtual memory size (VmSize)
        uint64_t resident_size_kb = 0;   // Resident set size (VmRSS)
        uint64_t shared_kb = 0;          // Shared memory
        uint64_t data_kb = 0;            // Data segment size (VmData)
        uint64_t peak_resident_kb = 0;   // Peak resident size (VmHWM)
    };

    /**
     * Get current process memory usage from /proc/self/status
     */
    static MemoryInfo GetProcessMemory() {
        MemoryInfo info;
#ifdef __linux__
        std::ifstream status("/proc/self/status");
        if (!status.is_open()) return info;

        std::string line;
        while (std::getline(status, line)) {
            if (line.compare(0, 7, "VmSize:") == 0) {
                sscanf(line.c_str(), "VmSize: %lu", &info.virtual_size_kb);
            } else if (line.compare(0, 6, "VmRSS:") == 0) {
                sscanf(line.c_str(), "VmRSS: %lu", &info.resident_size_kb);
            } else if (line.compare(0, 6, "VmHWM:") == 0) {
                sscanf(line.c_str(), "VmHWM: %lu", &info.peak_resident_kb);
            } else if (line.compare(0, 7, "VmData:") == 0) {
                sscanf(line.c_str(), "VmData: %lu", &info.data_kb);
            } else if (line.compare(0, 8, "RssAnon:") == 0) {
                // Shared is approximated as RSS - RssAnon
                uint64_t anon = 0;
                sscanf(line.c_str(), "RssAnon: %lu", &anon);
                info.shared_kb = info.resident_size_kb > anon ?
                                 info.resident_size_kb - anon : 0;
            }
        }
#endif
        return info;
    }

    /**
     * Get available system memory in KB
     */
    static uint64_t GetAvailableMemoryKB() {
#ifdef __linux__
        std::ifstream meminfo("/proc/meminfo");
        if (!meminfo.is_open()) return 0;

        std::string line;
        uint64_t available = 0;
        while (std::getline(meminfo, line)) {
            if (line.compare(0, 13, "MemAvailable:") == 0) {
                sscanf(line.c_str(), "MemAvailable: %lu", &available);
                break;
            }
        }
        return available;
#else
        return 0;
#endif
    }

    /**
     * Get total system memory in KB
     */
    static uint64_t GetTotalMemoryKB() {
#ifdef __linux__
        std::ifstream meminfo("/proc/meminfo");
        if (!meminfo.is_open()) return 0;

        std::string line;
        uint64_t total = 0;
        while (std::getline(meminfo, line)) {
            if (line.compare(0, 9, "MemTotal:") == 0) {
                sscanf(line.c_str(), "MemTotal: %lu", &total);
                break;
            }
        }
        return total;
#else
        return 0;
#endif
    }

    /**
     * Check if current memory usage exceeds the limit
     * @param limit_mb Memory limit in MB (0 = no limit)
     * @return true if over limit
     */
    static bool IsOverMemoryLimit(uint64_t limit_mb) {
        if (limit_mb == 0) return false;
        MemoryInfo info = GetProcessMemory();
        return (info.resident_size_kb / 1024) > limit_mb;
    }

    /**
     * Get memory usage as a formatted string
     */
    static std::string GetMemoryUsageString() {
        MemoryInfo info = GetProcessMemory();
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RSS: %lu MB, Peak: %lu MB, Virtual: %lu MB",
                 info.resident_size_kb / 1024,
                 info.peak_resident_kb / 1024,
                 info.virtual_size_kb / 1024);
        return std::string(buf);
    }

    /**
     * Get resident set size in MB
     */
    static uint64_t GetRSSMB() {
        MemoryInfo info = GetProcessMemory();
        return info.resident_size_kb / 1024;
    }

    /**
     * Apply a hard address-space limit for the current process.
     * Note: this limits virtual memory (RLIMIT_AS); allocations may fail with std::bad_alloc.
     */
    static bool ApplyMemoryLimitMB(uint64_t limit_mb) {
#ifdef __linux__
        if (limit_mb == 0) return true;
        struct rlimit rl;
        rl.rlim_cur = rl.rlim_max = static_cast<rlim_t>(limit_mb) * 1024 * 1024;
        return setrlimit(RLIMIT_AS, &rl) == 0;
#else
        (void)limit_mb;
        return false;
#endif
    }

    /**
     * Calculate recommended queue capacity based on available memory
     * @param block_size_bytes Size of each block in bytes
     * @param max_memory_mb Maximum memory to use (0 = auto)
     * @param memory_fraction Fraction of available memory to use (0.0-1.0)
     */
    static size_t CalculateQueueCapacity(size_t block_size_bytes,
                                          uint64_t max_memory_mb = 0,
                                          double memory_fraction = 0.3) {
        uint64_t available_kb;
        if (max_memory_mb > 0) {
            available_kb = max_memory_mb * 1024;
        } else {
            available_kb = GetAvailableMemoryKB();
            if (available_kb == 0) {
                // Fallback if we can't read memory info
                return 8;
            }
        }

        // Use specified fraction of available memory
        uint64_t queue_memory_kb = static_cast<uint64_t>(available_kb * memory_fraction);
        size_t block_size_kb = (block_size_bytes + 1023) / 1024;

        if (block_size_kb == 0) return 8;

        size_t capacity = queue_memory_kb / block_size_kb;
        return std::max((size_t)4, std::min(capacity, (size_t)64));
    }
};

} // namespace gsc
