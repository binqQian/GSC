#pragma once

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <thread>
#include <string>
#include "gsc_params.h"
#include "memory_monitor.h"
#include "logger.h"

namespace gsc {

/**
 * System information collected at startup
 */
struct SystemInfo {
    uint32_t cpu_cores = 1;
    uint64_t available_memory_mb = 0;
    uint64_t total_memory_mb = 0;

    static SystemInfo Detect() {
        SystemInfo info;
        info.cpu_cores = std::thread::hardware_concurrency();
        if (info.cpu_cores == 0) info.cpu_cores = 1;
        info.available_memory_mb = MemoryMonitor::GetAvailableMemoryKB() / 1024;
        info.total_memory_mb = MemoryMonitor::GetTotalMemoryKB() / 1024;
        return info;
    }
};

/**
 * Workload characteristics for resource estimation
 */
struct WorkloadEstimate {
    uint32_t n_samples = 0;
    uint32_t ploidy = 2;
    uint32_t max_block_rows = 10000;
    uint32_t n_col_blocks = 1;
    compression_backend_t backend = compression_backend_t::bsc;
};

/**
 * Thread configuration result
 */
struct ThreadConfig {
    uint32_t no_threads = 4;       // General compression threads
    uint32_t no_gt_threads = 1;    // GT processing threads
    uint32_t no_parse_threads = 1; // VCF parsing threads
    size_t queue_capacity = 8;     // GT block queue capacity
};

/**
 * Memory estimate for a single block
 */
struct MemoryEstimate {
    size_t gt_data = 0;            // GT double buffer
    size_t variant_metadata = 0;   // variant_desc_t per row
    size_t fixed_fields = 0;       // Fixed field buffers
    size_t bit_vectors = 0;        // Bit vectors and rank support
    size_t compression_buffer = 0; // Temporary compression buffers

    size_t total() const {
        return gt_data + variant_metadata + fixed_fields +
               bit_vectors + compression_buffer;
    }

    size_t total_mb() const {
        return (total() + 1024 * 1024 - 1) / (1024 * 1024);
    }
};

/**
 * Memory budget allocation
 */
struct MemoryBudget {
    uint64_t total_budget_mb = 0;
    uint64_t backend_reserve_mb = 0;
    uint64_t queue_budget_mb = 0;
    uint64_t overhead_mb = 100;
    size_t per_block_estimate_mb = 0;
    bool user_specified = false;
};

/**
 * Resource Manager - Central resource allocation and scheduling
 *
 * Provides automatic configuration for threads and memory based on:
 * - System resources (CPU cores, available memory)
 * - Workload characteristics (samples, variants, compression backend)
 * - User overrides (command line parameters)
 *
 * Usage:
 *   ResourceManager rm;
 *   rm.SetSystemInfo(SystemInfo::Detect());
 *   rm.SetWorkload(workload);
 *   rm.SetUserLimits(params);
 *   auto config = rm.Configure();
 */
class ResourceManager {
public:
    ResourceManager() = default;

    // Set system information
    void SetSystemInfo(const SystemInfo& info) {
        sys_info_ = info;
    }

    // Set workload characteristics
    void SetWorkload(const WorkloadEstimate& work) {
        workload_ = work;
    }

    // Set user limits from params (0 = auto)
    void SetUserLimits(const GSC_Params& params) {
        user_no_threads_ = params.no_threads;
        user_no_gt_threads_ = params.no_gt_threads;
        user_no_parse_threads_ = params.no_parse_threads;
        user_queue_capacity_ = params.queue_capacity;
        user_max_memory_mb_ = params.max_memory_mb;
    }

    // Calculate and return optimal configuration
    ThreadConfig Configure();

    // Get the calculated memory budget
    const MemoryBudget& GetMemoryBudget() const { return memory_budget_; }

    // Get the block memory estimate
    const MemoryEstimate& GetBlockEstimate() const { return block_estimate_; }

    // Log the resource plan
    void LogResourcePlan() const;

    // Static helper: estimate backend memory requirement
    static size_t GetBackendMemoryReserve(compression_backend_t backend);

    // Static helper: estimate memory per block
    static MemoryEstimate EstimateBlockMemory(
        size_t n_samples,
        size_t ploidy,
        size_t max_block_rows,
        size_t n_col_blocks);

private:
    void CalculateMemoryBudget();
    ThreadConfig AutoConfigureThreads();

    SystemInfo sys_info_;
    WorkloadEstimate workload_;
    MemoryBudget memory_budget_;
    MemoryEstimate block_estimate_;
    ThreadConfig final_config_;

    // User overrides (0 = auto)
    uint32_t user_no_threads_ = 0;
    uint32_t user_no_gt_threads_ = 0;
    uint32_t user_no_parse_threads_ = 0;
    uint32_t user_queue_capacity_ = 0;
    uint64_t user_max_memory_mb_ = 0;
};

} // namespace gsc
