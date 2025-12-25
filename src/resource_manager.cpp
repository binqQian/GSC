#include "resource_manager.h"
#include <algorithm>

namespace gsc {

size_t ResourceManager::GetBackendMemoryReserve(compression_backend_t backend) {
    switch (backend) {
        case compression_backend_t::bsc:
            // BSC: BWT needs ~5x block size, plus LZP hash table
            // block_size=25 -> 33MB * 5 + hash = ~180MB
            return 200 * 1024 * 1024;

        case compression_backend_t::zstd:
            // ZSTD level 10: internal encoder state ~300MB
            return 350 * 1024 * 1024;

        case compression_backend_t::brotli:
            // Brotli: 16MB window + encoder state
            return 30 * 1024 * 1024;

        default:
            return 100 * 1024 * 1024;
    }
}

MemoryEstimate ResourceManager::EstimateBlockMemory(
    size_t n_samples,
    size_t ploidy,
    size_t max_block_rows,
    size_t n_col_blocks)
{
    MemoryEstimate est;

    // Vector length in bytes
    size_t n_haplotypes = n_samples * ploidy;
    size_t vec_len = (n_haplotypes + 7) / 8;

    // 1. GT data double buffer (for all column blocks)
    est.gt_data = vec_len * max_block_rows * 2 * n_col_blocks;

    // 2. variant_desc_t metadata (~300 bytes per variant)
    // Includes: chrom, pos, id, ref, alt, qual, filter, info, format strings
    est.variant_metadata = max_block_rows * 300;

    // 3. Fixed field blocks (7 fields * 8KB each)
    est.fixed_fields = 7 * 8 * 1024;

    // 4. Bit vectors (zeros_only, copies, unique + rank support)
    // Each bit vector: max_block_rows / 8 bytes
    // With rank support overhead: ~1.25x
    est.bit_vectors = (max_block_rows / 8) * 5 * 2;  // 5 vectors, 2x for rank

    // 5. Compression temporary buffers (approximately input size)
    est.compression_buffer = est.gt_data;

    return est;
}

void ResourceManager::CalculateMemoryBudget() {
    // Calculate per-block memory estimate first
    block_estimate_ = EstimateBlockMemory(
        workload_.n_samples,
        workload_.ploidy,
        workload_.max_block_rows,
        workload_.n_col_blocks);

    // Total budget: user limit or 70% of available memory
    if (user_max_memory_mb_ > 0) {
        memory_budget_.total_budget_mb = user_max_memory_mb_;
        memory_budget_.user_specified = true;
    } else {
        memory_budget_.total_budget_mb = sys_info_.available_memory_mb * 70 / 100;
        memory_budget_.user_specified = false;

        // Ensure a sane minimum for auto planning (not a hard limit).
        if (memory_budget_.total_budget_mb < 500) {
            memory_budget_.total_budget_mb = 500;
        }
    }

    // 1. Backend memory reserve (fixed based on backend type)
    memory_budget_.backend_reserve_mb =
        GetBackendMemoryReserve(workload_.backend) / (1024 * 1024);

    // 2. Fixed overhead (~100MB)
    memory_budget_.overhead_mb = 100;

    // 3. Queue budget (remaining)
    uint64_t reserved = memory_budget_.backend_reserve_mb + memory_budget_.overhead_mb;
    if (memory_budget_.total_budget_mb > reserved) {
        memory_budget_.queue_budget_mb = memory_budget_.total_budget_mb - reserved;
    } else {
        // Not enough memory; queue budget becomes 0 for planning purposes.
        memory_budget_.queue_budget_mb = 0;
    }

    // Store per-block estimate in MB
    memory_budget_.per_block_estimate_mb = block_estimate_.total_mb();
    if (memory_budget_.per_block_estimate_mb == 0) {
        memory_budget_.per_block_estimate_mb = 1;
    }
}

ThreadConfig ResourceManager::AutoConfigureThreads() {
    ThreadConfig cfg;

    // Available cores for work (reserve 1 for system)
    int available_cores = std::max(1, (int)sys_info_.cpu_cores - 1);

    // ========== GT Thread Count ==========
    // Considerations:
    // - Column blocks can be processed in parallel
    // - Memory limit per thread
    // - CPU cores available

    int gt_by_cols = std::max(1, (int)workload_.n_col_blocks);
    int gt_by_cores = std::max(1, available_cores / 2);
    int gt_by_memory = 1;

    if (memory_budget_.per_block_estimate_mb > 0 && memory_budget_.queue_budget_mb > 0) {
        // How many blocks can fit in queue memory?
        size_t blocks_fit = memory_budget_.queue_budget_mb / memory_budget_.per_block_estimate_mb;
        // Each GT thread needs ~2 blocks in flight
        gt_by_memory = std::max(1, (int)(blocks_fit / 2));
    }

    int auto_gt_threads = std::min({gt_by_cols, gt_by_cores, gt_by_memory});
    cfg.no_gt_threads = std::max(1u, std::min((uint32_t)auto_gt_threads, 16u));

    // ========== Compression Thread Count ==========
    // Use remaining cores for compression
    int remaining_cores = std::max(1, available_cores - (int)cfg.no_gt_threads);
    cfg.no_threads = std::max(1u, std::min((uint32_t)remaining_cores, 16u));

    // ========== Parse Thread Count ==========
    // Usually 1-2 is enough (I/O bound)
    cfg.no_parse_threads = std::min(2u, (uint32_t)(available_cores / 4 + 1));
    if (cfg.no_parse_threads == 0) cfg.no_parse_threads = 1;

    // ========== Queue Capacity ==========
    // Based on memory budget and block size
    size_t blocks_fit = 8;  // Default
    if (memory_budget_.per_block_estimate_mb > 0) {
        blocks_fit = memory_budget_.queue_budget_mb / memory_budget_.per_block_estimate_mb;
    }
    cfg.queue_capacity = std::max((size_t)4, std::min(blocks_fit, (size_t)64));

    // Ensure queue can hold at least gt_threads * 2 blocks
    size_t min_queue = (size_t)(cfg.no_gt_threads * 2);
    cfg.queue_capacity = std::max(cfg.queue_capacity, min_queue);

    return cfg;
}

ThreadConfig ResourceManager::Configure() {
    // Calculate memory budget first
    CalculateMemoryBudget();

    // Get auto-configured values
    ThreadConfig auto_cfg = AutoConfigureThreads();

    // Apply user overrides
    final_config_.no_threads = (user_no_threads_ > 0) ?
        user_no_threads_ : auto_cfg.no_threads;

    final_config_.no_gt_threads = (user_no_gt_threads_ > 0) ?
        user_no_gt_threads_ : auto_cfg.no_gt_threads;

    final_config_.no_parse_threads = (user_no_parse_threads_ > 0) ?
        user_no_parse_threads_ : auto_cfg.no_parse_threads;

    final_config_.queue_capacity = (user_queue_capacity_ > 0) ?
        user_queue_capacity_ : auto_cfg.queue_capacity;

    return final_config_;
}

void ResourceManager::LogResourcePlan() const {
    auto logger = LogManager::Instance().Logger();

    logger->info("=== Resource Configuration ===");
    logger->info("System: {} cores, {} MB available, {} MB total",
                 sys_info_.cpu_cores,
                 sys_info_.available_memory_mb,
                 sys_info_.total_memory_mb);

    // Backend name
    std::string backend_name;
    switch (workload_.backend) {
        case compression_backend_t::bsc: backend_name = "BSC"; break;
        case compression_backend_t::zstd: backend_name = "ZSTD"; break;
        case compression_backend_t::brotli: backend_name = "Brotli"; break;
        default: backend_name = "Unknown"; break;
    }
    logger->info("Backend: {} (reserve {} MB)", backend_name, memory_budget_.backend_reserve_mb);

    logger->info("Memory budget: {} MB{} (queue: {} MB, overhead: {} MB)",
                 memory_budget_.total_budget_mb,
                 memory_budget_.user_specified ? " (user specified)" : " (auto: 70% available)",
                 memory_budget_.queue_budget_mb,
                 memory_budget_.overhead_mb);

    logger->info("Workload: {} samples, ploidy={}, {} variants/block, {} col blocks",
                 workload_.n_samples,
                 workload_.ploidy,
                 workload_.max_block_rows,
                 workload_.n_col_blocks);

    logger->info("Per-block estimate: {} MB (GT:{} + meta:{} + fixed:{} + bits:{} + comp:{})",
                 block_estimate_.total_mb(),
                 block_estimate_.gt_data / (1024 * 1024),
                 block_estimate_.variant_metadata / (1024 * 1024),
                 block_estimate_.fixed_fields / (1024 * 1024),
                 block_estimate_.bit_vectors / (1024 * 1024),
                 block_estimate_.compression_buffer / (1024 * 1024));

    logger->info("Threads: compress={}{}, gt={}{}, parse={}{}",
                 final_config_.no_threads,
                 (user_no_threads_ > 0) ? " (user)" : " (auto)",
                 final_config_.no_gt_threads,
                 (user_no_gt_threads_ > 0) ? " (user)" : " (auto)",
                 final_config_.no_parse_threads,
                 (user_no_parse_threads_ > 0) ? " (user)" : " (auto)");

    logger->info("Queue capacity: {} blocks{}",
                 final_config_.queue_capacity,
                 (user_queue_capacity_ > 0) ? " (user)" : " (auto)");

    logger->info("==============================");
}

} // namespace gsc
