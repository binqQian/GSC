#ifndef _PARALLEL_VCF_READER_H
#define _PARALLEL_VCF_READER_H

#include <htslib/bgzf.h>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <memory>
#include <map>
#include <string>
#include "logger.h"
#include "defs.h"

namespace gsc {

// Constants for parallel parsing
const size_t VCF_PARSE_BATCH_SIZE = 50000;  // Variants per batch
const size_t VCF_PARSE_QUEUE_SIZE = 12;     // Line queue capacity

/**
 * @brief Parallel VCF/BCF reader with multi-threaded parsing
 *
 * This class provides accelerated VCF/BCF reading through:
 * 1. BGZF multi-threaded decompression for .vcf.gz files
 * 2. Parallel text parsing using multiple worker threads
 * 3. Ordered result delivery to maintain variant sequence
 *
 * For BCF files, falls back to single-threaded bcf_read() since
 * binary parsing is already fast and doesn't benefit from parallelization.
 */
class ParallelVCFReader {
private:
    // VCF text line with sequence number for ordered delivery
    struct VCFLine {
        uint64_t seq;
        std::string line;
        bool is_valid;

        VCFLine() : seq(0), is_valid(false) {}
    };

    // Thread-safe queue template with capacity limit
    template<typename T>
    class ThreadSafeQueue {
    private:
        std::queue<T> queue_;
        mutable std::mutex mutex_;
        std::condition_variable not_empty_;
        std::condition_variable not_full_;
        size_t capacity_;
        bool finished_;

    public:
        explicit ThreadSafeQueue(size_t max_size)
            : capacity_(max_size), finished_(false) {}

        void push(const T& value) {
            std::unique_lock<std::mutex> lock(mutex_);
            not_full_.wait(lock, [this]() {
                return queue_.size() < capacity_ || finished_;
            });
            if (!finished_) {
                queue_.push(value);
                lock.unlock();
                not_empty_.notify_one();
            }
        }

        bool pop(T& value) {
            std::unique_lock<std::mutex> lock(mutex_);
            not_empty_.wait(lock, [this]() {
                return !queue_.empty() || finished_;
            });
            if (!queue_.empty()) {
                value = std::move(queue_.front());
                queue_.pop();
                lock.unlock();
                not_full_.notify_one();
                return true;
            }
            return false;
        }

        void finish() {
            std::unique_lock<std::mutex> lock(mutex_);
            finished_ = true;
            lock.unlock();
            not_empty_.notify_all();
            not_full_.notify_all();
        }

        bool is_finished() const {
            std::lock_guard<std::mutex> lock(mutex_);
            return finished_ && queue_.empty();
        }
    };

    // Parser thread context with independent HTSlib structures
    struct ParserContext {
        bcf_hdr_t* hdr;
        bcf1_t* rec;
        kstring_t kstr;

        explicit ParserContext(bcf_hdr_t* header) {
            hdr = bcf_hdr_dup(header);
            rec = bcf_init();
            kstr = {0, 0, nullptr};
        }

        ~ParserContext() {
            if (kstr.s) free(kstr.s);
            bcf_destroy(rec);
            bcf_hdr_destroy(hdr);
        }

        // Prevent copy
        ParserContext(const ParserContext&) = delete;
        ParserContext& operator=(const ParserContext&) = delete;
    };

    // Member variables
    htsFile* fp_hts_;
    bcf_hdr_t* hdr_;
    int num_threads_;
    bool initialized_;
    bool use_parallel_parsing_;

    // Parallel parsing infrastructure
    ThreadSafeQueue<VCFLine> line_queue_;
    std::vector<std::thread> parser_threads_;
    std::mutex results_mutex_;
    std::condition_variable results_cv_;
    std::map<uint64_t, bcf1_t*> parsed_results_;
    uint64_t next_seq_to_process_;
    std::atomic<uint64_t> total_parsed_;
    bool parsing_finished_;

    // Parser worker thread function
    void ParserThread();

public:
    /**
     * @brief Construct parallel VCF reader
     * @param threads Number of parsing threads (ignored for BCF)
     */
    explicit ParallelVCFReader(int threads = 1);

    ~ParallelVCFReader();

    /**
     * @brief Initialize reader and open file
     * @param filename Path to VCF/BCF/VCF.gz file
     * @return true on success, false on failure
     */
    bool Initialize(const char* filename);

    /**
     * @brief Get VCF/BCF header
     * @return Pointer to bcf_hdr_t (owned by reader)
     */
    bcf_hdr_t* GetHeader() const { return hdr_; }

    /**
     * @brief Parse next batch of variants
     * @return Vector of bcf1_t* records (caller must free with bcf_destroy)
     */
    std::vector<bcf1_t*> ParseNextBatch();

    /**
     * @brief Get total number of parsed variants
     * @return Variant count
     */
    uint64_t GetTotalParsed() const { return total_parsed_.load(); }

    /**
     * @brief Check if parallel parsing is enabled
     * @return true for VCF.gz, false for BCF
     */
    bool IsParallelMode() const { return use_parallel_parsing_; }

    /**
     * @brief Cleanup resources and close file
     */
    void Cleanup();
};

} // namespace gsc

#endif // _PARALLEL_VCF_READER_H
