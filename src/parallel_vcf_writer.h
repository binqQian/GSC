#ifndef _PARALLEL_VCF_WRITER_H
#define _PARALLEL_VCF_WRITER_H

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <map>
#include <atomic>
#include <memory>
#include "logger.h"
#include "defs.h"

namespace gsc {

// Constants for parallel writing
const size_t VCF_WRITE_QUEUE_SIZE = 1000;  // Max queued records

/**
 * @brief Parallel VCF/BCF writer with dedicated write thread
 *
 * This class provides accelerated VCF/BCF writing by:
 * 1. Offloading bcf_write1() calls to a dedicated writer thread
 * 2. Using a thread-safe queue with ordered delivery
 * 3. Allowing main thread to continue processing while writes happen in background
 *
 * Note: HTSlib's bcf_write1() must be called from a single thread, but we can
 * parallelize the queueing and preparation of records.
 */
class ParallelVCFWriter {
private:
    // Write request with sequence number for ordering
    struct WriteRequest {
        uint64_t seq;
        bcf1_t* rec;
        bool is_sentinel;  // True for shutdown signal

        WriteRequest() : seq(0), rec(nullptr), is_sentinel(false) {}
        WriteRequest(uint64_t s, bcf1_t* r) : seq(s), rec(r), is_sentinel(false) {}
    };

    // Thread-safe priority queue for ordered writes
    class OrderedWriteQueue {
    private:
        std::map<uint64_t, bcf1_t*> pending_writes_;
        mutable std::mutex mutex_;
        std::condition_variable data_ready_;
        std::condition_variable space_available_;
        uint64_t next_seq_to_write_;
        size_t max_queue_size_;
        bool finished_;
        std::atomic<size_t> current_size_;

    public:
        explicit OrderedWriteQueue(size_t max_size)
            : next_seq_to_write_(0)
            , max_queue_size_(max_size)
            , finished_(false)
            , current_size_(0) {}

        // Submit a record for writing (blocks if queue is full)
        void push(uint64_t seq, bcf1_t* rec) {
            std::unique_lock<std::mutex> lock(mutex_);
            space_available_.wait(lock, [this]() {
                return current_size_.load() < max_queue_size_ || finished_;
            });

            if (!finished_) {
                pending_writes_[seq] = rec;
                current_size_++;
                lock.unlock();
                data_ready_.notify_one();
            }
            else if (rec) {
                bcf_destroy(rec);  // Cleanup if shutting down
            }
        }

        // Get next record in sequence (blocks if not ready)
        bool pop(bcf1_t*& rec) {
            std::unique_lock<std::mutex> lock(mutex_);
            data_ready_.wait(lock, [this]() {
                return pending_writes_.count(next_seq_to_write_) > 0 || finished_;
            });

            if (pending_writes_.count(next_seq_to_write_) > 0) {
                rec = pending_writes_[next_seq_to_write_];
                pending_writes_.erase(next_seq_to_write_);
                next_seq_to_write_++;
                current_size_--;
                lock.unlock();
                space_available_.notify_one();
                return true;
            }

            return false;  // Finished and no more data
        }

        void finish() {
            std::unique_lock<std::mutex> lock(mutex_);
            finished_ = true;
            lock.unlock();
            data_ready_.notify_all();
            space_available_.notify_all();
        }

        size_t pending_count() const {
            return current_size_.load();
        }
    };

    // Member variables
    htsFile* fp_out_;
    bcf_hdr_t* hdr_;
    bool initialized_;
    bool use_parallel_writing_;

    // Writing infrastructure
    OrderedWriteQueue write_queue_;
    std::thread writer_thread_;
    std::atomic<uint64_t> next_write_seq_;
    std::atomic<uint64_t> total_written_;

    // Writer thread function
    void WriterThread();

public:
    /**
     * @brief Construct parallel VCF writer
     */
    ParallelVCFWriter();

    ~ParallelVCFWriter();

    /**
     * @brief Initialize writer and create output file
     * @param filename Output file path (VCF/BCF/VCF.gz based on extension)
     * @param hdr VCF/BCF header
     * @param mode Write mode ("w" for VCF, "wb" for BCF, "wz" for VCF.gz)
     * @param use_parallel Enable parallel writing (default: true)
     * @return true on success, false on failure
     */
    bool Initialize(const char* filename, bcf_hdr_t* hdr, const char* mode = "wz",
                    bool use_parallel = true);

    /**
     * @brief Initialize writer with already-opened file handle (no header write)
     * @param fp Already-opened htsFile* (takes ownership, will close on Finalize)
     * @param hdr VCF/BCF header (must be already written to fp)
     * @param use_parallel Enable parallel writing (default: true)
     * @return true on success, false on failure
     */
    bool InitializeWithHandle(htsFile* fp, bcf_hdr_t* hdr, bool use_parallel = true);

    /**
     * @brief Write a variant record (takes ownership of rec)
     * @param rec bcf1_t* record to write (will be freed by writer)
     * @return true on success, false on failure
     */
    bool WriteRecord(bcf1_t* rec);

    /**
     * @brief Get total number of written variants
     * @return Variant count
     */
    uint64_t GetTotalWritten() const { return total_written_.load(); }

    /**
     * @brief Get number of pending writes in queue
     * @return Pending count
     */
    size_t GetPendingCount() const { return write_queue_.pending_count(); }

    /**
     * @brief Flush all pending writes and close file
     */
    void Finalize();
};

} // namespace gsc

#endif // _PARALLEL_VCF_WRITER_H
