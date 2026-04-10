#include "parallel_vcf_writer.h"
#include "logger.h"

namespace gsc {

// writer 默认是关闭状态；是否真的并行写由 Initialize() 决定。
ParallelVCFWriter::ParallelVCFWriter()
    : fp_out_(nullptr)
    , hdr_(nullptr)
    , initialized_(false)
    , use_parallel_writing_(false)
    , write_queue_(VCF_WRITE_QUEUE_SIZE)
    , next_write_seq_(0)
    , total_written_(0)
    , max_pool_size_(VCF_WRITE_QUEUE_SIZE * 2)
{
}

ParallelVCFWriter::~ParallelVCFWriter()
{
    Finalize();
}

bool ParallelVCFWriter::Initialize(const char* filename, bcf_hdr_t* hdr,
                                    const char* mode, bool use_parallel)
{
    if (initialized_) {
        LogManager::Instance().Logger()->warn("ParallelVCFWriter already initialized");
        return false;
    }

    if (!hdr) {
        LogManager::Instance().Logger()->error("Invalid header for ParallelVCFWriter");
        return false;
    }

    // 打开输出文件并写 header；真正的数据写出可以放到后台线程。
    fp_out_ = hts_open(filename, mode);
    if (!fp_out_) {
        LogManager::Instance().Logger()->error("Failed to open output file: {}", filename);
        return false;
    }

    // Duplicate header for thread safety
    hdr_ = bcf_hdr_dup(hdr);
    if (!hdr_) {
        LogManager::Instance().Logger()->error("Failed to duplicate VCF header");
        hts_close(fp_out_);
        fp_out_ = nullptr;
        return false;
    }

    // Write header
    if (bcf_hdr_write(fp_out_, hdr_) != 0) {
        LogManager::Instance().Logger()->error("Failed to write VCF header");
        bcf_hdr_destroy(hdr_);
        hdr_ = nullptr;
        hts_close(fp_out_);
        fp_out_ = nullptr;
        return false;
    }

    use_parallel_writing_ = use_parallel;

    // 并行模式下，writer 线程独占所有 htslib 写调用。
    if (use_parallel_writing_) {
        writer_thread_ = std::thread(&ParallelVCFWriter::WriterThread, this);
        LogManager::Instance().Logger()->debug("Parallel VCF writer initialized");
    }
    else {
        LogManager::Instance().Logger()->debug("Serial VCF writer initialized");
    }

    initialized_ = true;
    return true;
}

bool ParallelVCFWriter::InitializeWithHandle(htsFile* fp, bcf_hdr_t* hdr, bool use_parallel)
{
    if (initialized_) {
        LogManager::Instance().Logger()->warn("ParallelVCFWriter already initialized");
        return false;
    }

    if (!fp) {
        LogManager::Instance().Logger()->error("Invalid file handle for ParallelVCFWriter");
        return false;
    }

    if (!hdr) {
        LogManager::Instance().Logger()->error("Invalid header for ParallelVCFWriter");
        return false;
    }

    // 这里直接接管已打开的句柄；header 视为外部已经写好。
    fp_out_ = fp;

    // Duplicate header for thread safety
    hdr_ = bcf_hdr_dup(hdr);
    if (!hdr_) {
        LogManager::Instance().Logger()->error("Failed to duplicate VCF header");
        fp_out_ = nullptr;  // Don't close since we don't own it yet
        return false;
    }

    use_parallel_writing_ = use_parallel;

    // Launch writer thread if parallel writing is enabled
    if (use_parallel_writing_) {
        writer_thread_ = std::thread(&ParallelVCFWriter::WriterThread, this);
        LogManager::Instance().Logger()->debug("Parallel VCF writer initialized with existing handle");
    }
    else {
        LogManager::Instance().Logger()->debug("Serial VCF writer initialized with existing handle");
    }

    initialized_ = true;
    return true;
}

void ParallelVCFWriter::WriterThread()
{
    // 唯一真正调用 bcf_write1() 的线程。
    // 主线程只负责排队，不直接碰 htslib 写接口。
    bcf1_t* rec = nullptr;

    while (write_queue_.pop(rec)) {
        if (!rec) continue;

        // 串行写出一条记录；失败只记日志，不中断整个 writer。
        if (bcf_write1(fp_out_, hdr_, rec) < 0) {
            LogManager::Instance().Logger()->error(
                "Failed to write variant at position {}",
                rec->pos);
        }
        else {
            total_written_++;
        }

        // 写完后把 record 回收到对象池，降低分配成本。
        ReleaseRecord(rec);

        // Log progress periodically
        if (total_written_ % 100000 == 0) {
            LogManager::Instance().Logger()->debug(
                "Written {} variants, {} pending",
                total_written_.load(), write_queue_.pending_count());
        }
    }

    LogManager::Instance().Logger()->debug("Writer thread finished, total written: {}",
                                           total_written_.load());
}

bcf1_t* ParallelVCFWriter::AcquireRecord()
{
    // 优先复用对象池里的 bcf1_t；没有时再新建。
    std::lock_guard<std::mutex> lock(pool_mutex_);
    if (!record_pool_.empty())
    {
        bcf1_t* rec = record_pool_.back();
        record_pool_.pop_back();
        bcf_clear(rec);
        return rec;
    }
    return bcf_init();
}

void ParallelVCFWriter::ReleaseRecord(bcf1_t* rec)
{
    // 池满了就直接销毁，避免对象池无限增长。
    if (!rec)
        return;
    std::lock_guard<std::mutex> lock(pool_mutex_);
    if (record_pool_.size() >= max_pool_size_)
    {
        bcf_destroy(rec);
        return;
    }
    bcf_clear(rec);
    record_pool_.push_back(rec);
}

bool ParallelVCFWriter::WriteRecord(bcf1_t* rec)
{
    if (!initialized_) {
        LogManager::Instance().Logger()->error("ParallelVCFWriter not initialized");
        if (rec) bcf_destroy(rec);
        return false;
    }

    if (!rec) {
        LogManager::Instance().Logger()->warn("Attempted to write null record");
        return false;
    }

    if (use_parallel_writing_) {
        // 并行模式下只负责分配 seq 并投递给 writer 线程。
        uint64_t seq = next_write_seq_.fetch_add(1);
        write_queue_.push(seq, rec);
    }
    else {
        // 串行模式直接在当前线程写出。
        if (bcf_write1(fp_out_, hdr_, rec) < 0) {
            LogManager::Instance().Logger()->error(
                "Failed to write variant at position {}", rec->pos);
            ReleaseRecord(rec);
            return false;
        }
        total_written_++;
        ReleaseRecord(rec);
    }

    return true;
}

void ParallelVCFWriter::Finalize()
{
    if (!initialized_) return;

    LogManager::Instance().Logger()->info("Finalizing VCF writer, total written: {}",
                                          total_written_.load());

    // 先停 writer，再关文件，最后清对象池。
    if (use_parallel_writing_ && writer_thread_.joinable()) {
        write_queue_.finish();
        writer_thread_.join();
        LogManager::Instance().Logger()->debug("Writer thread joined");
    }

    // Close file
    if (fp_out_) {
        hts_close(fp_out_);
        fp_out_ = nullptr;
    }

    // Clean up header
    if (hdr_) {
        bcf_hdr_destroy(hdr_);
        hdr_ = nullptr;
    }

    // Clean up pooled records
    {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        for (auto *rec : record_pool_)
            bcf_destroy(rec);
        record_pool_.clear();
    }

    initialized_ = false;
    LogManager::Instance().Logger()->info(
        "ParallelVCFWriter finalized, total records written: {}",
        total_written_.load());
}

} // namespace gsc
