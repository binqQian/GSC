#include "parallel_vcf_reader.h"
#include "logger.h"
#include <algorithm>

namespace gsc {

// 构造后先处于未初始化状态，真正打开文件在 Initialize() 里完成。
ParallelVCFReader::ParallelVCFReader(int threads)
    : fp_hts_(nullptr)
    , hdr_(nullptr)
    , num_threads_(threads)
    , initialized_(false)
    , use_parallel_parsing_(false)
    , line_queue_(VCF_PARSE_QUEUE_SIZE)
    , next_seq_to_process_(0)
    , total_parsed_(0)
    , parsing_finished_(false)
{
}

ParallelVCFReader::~ParallelVCFReader()
{
    Cleanup();
}

bool ParallelVCFReader::Initialize(const char* filename)
{
    if (initialized_) {
        LogManager::Instance().Logger()->warn("ParallelVCFReader already initialized");
        return false;
    }

    // 统一入口：hts_open 会自动识别 VCF / BCF / 压缩格式。
    fp_hts_ = hts_open(filename, "r");
    if (!fp_hts_) {
        LogManager::Instance().Logger()->error("Failed to open file: {}", filename);
        return false;
    }

    // Read VCF/BCF header
    hdr_ = bcf_hdr_read(fp_hts_);
    if (!hdr_) {
        LogManager::Instance().Logger()->error("Failed to read VCF/BCF header from: {}", filename);
        hts_close(fp_hts_);
        fp_hts_ = nullptr;
        return false;
    }

    // 解析策略选择：
    // - BCF 直接单线程 bcf_read()
    // - VCF/VCF.gz 走“主线程读文本 + worker 线程解析”
    if (fp_hts_->format.format == htsExactFormat::bcf) {
        // BCF files: use single-threaded bcf_read() (binary parsing is fast)
        use_parallel_parsing_ = false;
        LogManager::Instance().Logger()->debug("BCF format detected, using single-threaded parsing");
    }
    else if (fp_hts_->format.format == htsExactFormat::vcf) {
        // VCF files: enable parallel text parsing
        use_parallel_parsing_ = true;

        // Enable BGZF multi-threaded decompression for .vcf.gz
        if (fp_hts_->format.compression == htsCompression::bgzf) {
            if (bgzf_mt(fp_hts_->fp.bgzf, num_threads_, 256) != 0) {
                LogManager::Instance().Logger()->error("Failed to enable BGZF multi-threading");
                Cleanup();
                return false;
            }
            LogManager::Instance().Logger()->debug(
                "VCF.gz format detected, using {} threads for BGZF decompression + parallel parsing",
                num_threads_);
        }
        else {
            LogManager::Instance().Logger()->debug(
                "Uncompressed VCF detected, using {} threads for parallel parsing",
                num_threads_);
        }
    }
    else {
        LogManager::Instance().Logger()->error("Unsupported file format");
        Cleanup();
        return false;
    }

    // 只有文本 VCF 才需要启动解析线程。
    if (use_parallel_parsing_) {
        parser_threads_.reserve(num_threads_);
        for (int i = 0; i < num_threads_; ++i) {
            parser_threads_.emplace_back(&ParallelVCFReader::ParserThread, this);
        }
        LogManager::Instance().Logger()->debug("Spawned {} parser threads", num_threads_);
    }

    initialized_ = true;
    return true;
}

void ParallelVCFReader::ParserThread()
{
    // 每个线程独占一套 htslib 上下文，避免跨线程共享 rec/hdr。
    std::unique_ptr<ParserContext> ctx(new ParserContext(hdr_));
    VCFLine line;

    while (line_queue_.pop(line)) {
        if (!line.is_valid) continue;

        // 先把 std::string 行内容拷回 htslib 需要的 kstring。
        ctx->kstr.l = line.line.size();
        if (ctx->kstr.m < line.line.size() + 1) {
            ctx->kstr.s = (char*)realloc(ctx->kstr.s, line.line.size() + 1);
            ctx->kstr.m = line.line.size() + 1;
        }
        memcpy(ctx->kstr.s, line.line.c_str(), line.line.size());
        ctx->kstr.s[line.line.size()] = '\0';

        // 把一行 VCF 文本解析成 bcf1_t，再复制出一个可交给主线程的对象。
        bcf1_t* rec_copy = nullptr;
        if (vcf_parse1(&ctx->kstr, ctx->hdr, ctx->rec) >= 0) {
            rec_copy = bcf_dup(ctx->rec);
            if (rec_copy) {
                total_parsed_++;
            }
            else {
                LogManager::Instance().Logger()->error("Failed to duplicate bcf1_t record at sequence {}", line.seq);
            }
        }
        else {
            LogManager::Instance().Logger()->warn(
                "Failed to parse VCF line at sequence {}: {}",
                line.seq, line.line.substr(0, 50));
        }

        // 关键点：即使解析失败也要插入 nullptr 结果。
        // 否则主线程会一直等这个 seq，形成死锁。
        {
            std::lock_guard<std::mutex> lock(results_mutex_);
            parsed_results_[line.seq] = rec_copy;
            results_cv_.notify_one();
        }
    }
}

std::vector<bcf1_t*> ParallelVCFReader::ParseNextBatch()
{
    if (!initialized_) {
        LogManager::Instance().Logger()->error("ParallelVCFReader not initialized");
        return std::vector<bcf1_t*>();
    }

    std::vector<bcf1_t*> batch_records;
    batch_records.reserve(VCF_PARSE_BATCH_SIZE);

    // BCF 直接顺序读即可，省掉线程同步开销。
    if (!use_parallel_parsing_) {
        for (size_t i = 0; i < VCF_PARSE_BATCH_SIZE; ++i) {
            bcf1_t* rec = bcf_init();
            if (bcf_read(fp_hts_, hdr_, rec) == 0) {
                batch_records.push_back(rec);
                total_parsed_++;
            }
            else {
                bcf_destroy(rec);
                break; // EOF or error
            }
        }
        return batch_records;
    }

    // VCF 模式：
    // 1. 主线程顺序读文本行
    // 2. worker 并行解析
    // 3. 主线程按 seq 顺序回收结果
    uint64_t batch_count = 0;
    kstring_t kstr = {0, 0, nullptr};

    // 先把本批文本行投递出去。
    while (batch_count < VCF_PARSE_BATCH_SIZE) {
        int ret;
        if (fp_hts_->format.compression == htsCompression::bgzf) {
            ret = bgzf_getline(fp_hts_->fp.bgzf, '\n', &kstr);
        }
        else {
            ret = hts_getline(fp_hts_, KS_SEP_LINE, &kstr);
        }

        if (ret < 0) break; // EOF

        // Skip empty lines and header lines
        if (kstr.l == 0 || kstr.s[0] == '#') continue;

        VCFLine vcf_line;
        vcf_line.seq = next_seq_to_process_ + batch_count;
        vcf_line.line = std::string(kstr.s, kstr.l);
        vcf_line.is_valid = true;

        line_queue_.push(vcf_line);
        batch_count++;
    }

    // No more lines read
    if (batch_count == 0) {
        if (kstr.s) free(kstr.s);
        return batch_records;
    }

    // 再严格按 seq 顺序收集本批解析结果，保证变体顺序不变。
    for (uint64_t i = 0; i < batch_count; ++i) {
        uint64_t target_seq = next_seq_to_process_ + i;

        bcf1_t* rec = nullptr;
        {
            std::unique_lock<std::mutex> lock(results_mutex_);
            results_cv_.wait(lock, [this, target_seq]() {
                return parsed_results_.find(target_seq) != parsed_results_.end();
            });

            rec = parsed_results_[target_seq];
            parsed_results_.erase(target_seq);
        }

        // Only add non-null records (skip parse failures)
        if (rec != nullptr) {
            batch_records.push_back(rec);
        }
        else {
            LogManager::Instance().Logger()->debug(
                "Skipping failed variant at sequence {}", target_seq);
        }
    }

    next_seq_to_process_ += batch_count;
    if (kstr.s) free(kstr.s);

    return batch_records;
}

void ParallelVCFReader::Cleanup()
{
    if (!initialized_) return;

    // 先通知 worker 停止，再 join，最后清理残留结果和 htslib 资源。
    if (use_parallel_parsing_) {
        line_queue_.finish();
        for (auto& thread : parser_threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        parser_threads_.clear();
        LogManager::Instance().Logger()->debug("Parser threads joined");
    }

    // Clean up unparsed results
    for (auto& pair : parsed_results_) {
        bcf_destroy(pair.second);
    }
    parsed_results_.clear();

    // Close HTSlib resources
    if (hdr_) {
        bcf_hdr_destroy(hdr_);
        hdr_ = nullptr;
    }
    if (fp_hts_) {
        hts_close(fp_hts_);
        fp_hts_ = nullptr;
    }

    initialized_ = false;
    LogManager::Instance().Logger()->info(
        "ParallelVCFReader cleanup complete, total parsed: {}", total_parsed_.load());
}

} // namespace gsc
