#pragma once

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

// 常量定义
const unsigned long long BATCH_SIZE = 50000;      // 后续根据样本数自适应修改
const int BUFFER_SIZE = 12;

// VCF/BCF解析器类
class vcf_Parser {
private:
    // 结构体定义
    struct VcfLine {
        unsigned long long seq;
        std::string line;
        bool is_valid;
        
        VcfLine() : seq(0), is_valid(false) {}
    };

    // 线程安全队列模板
    template<typename T>
    class ThreadSafeQueue {
    private:
        std::queue<T> queue;
        mutable std::mutex mutex;
        std::condition_variable not_empty;
        std::condition_variable not_full;
        size_t capacity;
        bool finished;

    public:
        ThreadSafeQueue(size_t max_size) : capacity(max_size), finished(false) {}

        void push(const T& value) {
            std::unique_lock<std::mutex> lock(mutex);
            not_full.wait(lock, [this]() { return queue.size() < capacity || finished; });
            if (!finished) {
                queue.push(value);
                lock.unlock();
                not_empty.notify_one();
            }
        }

        bool pop(T& value) {
            std::unique_lock<std::mutex> lock(mutex);
            not_empty.wait(lock, [this]() { return !queue.empty() || finished; });
            if (!queue.empty()) {
                value = std::move(queue.front());
                queue.pop();
                lock.unlock();
                not_full.notify_one();
                return true;
            }
            return false;
        }

        void finish() {
            std::unique_lock<std::mutex> lock(mutex);
            finished = true;
            lock.unlock();
            not_empty.notify_all();
            not_full.notify_all();
        }

        bool is_finished() const {
            std::lock_guard<std::mutex> lock(mutex);
            return finished && queue.empty();
        }
    };

    // 解析线程的上下文结构
    struct ParserContext {
        bcf_hdr_t* hdr;
        bcf1_t* rec;
        kstring_t kstr;
        
        ParserContext(bcf_hdr_t* header) {
            hdr = bcf_hdr_dup(header);
            rec = bcf_init();
            kstr = {0, 0, NULL};
        }
        
        ~ParserContext() {
            if (kstr.s) free(kstr.s);
            bcf_destroy(rec);
            bcf_hdr_destroy(hdr);
        }
        
        ParserContext(const ParserContext&) = delete;
        ParserContext& operator=(const ParserContext&) = delete;
    };

    // 私有成员变量
    std::mutex results_mutex;
    std::condition_variable results_cv;
    std::map<unsigned long long, bcf1_t*> parsed_results;
    unsigned long long next_seq_to_process = 0;
    bool parsing_finished = false;
    std::atomic<unsigned long long> total_parsed{0};
    htsFile* fp_hts;
    bcf_hdr_t* hdr;
    ThreadSafeQueue<VcfLine> line_queue;
    std::vector<std::thread> parser_threads;
    int num_threads;
    bool initialized = false;
    bool use_parallel_parsing = false;

    // 解析线程函数
    void parser_thread(ThreadSafeQueue<VcfLine>& input_queue, bcf_hdr_t* header) {
        std::unique_ptr<ParserContext> ctx(new ParserContext(header));
        VcfLine line;
        
        while (input_queue.pop(line)) {
            if (!line.is_valid) continue;
            
            ctx->kstr.l = line.line.size();
            if (ctx->kstr.m < line.line.size() + 1) {
                ctx->kstr.s = (char*)realloc(ctx->kstr.s, line.line.size() + 1);
                ctx->kstr.m = line.line.size() + 1;
            }
            memcpy(ctx->kstr.s, line.line.c_str(), line.line.size());
            ctx->kstr.s[line.line.size()] = '\0';

            if (vcf_parse1(&ctx->kstr, ctx->hdr, ctx->rec) >= 0) {
                bcf1_t* rec_copy = bcf_dup(ctx->rec);
                if (rec_copy) {
                    total_parsed++;
                    std::lock_guard<std::mutex> lock(results_mutex);
                    parsed_results[line.seq] = rec_copy;
                    results_cv.notify_one();
                }
            }
        }
    }

public:
    vcf_Parser(int threads = 1) : line_queue(BUFFER_SIZE), num_threads(threads) {}
    
    ~vcf_Parser() {
        cleanup();
    }

    // 初始化解析器
    bool initialize(const char* filename) {
        if (initialized) return false;
        
        // 使用 hts_open 自动检测并打开不同类型的文件
        fp_hts = hts_open(filename, "r");
        if (!fp_hts) return false;

        // 读取 VCF/BCF 头信息
        hdr = bcf_hdr_read(fp_hts);
        if (!hdr) {
            cleanup();
            return false;
        }

        // 根据文件类型决定解析策略
        if (fp_hts->format.format == htsExactFormat::bcf) {
            // BCF 文件不需要文本解析加速
            use_parallel_parsing = false;
        } else if (fp_hts->format.format == htsExactFormat::vcf) {
            // VCF 文件使用并行解析
            use_parallel_parsing = true;
            
            // 如果是压缩 VCF，设置 BGZF 多线程解压
            if (fp_hts->format.compression == htsCompression::bgzf) {
                if (bgzf_mt(fp_hts->fp.bgzf, num_threads, 256) != 0) {
                    cleanup();
                    return false;
                }
            }
        } else {
            // 不支持的文件格式
            cleanup();
            return false;
        }

        // 仅在需要时启动解析线程
        if (use_parallel_parsing) {
            for (int i = 0; i < num_threads; ++i) {
                parser_threads.emplace_back(&vcf_Parser::parser_thread, this, 
                                          std::ref(line_queue), hdr);
            }
        }

        initialized = true;
        return true;
    }

    // 获取VCF/BCF header
    bcf_hdr_t* get_header() const {
        return hdr;
    }

    // 解析下一批变体
    std::vector<bcf1_t*> parse_next_batch() {
        if (!initialized) return std::vector<bcf1_t*>();
        
        std::vector<bcf1_t*> batch_records;
        
        if (!use_parallel_parsing) {
            // 对 BCF 或其他不需要并行解析的格式，直接使用 htslib 函数
            for (unsigned long long i = 0; i < BATCH_SIZE; ++i) {
                bcf1_t* rec = bcf_init();
                if (bcf_read(fp_hts, hdr, rec) == 0) {
                    batch_records.push_back(rec);
                    total_parsed++;
                } else {
                    bcf_destroy(rec);
                    break;
                }
            }
            return batch_records;
        }
        
        // 以下是并行解析 VCF 的代码
        unsigned long long batch_count = 0;
        kstring_t kstr = {0, 0, NULL};
        
        // 读取并提交新的变体进行解析
        while (batch_count < BATCH_SIZE) {
            // 使用 htslib 函数读取文本行
            int ret;
            if (fp_hts->format.compression == htsCompression::bgzf) {
                ret = bgzf_getline(fp_hts->fp.bgzf, '\n', &kstr);
            } else {
                ret = hts_getline(fp_hts, KS_SEP_LINE, &kstr);
            }
            
            if (ret < 0) break; // 文件结束
            
            if (kstr.l == 0 || kstr.s[0] == '#') continue;

            VcfLine vcf_line;
            vcf_line.seq = next_seq_to_process + batch_count;
            vcf_line.line = std::string(kstr.s, kstr.l);
            vcf_line.is_valid = true;
            
            line_queue.push(vcf_line);
            batch_count++;
        }

        // 如果文件已读完且没有新的变体
        if (batch_count == 0) {
            if (kstr.s) free(kstr.s);
            return batch_records;
        }

        // 收集解析结果
        for (unsigned long long i = 0; i < batch_count; ++i) {
            std::unique_lock<std::mutex> lock(results_mutex);
            results_cv.wait(lock, [this, i]() {
                return parsed_results.find(next_seq_to_process + i) != parsed_results.end();
            });

            batch_records.push_back(parsed_results[next_seq_to_process + i]);
            parsed_results.erase(next_seq_to_process + i);
        }

        next_seq_to_process += batch_count;
        if (kstr.s) free(kstr.s);
        return batch_records;
    }

    // 获取已解析变体数量
    unsigned long long get_total_parsed() const {
        return total_parsed;
    }

    // 清理资源
    void cleanup() {
        if (!initialized) return;

        if (use_parallel_parsing) {
            line_queue.finish();
            for (auto& thread : parser_threads) {
                if (thread.joinable()) thread.join();
            }
        }

        // 清理未处理的结果
        for (auto& pair : parsed_results) {
            bcf_destroy(pair.second);
        }
        parsed_results.clear();

        if (hdr) bcf_hdr_destroy(hdr);
        if (fp_hts) hts_close(fp_hts);

        initialized = false;
    }
};