#pragma once

#include <memory>
#include <vector>

#include "bsc.h"
#include "gsc_params.h"
#include "logger.h"
#include "zstd_compress.h"

// 后端压缩抽象层：
// 主流程只依赖这组接口，不直接关心底层是 bsc / zstd / brotli。
class CompressionStrategy {
public:
    virtual ~CompressionStrategy() = default;
    virtual bool Compress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) = 0;
    virtual bool Decompress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) = 0;
    // 允许直接从 mmap / 文件缓冲指针解压，避免先构造一个临时 vector。
    virtual bool DecompressFromPtr(const uint8_t* src, size_t src_size, std::vector<uint8_t>& output) {
        // 默认实现会先复制一份；真正高频路径可在具体后端里覆写。
        std::vector<uint8_t> tmp(src, src + src_size);
        return Decompress(tmp, output);
    }
};

// 初始化底层后端库的全局状态。
void InitializeCompressionBackend(compression_backend_t backend);
// 根据参数创建一个具体策略实例。
std::unique_ptr<CompressionStrategy> MakeCompressionStrategy(compression_backend_t backend, const bsc_params_t& bsc_params);
