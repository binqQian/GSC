#pragma once

#include <memory>
#include <vector>

#include "bsc.h"
#include "gsc_params.h"
#include "logger.h"
#include "zstd_compress.h"

// Abstract strategy for compress/decompress byte buffers.
class CompressionStrategy {
public:
    virtual ~CompressionStrategy() = default;
    virtual bool Compress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) = 0;
    virtual bool Decompress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) = 0;
    // Decompress directly from pointer, avoiding intermediate copy
    virtual bool DecompressFromPtr(const uint8_t* src, size_t src_size, std::vector<uint8_t>& output) {
        // Default implementation: copy to vector and call Decompress
        std::vector<uint8_t> tmp(src, src + src_size);
        return Decompress(tmp, output);
    }
};

void InitializeCompressionBackend(compression_backend_t backend);
std::unique_ptr<CompressionStrategy> MakeCompressionStrategy(compression_backend_t backend, const bsc_params_t& bsc_params);
