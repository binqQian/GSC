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
};

void InitializeCompressionBackend(compression_backend_t backend);
std::unique_ptr<CompressionStrategy> MakeCompressionStrategy(compression_backend_t backend, const bsc_params_t& bsc_params);
