#include "compression_strategy.h"

#include <brotli/decode.h>
#include <brotli/encode.h>
#include <algorithm>
#include <mutex>

namespace {

std::once_flag g_bsc_init_flag;

class BscCompressionStrategy : public CompressionStrategy {
public:
    explicit BscCompressionStrategy(const bsc_params_t& params) {
        wrapper_.InitCompress(params);
    }

    bool Compress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) override {
        return wrapper_.Compress(input, output);
    }

    bool Decompress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) override {
        std::vector<uint8_t> tmp = input;
        return wrapper_.Decompress(tmp, output);
    }

private:
    CBSCWrapper wrapper_;
};

class ZstdCompressionStrategy : public CompressionStrategy {
public:
    bool Compress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) override {
        return zstd::zstd_compress(input, output);
    }

    bool Decompress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) override {
        return zstd::zstd_decompress(input, output);
    }
};

class BrotliCompressionStrategy : public CompressionStrategy {
public:
    bool Compress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) override {
        size_t max_out = BrotliEncoderMaxCompressedSize(input.size());
        if (max_out == 0) return false;
        output.resize(max_out);
        size_t encoded = max_out;
        const bool ok = BrotliEncoderCompress(
            BROTLI_DEFAULT_QUALITY,
            BROTLI_DEFAULT_WINDOW,
            BROTLI_MODE_GENERIC,
            input.size(),
            input.data(),
            &encoded,
            output.data());
        if (!ok) return false;
        output.resize(encoded);
        return true;
    }

    bool Decompress(const std::vector<uint8_t>& input, std::vector<uint8_t>& output) override {
        BrotliDecoderState* state = BrotliDecoderCreateInstance(nullptr, nullptr, nullptr);
        if (!state) return false;

        size_t available_in = input.size();
        const uint8_t* next_in = input.data();
        std::vector<uint8_t> buffer(1024);
        BrotliDecoderResult result = BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT;
        output.clear();

        while (result == BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT || result == BROTLI_DECODER_RESULT_NEEDS_MORE_INPUT) {
            size_t available_out = buffer.size();
            uint8_t* next_out = buffer.data();
            result = BrotliDecoderDecompressStream(state, &available_in, &next_in, &available_out, &next_out, nullptr);
            const size_t produced = buffer.size() - available_out;
            output.insert(output.end(), buffer.data(), buffer.data() + produced);
            if (result == BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT) {
                buffer.resize(buffer.size() * 2);
            }
        }

        BrotliDecoderDestroyInstance(state);
        return result == BROTLI_DECODER_RESULT_SUCCESS;
    }
};

}  // namespace

void InitializeCompressionBackend(compression_backend_t backend) {
    if (backend == compression_backend_t::bsc) {
        std::call_once(g_bsc_init_flag, []() { CBSCWrapper::InitLibrary(p_bsc_features); });
    }
}

std::unique_ptr<CompressionStrategy> MakeCompressionStrategy(compression_backend_t backend, const bsc_params_t& bsc_params) {
    InitializeCompressionBackend(backend);
    switch (backend) {
    case compression_backend_t::bsc:
        return std::make_unique<BscCompressionStrategy>(bsc_params);
    case compression_backend_t::zstd:
        return std::make_unique<ZstdCompressionStrategy>();
    case compression_backend_t::brotli:
        return std::make_unique<BrotliCompressionStrategy>();
    default:
        return nullptr;
    }
}
