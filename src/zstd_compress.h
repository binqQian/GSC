#pragma once

#include <iostream>
#include <vector>
// #include <zstd.h>

// #include <thread>
#include "../include/zstd-1.5.2/lib/zstd.h"
using namespace std;
namespace zstd {
    // zstd 薄封装：统一项目里 vector / 指针两种调用方式。
    bool zstd_compress(const std::vector<uint8_t>& srcContent, std::vector<uint8_t>& cBuff);
    bool zstd_compress(vector<uint32_t>& v_text, std::vector<uint8_t>& cBuff, size_t& cSizeActual);
    bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint8_t>& dBuff);
    bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint32_t>& v_text, size_t& dSizeActual);
    // 直接从指针解压，主要给 mmap / on-disk payload 的高频路径使用。
    bool zstd_decompress_ptr(const uint8_t* src, size_t src_size, std::vector<uint8_t>& dBuff);
}
