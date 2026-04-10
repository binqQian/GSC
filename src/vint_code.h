#pragma once
#include <vector>
#include <cstdint>
#include "defs.h"

// 这是主压缩链路使用的 vint 编码。
// 注意：它和 fmt_compress/vint_codec.h 不是同一种格式。
// 这里把数值 0 编码成单个 '\0'，在 GT 稀疏索引里常被当作行终止符。
namespace vint_code {
    // 从 vector 缓冲读取一个 vint；读到 '\0' 时返回 0。
    uint32_t ReadVint(std::vector<uint8_t>& buffer, size_t& pos);
    // 把一个整数写成 vint；value==0 时写入 '\0'。
    size_t WriteVint(uint32_t value, std::vector<uint8_t>& buffer);
    
    // 批量编码 / 解码接口。
    std::vector<uint8_t> EncodeArray(const std::vector<uint32_t>& arr);
    void EncodeArray(const std::vector<uint32_t>& arr, std::vector<uint8_t>& out);
    std::vector<uint32_t> DecodeArray(std::vector<uint8_t>& buffer);

    // 指针版读取接口，主要给解压端直接处理 mmap / payload 使用。
    uint32_t ReadVint(const uint8_t* buffer, size_t size, size_t& pos);
    void DecodeArray(const uint8_t* buffer, size_t size, std::vector<uint32_t>& out);
};
