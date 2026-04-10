#include "vint_code.h"
#include <iostream>

namespace vint_code {


// 读取一个 vint；单个 '\0' 代表数值 0。
uint32_t ReadVint(std::vector<uint8_t>& buffer, size_t& pos)
{
    if(buffer[pos] == '\0')
    {
        pos++;
        return 0;
    }
    uint8_t firstByte = buffer[pos++];
    uint8_t mask = 0x80;
    uint32_t value = firstByte & 0x7f;
    int shift = 7;
    while (firstByte & mask)
    {
        firstByte = buffer[pos++];
        // value = (value << 7) | (firstByte & 0x7f);
        value = ((firstByte & 0x7f) << shift) | value ;
        shift += 7;
        // mask <<= 7;
    }

    return value;
}

uint32_t ReadVint(const uint8_t* buffer, size_t size, size_t& pos)
{
    // 指针版实现主要用于直接读取压缩 payload，避免中间复制。
    if (pos >= size)
        return 0;
    if (buffer[pos] == '\0')
    {
        pos++;
        return 0;
    }
    uint8_t firstByte = buffer[pos++];
    uint8_t mask = 0x80;
    uint32_t value = firstByte & 0x7f;
    int shift = 7;
    while (firstByte & mask)
    {
        if (pos >= size)
            break;
        firstByte = buffer[pos++];
        value = ((firstByte & 0x7f) << shift) | value;
        shift += 7;
    }
    return value;
}


size_t WriteVint(uint32_t value, std::vector<uint8_t>& buffer)
{
    // 这里采用 7 bit 一组的小端式拼接；
    // value==0 会被编码成 '\0'，这是调用方依赖的重要约定。
    size_t size = 0;

    while (value > 0x7f)
    {
        buffer.push_back((value & 0x7f) | 0x80);
        value >>= 7;
        size++;
    }
    if(value)
        buffer.push_back(value & 0x7f);
    else
        buffer.push_back('\0');
    size++;

    return size;
}
std::vector<uint8_t> EncodeArray(const std::vector<uint32_t>& arr)
{
    // 批量编码时不额外写长度，解码端直接顺序读到 buffer 结束。
    std::vector<uint8_t> buffer;

    buffer.reserve(arr.size() * 4);

    for (const auto& value : arr)
    {
        WriteVint(value, buffer);
    }

    return buffer;
}

void EncodeArray(const std::vector<uint32_t>& arr, std::vector<uint8_t>& out)
{
    out.clear();
    out.reserve(arr.size() * 4);
    for (const auto& value : arr)
        WriteVint(value, out);
}
std::vector<uint32_t> DecodeArray(std::vector<uint8_t> &buffer)
{
    size_t size = buffer.size();
    std::vector<uint32_t> arr;
    arr.reserve(size / 4);
    size_t pos = 0;

    while (pos < size)
    {
        arr.push_back(ReadVint(buffer, pos));
    }

    return arr;
}

void DecodeArray(const uint8_t* buffer, size_t size, std::vector<uint32_t>& out)
{
    out.clear();
    if (!buffer || size == 0)
        return;
    out.reserve(size / 4);
    size_t pos = 0;
    while (pos < size)
        out.push_back(ReadVint(buffer, size, pos));
}
};
