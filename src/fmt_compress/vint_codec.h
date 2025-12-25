#pragma once

#include <cstdint>
#include <cstring>

// EBML Variable Integer (Vint) Codec
// Big-endian variable-length integer encoding
// Format:
// 1xxx xxxx                                                    - 1 byte,  value 0 to 2^7-2
// 01xx xxxx xxxx xxxx                                          - 2 bytes, value 0 to 2^14-2
// 001x xxxx xxxx xxxx xxxx xxxx                                - 3 bytes, value 0 to 2^21-2
// 0001 xxxx xxxx xxxx xxxx xxxx xxxx xxxx                      - 4 bytes, value 0 to 2^28-2
// ...

namespace fmt_compress {

// Vint limits for each byte length
constexpr uint64_t kVintLimit[10] = {
    0, 0x7E, 0x3FFE, 0x1FFFFE, 0xFFFFFFE, 0x7FFFFFFFE,
    0x3FFFFFFFFFE, 0x1FFFFFFFFFFFE, 0xFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF
};

// Vint mask for each byte length
constexpr uint64_t kVintMask[10] = {
    0, 0x80, 0x4000, 0x200000, 0x10000000, 0x800000000,
    0x40000000000, 0x2000000000000, 0x100000000000000, 0xFFFFFFFFFFFFFFFF
};

class VintCodec {
public:
    // Get the number of bytes needed to encode a value
    static uint8_t getEncodedSize(uint64_t val) {
        uint8_t i = 1;
        while (val > kVintLimit[i]) {
            ++i;
        }
        return i;
    }

    // Encode a value to buffer, returns number of bytes written
    static uint8_t encode(uint64_t val, uint8_t* buf) {
        uint8_t len = getEncodedSize(val);
        uint64_t ret = kVintMask[len] | val;
        toBigEndian(len, ret, buf);
        return len;
    }

    // Encode with known length
    static void encode(uint64_t val, uint8_t len, uint8_t* buf) {
        uint64_t ret = kVintMask[len] | val;
        toBigEndian(len, ret, buf);
    }

    // Decode from buffer, returns number of bytes read
    static uint8_t decode(const uint8_t* buf, uint64_t& val) {
        uint8_t len = 0;
        val = 0;

        if (buf[0] & 0x80) {
            len = 1;
            val = buf[0] & 0x7f;
        } else if (buf[0] & 0x40) {
            len = 2;
            val = (static_cast<uint64_t>(buf[0] & 0x3f) << 8) | buf[1];
        } else if (buf[0] & 0x20) {
            len = 3;
            val = (static_cast<uint64_t>(buf[0] & 0x1f) << 16) |
                  (static_cast<uint64_t>(buf[1]) << 8) | buf[2];
        } else if (buf[0] & 0x10) {
            len = 4;
            val = (static_cast<uint64_t>(buf[0] & 0x0f) << 24) |
                  (static_cast<uint64_t>(buf[1]) << 16) |
                  (static_cast<uint64_t>(buf[2]) << 8) | buf[3];
        } else if (buf[0] & 0x08) {
            len = 5;
            val = (static_cast<uint64_t>(buf[0] & 0x07) << 32) |
                  (static_cast<uint64_t>(buf[1]) << 24) |
                  (static_cast<uint64_t>(buf[2]) << 16) |
                  (static_cast<uint64_t>(buf[3]) << 8) | buf[4];
        } else if (buf[0] & 0x04) {
            len = 6;
            val = (static_cast<uint64_t>(buf[0] & 0x03) << 40) |
                  (static_cast<uint64_t>(buf[1]) << 32) |
                  (static_cast<uint64_t>(buf[2]) << 24) |
                  (static_cast<uint64_t>(buf[3]) << 16) |
                  (static_cast<uint64_t>(buf[4]) << 8) | buf[5];
        } else if (buf[0] & 0x02) {
            len = 7;
            val = (static_cast<uint64_t>(buf[0] & 0x01) << 48) |
                  (static_cast<uint64_t>(buf[1]) << 40) |
                  (static_cast<uint64_t>(buf[2]) << 32) |
                  (static_cast<uint64_t>(buf[3]) << 24) |
                  (static_cast<uint64_t>(buf[4]) << 16) |
                  (static_cast<uint64_t>(buf[5]) << 8) | buf[6];
        } else if (buf[0] & 0x01) {
            len = 8;
            val = (static_cast<uint64_t>(buf[1]) << 48) |
                  (static_cast<uint64_t>(buf[2]) << 40) |
                  (static_cast<uint64_t>(buf[3]) << 32) |
                  (static_cast<uint64_t>(buf[4]) << 24) |
                  (static_cast<uint64_t>(buf[5]) << 16) |
                  (static_cast<uint64_t>(buf[6]) << 8) | buf[7];
        }
        return len;
    }

private:
    // Write value in big-endian format
    static void toBigEndian(uint8_t len, uint64_t val, uint8_t* buf) {
        switch (len) {
            case 1:
                buf[0] = val & 0xff;
                break;
            case 2:
                buf[0] = (val >> 8) & 0xff;
                buf[1] = val & 0xff;
                break;
            case 3:
                buf[0] = (val >> 16) & 0xff;
                buf[1] = (val >> 8) & 0xff;
                buf[2] = val & 0xff;
                break;
            case 4:
                buf[0] = (val >> 24) & 0xff;
                buf[1] = (val >> 16) & 0xff;
                buf[2] = (val >> 8) & 0xff;
                buf[3] = val & 0xff;
                break;
            case 5:
                buf[0] = (val >> 32) & 0xff;
                buf[1] = (val >> 24) & 0xff;
                buf[2] = (val >> 16) & 0xff;
                buf[3] = (val >> 8) & 0xff;
                buf[4] = val & 0xff;
                break;
            case 6:
                buf[0] = (val >> 40) & 0xff;
                buf[1] = (val >> 32) & 0xff;
                buf[2] = (val >> 24) & 0xff;
                buf[3] = (val >> 16) & 0xff;
                buf[4] = (val >> 8) & 0xff;
                buf[5] = val & 0xff;
                break;
            case 7:
                buf[0] = (val >> 48) & 0xff;
                buf[1] = (val >> 40) & 0xff;
                buf[2] = (val >> 32) & 0xff;
                buf[3] = (val >> 24) & 0xff;
                buf[4] = (val >> 16) & 0xff;
                buf[5] = (val >> 8) & 0xff;
                buf[6] = val & 0xff;
                break;
            case 8:
                buf[0] = (val >> 56) & 0xff;
                buf[1] = (val >> 48) & 0xff;
                buf[2] = (val >> 40) & 0xff;
                buf[3] = (val >> 32) & 0xff;
                buf[4] = (val >> 24) & 0xff;
                buf[5] = (val >> 16) & 0xff;
                buf[6] = (val >> 8) & 0xff;
                buf[7] = val & 0xff;
                break;
            default:
                break;
        }
    }
};

}  // namespace fmt_compress
