#pragma once

#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace fmt_compress {

// PL count lookup table: getplcnt[allele_count] = PL element count
// PL elements = (n+1)*n/2 where n = allele_count
constexpr uint16_t kPlCountTable[12] = {0, 0, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66};

// Missing value marker
constexpr uint32_t kMissingValue = 0xFFFFFFFF;
constexpr uint16_t kMissingValue16 = 0xFFFE;

// AD dictionary item: max 48 bytes
// data_[0]: data length, data_[1]: type, data_[2~47]: payload
struct ADItem {
    uint8_t data_[48];

    bool operator==(const ADItem& other) const {
        return memcmp(data_, other.data_, data_[0] + 1) == 0;
    }
};

// PL dictionary item: max 128 bytes
// data_[0]: data length, data_[1]: type, data_[2~127]: payload
struct PLItem {
    uint8_t data_[128];

    bool operator==(const PLItem& other) const {
        return memcmp(data_, other.data_, data_[0] + 1) == 0;
    }
};

// PID dictionary item: max 512 bytes
// data_[0~1]: length (uint16), data_[2~511]: payload
struct PIDItem {
    uint8_t data_[512];

    bool operator==(const PIDItem& other) const {
        uint16_t len = 0;
        memcpy(&len, data_, 2);
        return memcmp(data_, other.data_, len + 2) == 0;
    }
};

// Abnormal item for storing exceptions (pos, value)
struct AbnormalItem {
    uint32_t val = 0;
    uint32_t pos = 0;

    AbnormalItem() = default;
    AbnormalItem(uint32_t v, uint32_t p) : val(v), pos(p) {}
};

// Parse comma-separated integers, "." -> kMissingValue
inline void splitString2uint(const char* str, int slen, char delimiter,
                             std::vector<uint32_t>& vec) {
    vec.clear();
    if (!str || slen <= 0) return;

    const char* p = str;
    const char* q = str;
    const char* end = str + slen;

    while (p < end) {
        if (*p == delimiter) {
            if (*q == '.') {
                vec.push_back(kMissingValue);
            } else {
                vec.push_back(static_cast<uint32_t>(atoi(q)));
            }
            q = p + 1;
        }
        p++;
    }

    if (q < end) {
        if (*q == '.') {
            vec.push_back(kMissingValue);
        } else {
            vec.push_back(static_cast<uint32_t>(atoi(q)));
        }
    }
}

// Parse and sum comma-separated integers
inline uint32_t splitString2uintSum(const char* str, int slen, char delimiter) {
    if (!str || slen <= 0) return kMissingValue;

    uint32_t sum = 0;
    const char* p = str;
    const char* q = str;
    const char* end = str + slen;

    while (p < end) {
        if (*p == delimiter) {
            if (*q == '.') {
                return kMissingValue;
            }
            sum += static_cast<uint32_t>(atoi(q));
            q = p + 1;
        }
        p++;
    }

    if (q < end) {
        if (*q == '.') {
            return kMissingValue;
        }
        sum += static_cast<uint32_t>(atoi(q));
    }

    return sum;
}

// Sum array elements (optimized for small arrays)
template <typename T>
inline uint32_t sumArray(uint32_t cnt, const T* arr) {
    uint32_t sum = 0;
    switch (cnt) {
        case 2: return arr[0] + arr[1];
        case 3: return arr[0] + arr[1] + arr[2];
        case 4: return arr[0] + arr[1] + arr[2] + arr[3];
        case 5: return arr[0] + arr[1] + arr[2] + arr[3] + arr[4];
        case 6: return arr[0] + arr[1] + arr[2] + arr[3] + arr[4] + arr[5];
        case 7: return arr[0] + arr[1] + arr[2] + arr[3] + arr[4] + arr[5] + arr[6];
        case 8: return arr[0] + arr[1] + arr[2] + arr[3] + arr[4] + arr[5] + arr[6] + arr[7];
        default:
            for (uint32_t i = 0; i < cnt; ++i) {
                sum += arr[i];
            }
            return sum;
    }
}

// Check PL pattern type:
// Type 0: all zeros
// Type 1: [0, a, b, a, b, b, a, b, b, b, ...] pattern
// Type 2: Type 1 with b == 15*a
// Type 3: no pattern (use dictionary)
template<typename T>
inline uint8_t checkPlPattern(const T* array, uint32_t size,
                              uint32_t& a_val, uint32_t& b_val) {
    // First element must be 0
    if (array[0] != 0) {
        return 3;
    }

    // Check if all zeros
    bool allZero = true;
    for (uint32_t i = 0; i < size; ++i) {
        if (array[i] != 0) {
            allZero = false;
            break;
        }
    }

    if (allZero) {
        a_val = 0;
        b_val = 0;
        return 0;
    }

    // Get a and b values from positions 1 and 2
    a_val = array[1];
    b_val = array[2];

    // Special cases for size 3 and 6
    if (size == 3) {
        return (a_val * 15 == b_val) ? 2 : 1;
    } else if (size == 6) {
        if (array[1] == array[3] && array[2] == array[4] && array[2] == array[5]) {
            return (a_val * 15 == b_val) ? 2 : 1;
        }
        return 3;
    }

    // Check pattern iteratively for larger arrays
    uint32_t pos = 3;
    uint32_t b_count = 2;

    while (pos < size) {
        // Should be a_val
        if (array[pos] != a_val) {
            return 3;
        }
        pos++;

        // Should be b_count consecutive b_val
        for (uint32_t i = 0; i < b_count && pos < size; ++i, ++pos) {
            if (array[pos] != b_val) {
                return 3;
            }
        }

        b_count++;
    }

    return (a_val * 15 == b_val) ? 2 : 1;
}

// Get max value and second minimum from array
template <typename T>
inline void getMaxAndSecondMin(const T* data, uint32_t cnt, T& max_val, T& second_min) {
    max_val = data[0];
    T smallest = data[0];
    second_min = std::numeric_limits<T>::max();

    for (uint32_t i = 1; i < cnt; ++i) {
        if (data[i] > max_val) {
            max_val = data[i];
        }

        if (data[i] < smallest) {
            second_min = smallest;
            smallest = data[i];
        } else if (data[i] < second_min && data[i] != smallest) {
            second_min = data[i];
        }
    }
}

// Get second minimum from array
template <typename T>
inline T getSecondMin(const T* data, uint32_t cnt) {
    T smallest = data[0];
    T second_min = std::numeric_limits<T>::max();

    for (uint32_t i = 1; i < cnt; ++i) {
        if (data[i] < smallest) {
            second_min = smallest;
            smallest = data[i];
        } else if (data[i] < second_min && data[i] != smallest) {
            second_min = data[i];
        }
    }

    return second_min;
}

// Pack 2-bit tips into byte array
inline uint32_t packTips(const std::vector<uint8_t>& tips, uint8_t* out) {
    uint32_t byte_pos = 0;
    memset(out, 0, (tips.size() + 7) / 8);

    for (size_t i = 0; i < tips.size(); i++) {
        byte_pos = i >> 3;
        uint8_t bit_pos = i & 7;
        if (tips[i]) {
            out[byte_pos] |= tips[i] << bit_pos;
        }
    }

    return byte_pos + 1;
}

// Unpack byte array to 2-bit tips
inline void unpackTips(const uint8_t* data, uint32_t cnt, std::vector<uint8_t>& tips) {
    tips.clear();
    tips.reserve(cnt);

    for (uint32_t i = 0; i < cnt; i++) {
        uint32_t byte_pos = i >> 3;
        uint8_t bit_pos = i & 7;
        tips.push_back((data[byte_pos] >> bit_pos) & 1);
    }
}

// Generate PL string from pattern (a, b)
inline void patternToPlString(uint32_t a, uint32_t b, uint8_t allele,
                              uint16_t pl_cnt, std::string& str) {
    char buf[512] = {0};

    switch (pl_cnt) {
        case 3:
            snprintf(buf, sizeof(buf), "0,%u,%u", a, b);
            str.append(buf);
            break;
        case 6:
            snprintf(buf, sizeof(buf), "0,%u,%u,%u,%u,%u", a, b, a, b, b);
            str.append(buf);
            break;
        case 10:
            snprintf(buf, sizeof(buf), "0,%u,%u,%u,%u,%u,%u,%u,%u,%u",
                     a, b, a, b, b, a, b, b, b);
            str.append(buf);
            break;
        case 15:
            snprintf(buf, sizeof(buf), "0,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u",
                     a, b, a, b, b, a, b, b, b, a, b, b, b, b);
            str.append(buf);
            break;
        default:
            // General case for larger allele counts
            str.append("0,");
            for (uint16_t i = 1; i < allele; i++) {
                str.append(std::to_string(a)).append(",");
                for (uint16_t j = 0; j < i; j++) {
                    str.append(std::to_string(b)).append(",");
                }
            }
            str.pop_back();  // Remove trailing comma
            break;
    }
}

}  // namespace fmt_compress

// Hash functions for dictionary items
namespace std {
    template<>
    struct hash<fmt_compress::ADItem> {
        size_t operator()(const fmt_compress::ADItem& item) const {
            // Simple FNV-1a hash
            size_t hash = 14695981039346656037ULL;
            for (uint8_t i = 0; i <= item.data_[0]; ++i) {
                hash ^= item.data_[i];
                hash *= 1099511628211ULL;
            }
            return hash;
        }
    };

    template<>
    struct hash<fmt_compress::PLItem> {
        size_t operator()(const fmt_compress::PLItem& item) const {
            size_t hash = 14695981039346656037ULL;
            for (uint8_t i = 0; i <= item.data_[0]; ++i) {
                hash ^= item.data_[i];
                hash *= 1099511628211ULL;
            }
            return hash;
        }
    };

    template<>
    struct hash<fmt_compress::PIDItem> {
        size_t operator()(const fmt_compress::PIDItem& item) const {
            uint16_t len = 0;
            memcpy(&len, item.data_, 2);
            size_t hash = 14695981039346656037ULL;
            for (uint16_t i = 0; i < len + 2; ++i) {
                hash ^= item.data_[i];
                hash *= 1099511628211ULL;
            }
            return hash;
        }
    };
}
