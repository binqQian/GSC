/**
 * @file gvcf_field_decompress.cpp
 * @brief Implementation of gVCF field-level decompression
 */

#include "gvcf_field_decompress.h"
#include "../logger.h"
#include <cstring>

namespace gvcf {

// ============================================================================
// FieldDecompressor Base Implementation
// ============================================================================

bool FieldDecompressor::ApplyBackendDecompression(const CompressedField& input,
                                                  std::vector<uint8_t>& output) {
    if (input.data.empty()) {
        output.clear();
        return true;
    }

    if (!backend_) {
        output = input.data;
        return true;
    }

    // First byte is compression flag
    if (input.data.size() < 1) {
        output.clear();
        return true;
    }

    uint8_t is_compressed = input.data[0];
    if (is_compressed == 0) {
        // Data is NOT compressed, just copy (skip flag byte)
        output.assign(input.data.begin() + 1, input.data.end());
        return true;
    } else {
        // Data is compressed, decompress (skip flag byte)
        return backend_->Decompress(input.data.data() + 1, input.data.size() - 1, output);
    }
}

// ============================================================================
// ChromDecompressor Implementation
// ============================================================================

bool ChromDecompressor::Decompress(const CompressedField& input,
                                   std::vector<std::string>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    RLEResult rle_result;
    if (!RLEEncoder::Deserialize(decompressed, rle_result)) {
        return false;
    }

    return RLEEncoder::Decompress(rle_result, output);
}

// ============================================================================
// PosDecompressor Implementation
// ============================================================================

bool PosDecompressor::Decompress(const CompressedField& input,
                                 std::vector<uint64_t>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    DeltaResult delta_result;
    if (!DeltaEncoder::Deserialize(decompressed, delta_result)) {
        return false;
    }

    return DeltaEncoder::Decompress(delta_result, output);
}

// ============================================================================
// EndDecompressor Implementation
// ============================================================================

bool EndDecompressor::Decompress(const CompressedField& input,
                                 const std::vector<uint64_t>& pos_values,
                                 std::vector<int64_t>& output) {
    output.clear();
    if (input.data.empty() || pos_values.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    // New format: inference-based compression
    // [exception_count] [delta_indices...] [values...]
    size_t pos = 0;
    uint64_t exception_count = VarIntUtil::ReadVarUint(decompressed.data(),
                                                       decompressed.size(), pos);

    // Read exception indices (delta encoded)
    std::vector<uint32_t> exception_indices;
    exception_indices.reserve(exception_count);
    uint32_t prev_idx = 0;
    for (uint64_t i = 0; i < exception_count; ++i) {
        uint32_t delta = static_cast<uint32_t>(
            VarIntUtil::ReadVarUint(decompressed.data(), decompressed.size(), pos));
        prev_idx += delta;
        exception_indices.push_back(prev_idx);
    }

    // Read exception values (END - POS differences)
    std::vector<int64_t> exception_values;
    exception_values.reserve(exception_count);
    for (uint64_t i = 0; i < exception_count; ++i) {
        int64_t val = VarIntUtil::ReadVarInt(decompressed.data(),
                                              decompressed.size(), pos);
        exception_values.push_back(val);
    }

    // Reconstruct END values
    // Most END values are inferred: END = next_POS - 1
    // Exceptions are stored explicitly
    output.resize(pos_values.size());
    size_t exc_idx = 0;

    for (size_t i = 0; i < pos_values.size(); ++i) {
        bool is_exception = (exc_idx < exception_indices.size() &&
                            exception_indices[exc_idx] == i);

        if (is_exception) {
            // Use stored value: END = POS + diff
            output[i] = static_cast<int64_t>(pos_values[i]) + exception_values[exc_idx];
            exc_idx++;
        } else {
            // Infer from next POS: END = next_POS - 1
            if (i + 1 < pos_values.size()) {
                output[i] = static_cast<int64_t>(pos_values[i + 1]) - 1;
            } else {
                // Last record without exception - shouldn't happen
                // Fall back to END = POS
                output[i] = static_cast<int64_t>(pos_values[i]);
            }
        }
    }

    return true;
}

bool EndDecompressor::DecompressStandalone(const CompressedField& input,
                                           std::vector<int64_t>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    DeltaResult delta_result;
    if (!DeltaEncoder::Deserialize(decompressed, delta_result)) {
        return false;
    }

    return DeltaEncoder::DecompressSigned(delta_result, output);
}

// ============================================================================
// IdDecompressor Implementation
// ============================================================================

bool IdDecompressor::Decompress(const CompressedField& input,
                                std::vector<std::string>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    if (input.method == FieldCompressionMethod::MASK) {
        MaskResult mask_result;
        if (!MaskEncoder::Deserialize(decompressed, mask_result)) {
            return false;
        }
        return MaskEncoder::Decompress(mask_result, output);
    } else {
        DictResult dict_result;
        if (!DictEncoder::Deserialize(decompressed, dict_result)) {
            return false;
        }
        return DictEncoder::Decompress(dict_result, output);
    }
}

// ============================================================================
// RefDecompressor Implementation
// ============================================================================

bool RefDecompressor::Decompress(const CompressedField& input,
                                 std::vector<std::string>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    // New 2-bit encoding format
    if (input.method == FieldCompressionMethod::MASK) {
        // 2-bit encoding: A=00, C=01, G=10, T=11
        static const char bits_to_base[] = {'A', 'C', 'G', 'T'};

        size_t pos = 0;
        uint64_t exception_count = VarIntUtil::ReadVarUint(decompressed.data(),
                                                           decompressed.size(), pos);

        // Read exception indices (delta encoded)
        std::vector<uint32_t> exception_indices;
        exception_indices.reserve(exception_count);
        uint32_t prev_idx = 0;
        for (uint64_t i = 0; i < exception_count; ++i) {
            uint32_t delta = static_cast<uint32_t>(
                VarIntUtil::ReadVarUint(decompressed.data(), decompressed.size(), pos));
            prev_idx += delta;
            exception_indices.push_back(prev_idx);
        }

        // Read exception values
        std::vector<std::string> exception_values;
        exception_values.reserve(exception_count);
        for (uint64_t i = 0; i < exception_count; ++i) {
            uint64_t len = VarIntUtil::ReadVarUint(decompressed.data(),
                                                    decompressed.size(), pos);
            std::string val(reinterpret_cast<const char*>(&decompressed[pos]), len);
            pos += len;
            exception_values.push_back(val);
        }

        // Read packed bases size and data
        uint64_t packed_size = VarIntUtil::ReadVarUint(decompressed.data(),
                                                        decompressed.size(), pos);

        // Reconstruct output
        output.resize(input.original_count);
        size_t exc_idx = 0;
        size_t base_idx = 0;  // Index into packed bases

        for (uint32_t i = 0; i < input.original_count; ++i) {
            if (exc_idx < exception_indices.size() && exception_indices[exc_idx] == i) {
                output[i] = exception_values[exc_idx];
                exc_idx++;
            } else {
                // Unpack 2-bit base
                size_t byte_idx = base_idx / 4;
                int bit_offset = (base_idx % 4) * 2;
                if (pos + byte_idx < decompressed.size()) {
                    int bits = (decompressed[pos + byte_idx] >> bit_offset) & 0x3;
                    output[i] = std::string(1, bits_to_base[bits]);
                }
                base_idx++;
            }
        }

        return true;
    }

    // Legacy formats
    if (input.method == FieldCompressionMethod::RLE) {
        RLEResult rle_result;
        if (!RLEEncoder::Deserialize(decompressed, rle_result)) {
            return false;
        }
        return RLEEncoder::Decompress(rle_result, output);
    } else {
        DictResult dict_result;
        if (!DictEncoder::Deserialize(decompressed, dict_result)) {
            return false;
        }
        return DictEncoder::Decompress(dict_result, output);
    }
}

// ============================================================================
// AltDecompressor Implementation
// ============================================================================

bool AltDecompressor::Decompress(const CompressedField& input,
                                 std::vector<std::vector<std::string>>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    // Decompress first allele using mask
    size_t pos = 0;
    MaskResult mask_result;
    if (!MaskEncoder::DeserializeWithPos(decompressed, pos, mask_result)) {
        return false;
    }

    std::vector<std::string> first_alts;
    if (!MaskEncoder::Decompress(mask_result, first_alts)) {
        return false;
    }

    // Convert to vector of vectors (single allele per position)
    output.reserve(first_alts.size());
    for (const auto& alt : first_alts) {
        if (alt.empty()) {
            output.push_back({});
        } else {
            output.push_back({alt});
        }
    }

    // Read multi-allelic data if present
    if (pos < decompressed.size()) {
        uint8_t has_extra = decompressed[pos++];
        if (has_extra && pos < decompressed.size()) {
            uint64_t extra_count = VarIntUtil::ReadVarUint(decompressed, pos);
            for (uint64_t i = 0; i < extra_count && pos < decompressed.size(); ++i) {
                uint64_t record_idx = VarIntUtil::ReadVarUint(decompressed, pos);
                uint64_t alt_count = VarIntUtil::ReadVarUint(decompressed, pos);
                if (record_idx < output.size()) {
                    for (uint64_t j = 0; j < alt_count && pos < decompressed.size(); ++j) {
                        uint64_t alt_len = VarIntUtil::ReadVarUint(decompressed, pos);
                        if (pos + alt_len <= decompressed.size()) {
                            std::string alt_str(decompressed.begin() + pos,
                                                decompressed.begin() + pos + alt_len);
                            output[record_idx].push_back(alt_str);
                            pos += alt_len;
                        }
                    }
                }
            }
        }
    }

    return true;
}

// ============================================================================
// QualDecompressor Implementation
// ============================================================================

bool QualDecompressor::Decompress(const CompressedField& input,
                                  std::vector<float>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    size_t pos = 0;
    uint64_t count = VarIntUtil::ReadVarUint(decompressed, pos);

    output.reserve(count);
    for (uint64_t i = 0; i < count && pos + 4 <= decompressed.size(); ++i) {
        uint32_t bits = decompressed[pos] |
                       (static_cast<uint32_t>(decompressed[pos + 1]) << 8) |
                       (static_cast<uint32_t>(decompressed[pos + 2]) << 16) |
                       (static_cast<uint32_t>(decompressed[pos + 3]) << 24);
        pos += 4;

        float val;
        std::memcpy(&val, &bits, sizeof(float));
        output.push_back(val);
    }

    return output.size() == count;
}

// ============================================================================
// FilterDecompressor Implementation
// ============================================================================

bool FilterDecompressor::Decompress(const CompressedField& input,
                                    std::vector<std::string>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    switch (input.method) {
        case FieldCompressionMethod::RLE: {
            RLEResult rle_result;
            if (!RLEEncoder::Deserialize(decompressed, rle_result)) {
                return false;
            }
            return RLEEncoder::Decompress(rle_result, output);
        }
        case FieldCompressionMethod::MASK: {
            MaskResult mask_result;
            if (!MaskEncoder::Deserialize(decompressed, mask_result)) {
                return false;
            }
            return MaskEncoder::Decompress(mask_result, output);
        }
        default: {
            DictResult dict_result;
            if (!DictEncoder::Deserialize(decompressed, dict_result)) {
                return false;
            }
            return DictEncoder::Decompress(dict_result, output);
        }
    }
}

// ============================================================================
// GTDecompressor Implementation
// ============================================================================

bool GTDecompressor::Decompress(const CompressedField& mask_input,
                                const CompressedField& patches_input,
                                const CompressedField& phase_input,
                                std::vector<GenotypeData>& output) {
    // First decompress GT strings
    std::vector<std::string> gt_strings;
    if (!DecompressStrings(mask_input, patches_input, phase_input, gt_strings)) {
        return false;
    }

    // Convert to GenotypeData
    output.clear();
    output.reserve(gt_strings.size());
    for (const auto& str : gt_strings) {
        output.push_back(GenotypeData::FromString(str));
    }

    return true;
}

bool GTDecompressor::DecompressSingle(const CompressedField& input,
                                      std::vector<GenotypeData>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    MaskResult mask_result;
    if (!MaskEncoder::Deserialize(decompressed, mask_result)) {
        return false;
    }

    std::vector<std::string> gt_strings;
    if (!MaskEncoder::Decompress(mask_result, gt_strings)) {
        return false;
    }

    output.reserve(gt_strings.size());
    for (const auto& str : gt_strings) {
        output.push_back(GenotypeData::FromString(str));
    }

    return true;
}

bool GTDecompressor::DecompressStrings(const CompressedField& mask_input,
                                       const CompressedField& patches_input,
                                       const CompressedField& phase_input,
                                       std::vector<std::string>& output) {
    output.clear();

    // Decompress mask (GT values)
    std::vector<uint8_t> mask_decompressed;
    if (!ApplyBackendDecompression(mask_input, mask_decompressed)) {
        return false;
    }

    MaskResult mask_result;
    if (!MaskEncoder::Deserialize(mask_decompressed, mask_result)) {
        return false;
    }

    if (!MaskEncoder::Decompress(mask_result, output)) {
        return false;
    }

    // Decompress phase and apply
    if (!phase_input.data.empty()) {
        std::vector<uint8_t> phase_decompressed;
        if (!ApplyBackendDecompression(phase_input, phase_decompressed)) {
            return false;
        }

        size_t pos = 0;
        uint64_t count = VarIntUtil::ReadVarUint(phase_decompressed, pos);

        RLEByteResult phase_rle;
        std::vector<uint8_t> rle_data(phase_decompressed.begin() + pos,
                                      phase_decompressed.end());
        if (!RLEByteEncoder::Deserialize(rle_data, phase_rle)) {
            return false;
        }

        std::vector<uint8_t> phase_bitmap;
        if (!RLEByteEncoder::Decompress(phase_rle, phase_bitmap)) {
            return false;
        }

        // Apply phase information
        for (size_t i = 0; i < output.size() && i < count; ++i) {
            bool is_phased = BitmapUtil::GetBit(phase_bitmap, static_cast<uint32_t>(i));
            if (is_phased) {
                // Replace '/' with '|'
                size_t sep_pos = output[i].find('/');
                if (sep_pos != std::string::npos) {
                    output[i][sep_pos] = '|';
                }
            }
        }
    }

    return true;
}

// ============================================================================
// DPDecompressor Implementation
// ============================================================================

bool DPDecompressor::Decompress(const CompressedField& input,
                                std::vector<int32_t>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    if (input.method == FieldCompressionMethod::RLE) {
        RLEIntResult rle_result;
        if (!RLEIntEncoder::Deserialize(decompressed, rle_result)) {
            return false;
        }
        return RLEIntEncoder::Decompress(rle_result, output);
    } else {
        DictIntResult dict_result;
        if (!DictIntEncoder::Deserialize(decompressed, dict_result)) {
            return false;
        }
        return DictIntEncoder::Decompress(dict_result, output);
    }
}

// ============================================================================
// MinDPDecompressor Implementation
// ============================================================================

bool MinDPDecompressor::Decompress(const CompressedField& diff_input,
                                   const std::vector<int32_t>& dp,
                                   std::vector<int32_t>& output) {
    output.clear();
    if (diff_input.data.empty() || dp.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(diff_input, decompressed)) {
        return false;
    }

    std::vector<int32_t> diffs;

    if (diff_input.method == FieldCompressionMethod::MASK) {
        MaskIntResult mask_result;
        if (!MaskIntEncoder::Deserialize(decompressed, mask_result)) {
            return false;
        }
        if (!MaskIntEncoder::Decompress(mask_result, diffs)) {
            return false;
        }
    } else {
        RLEIntResult rle_result;
        if (!RLEIntEncoder::Deserialize(decompressed, rle_result)) {
            return false;
        }
        if (!RLEIntEncoder::Decompress(rle_result, diffs)) {
            return false;
        }
    }

    // Calculate MIN_DP = DP - diff
    if (diffs.size() != dp.size()) {
        return false;
    }

    output.reserve(diffs.size());
    for (size_t i = 0; i < diffs.size(); ++i) {
        output.push_back(dp[i] - diffs[i]);
    }

    return true;
}

bool MinDPDecompressor::DecompressStandalone(const CompressedField& input,
                                             std::vector<int32_t>& output) {
    DPDecompressor dp_decomp(backend_);
    return dp_decomp.Decompress(input, output);
}

// ============================================================================
// GQDecompressor Implementation
// ============================================================================

bool GQDecompressor::Decompress(const CompressedField& input,
                                const std::vector<std::vector<int32_t>>& pl,
                                std::vector<int32_t>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    size_t pos = 0;

    // Read original count
    uint64_t count = VarIntUtil::ReadVarUint(decompressed, pos);

    // Read bitmask
    uint64_t bitmask_size = VarIntUtil::ReadVarUint(decompressed, pos);
    RLEByteResult bitmask_rle;
    RLEByteEncoder::Deserialize(decompressed.data() + pos, bitmask_size, bitmask_rle);
    std::vector<uint8_t> bitmask;
    RLEByteEncoder::Decompress(bitmask_rle, bitmask);
    pos += bitmask_size;

    // Read exception count
    uint64_t exception_count = VarIntUtil::ReadVarUint(decompressed, pos);

    // Read exception indices
    std::vector<uint32_t> exception_indices;
    if (exception_count > 0) {
        uint64_t idx_size = VarIntUtil::ReadVarUint(decompressed, pos);
        DeltaResult delta_indices;
        DeltaEncoder::Deserialize(decompressed.data() + pos, idx_size, delta_indices);
        std::vector<uint64_t> indices_u64;
        DeltaEncoder::Decompress(delta_indices, indices_u64);
        exception_indices.assign(indices_u64.begin(), indices_u64.end());
        pos += idx_size;
    }

    // Read exception values
    std::vector<int32_t> exception_values;
    exception_values.reserve(exception_count);
    for (uint64_t i = 0; i < exception_count; ++i) {
        exception_values.push_back(static_cast<int32_t>(
            VarIntUtil::ReadVarInt(decompressed, pos)));
    }

    // Reconstruct GQ values
    output.resize(count);
    size_t exception_idx = 0;

    for (uint32_t i = 0; i < count; ++i) {
        if (BitmapUtil::GetBit(bitmask, i)) {
            // Predicted from PL
            output[i] = (i < pl.size()) ? GQCompressor::PredictGQFromPL(pl[i]) : 0;
        } else {
            // Exception value
            if (exception_idx < exception_values.size()) {
                output[i] = exception_values[exception_idx++];
            } else {
                output[i] = 0;
            }
        }
    }

    return true;
}

bool GQDecompressor::DecompressStandalone(const CompressedField& input,
                                          std::vector<int32_t>& output) {
    DPDecompressor dp_decomp(backend_);
    return dp_decomp.Decompress(input, output);
}

// ============================================================================
// PLDecompressor Implementation
// ============================================================================

bool PLDecompressor::Decompress(const CompressedField& input,
                                std::vector<std::vector<int32_t>>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    size_t pos = 0;

    // Check version
    uint8_t version = decompressed[pos++];

    if (version == 2) {
        // Version 2: Optimized PL encoding
        uint64_t total_count = VarIntUtil::ReadVarUint(decompressed, pos);

        // Read standard pattern bitmask
        uint64_t bitmask_size = VarIntUtil::ReadVarUint(decompressed, pos);
        RLEByteResult bitmask_rle;
        RLEByteEncoder::Deserialize(decompressed.data() + pos, bitmask_size, bitmask_rle);
        std::vector<uint8_t> is_standard;
        RLEByteEncoder::Decompress(bitmask_rle, is_standard);
        pos += bitmask_size;

        // Read standard PL values
        uint64_t standard_count = VarIntUtil::ReadVarUint(decompressed, pos);
        std::vector<int32_t> standard_pl1, standard_pl2;
        standard_pl1.reserve(standard_count);
        standard_pl2.reserve(standard_count);
        for (uint64_t i = 0; i < standard_count; ++i) {
            standard_pl1.push_back(static_cast<int32_t>(VarIntUtil::ReadVarInt(decompressed, pos)));
        }
        for (uint64_t i = 0; i < standard_count; ++i) {
            standard_pl2.push_back(static_cast<int32_t>(VarIntUtil::ReadVarInt(decompressed, pos)));
        }

        // Read exception dictionary
        uint64_t exception_dict_size = VarIntUtil::ReadVarUint(decompressed, pos);
        std::vector<std::vector<int32_t>> exception_dict;
        exception_dict.reserve(exception_dict_size);
        for (uint64_t i = 0; i < exception_dict_size; ++i) {
            uint64_t pl_len = VarIntUtil::ReadVarUint(decompressed, pos);
            std::vector<int32_t> pl;
            pl.reserve(pl_len);
            for (uint64_t j = 0; j < pl_len; ++j) {
                pl.push_back(static_cast<int32_t>(VarIntUtil::ReadVarInt(decompressed, pos)));
            }
            exception_dict.push_back(std::move(pl));
        }

        // Read exception indices
        uint64_t exception_count = VarIntUtil::ReadVarUint(decompressed, pos);
        std::vector<int32_t> exception_indices;
        if (exception_count > 0) {
            uint64_t idx_size = VarIntUtil::ReadVarUint(decompressed, pos);
            RLEIntResult rle_indices;
            RLEIntEncoder::Deserialize(decompressed.data() + pos, idx_size, rle_indices);
            RLEIntEncoder::Decompress(rle_indices, exception_indices);
        }

        // Reconstruct PL values
        output.reserve(total_count);
        size_t std_idx = 0, exc_idx = 0;
        for (uint64_t i = 0; i < total_count; ++i) {
            bool is_std = (i / 8 < is_standard.size()) &&
                          ((is_standard[i / 8] >> (i % 8)) & 1);
            if (is_std && std_idx < standard_pl1.size()) {
                output.push_back({0, standard_pl1[std_idx], standard_pl2[std_idx]});
                std_idx++;
            } else if (exc_idx < exception_indices.size()) {
                int32_t idx = exception_indices[exc_idx++];
                if (idx >= 0 && static_cast<size_t>(idx) < exception_dict.size()) {
                    output.push_back(exception_dict[idx]);
                } else {
                    output.push_back({});
                }
            } else {
                output.push_back({});
            }
        }

        return output.size() == total_count;
    }

    // Version 1: Original dictionary encoding (fallback)
    // First byte was version, re-read as total_count for old format
    pos = 0;
    uint64_t total_count = VarIntUtil::ReadVarUint(decompressed, pos);
    uint64_t dict_size = VarIntUtil::ReadVarUint(decompressed, pos);

    // Read dictionary
    std::vector<std::vector<int32_t>> dictionary;
    dictionary.reserve(dict_size);
    for (uint64_t i = 0; i < dict_size; ++i) {
        uint64_t pl_len = VarIntUtil::ReadVarUint(decompressed, pos);
        std::vector<int32_t> pl;
        pl.reserve(pl_len);
        for (uint64_t j = 0; j < pl_len; ++j) {
            pl.push_back(static_cast<int32_t>(VarIntUtil::ReadVarInt(decompressed, pos)));
        }
        dictionary.push_back(std::move(pl));
    }

    // Read indices
    uint64_t idx_size = VarIntUtil::ReadVarUint(decompressed, pos);
    RLEIntResult rle_indices;
    RLEIntEncoder::Deserialize(decompressed.data() + pos, idx_size, rle_indices);
    std::vector<int32_t> indices;
    RLEIntEncoder::Decompress(rle_indices, indices);

    // Reconstruct PL values
    output.reserve(total_count);
    for (int32_t idx : indices) {
        if (idx >= 0 && static_cast<size_t>(idx) < dictionary.size()) {
            output.push_back(dictionary[idx]);
        } else {
            output.push_back({});
        }
    }

    return output.size() == total_count;
}

// ============================================================================
// ADDecompressor Implementation
// ============================================================================

bool ADDecompressor::Decompress(const CompressedField& input,
                                std::vector<std::vector<int32_t>>& output) {
    PLDecompressor pl_decomp(backend_);
    return pl_decomp.Decompress(input, output);
}

// ============================================================================
// GenericFieldDecompressor Implementation
// ============================================================================

bool GenericFieldDecompressor::Decompress(const CompressedField& input,
                                          std::vector<std::string>& output) {
    output.clear();
    if (input.data.empty()) {
        return true;
    }

    std::vector<uint8_t> decompressed;
    if (!ApplyBackendDecompression(input, decompressed)) {
        return false;
    }

    DictResult dict_result;
    if (!DictEncoder::Deserialize(decompressed, dict_result)) {
        return false;
    }

    return DictEncoder::Decompress(dict_result, output);
}

// ============================================================================
// GVCFBlockDecompressor Implementation
// ============================================================================

GVCFBlockDecompressor::GVCFBlockDecompressor(std::shared_ptr<CompressionBackend> backend)
    : backend_(backend)
    , chrom_decomp_(backend)
    , pos_decomp_(backend)
    , end_decomp_(backend)
    , id_decomp_(backend)
    , ref_decomp_(backend)
    , alt_decomp_(backend)
    , qual_decomp_(backend)
    , filter_decomp_(backend)
    , gt_decomp_(backend)
    , dp_decomp_(backend)
    , min_dp_decomp_(backend)
    , gq_decomp_(backend)
    , pl_decomp_(backend)
    , ad_decomp_(backend)
    , generic_decomp_(backend) {}

bool GVCFBlockDecompressor::Decompress(const CompressedGVCFBlock& input,
                                       GVCFBlock& output) {
    auto logger = LogManager::Instance().Logger();
    logger->debug("GVCFBlockDecompressor::Decompress - start, variants={}, samples={}",
                 input.variant_count, input.sample_count);

    output.Clear();
    output.variant_count = input.variant_count;
    output.sample_count = input.sample_count;

    // Decompress position fields
    logger->debug("Decompressing CHROM...");
    if (!chrom_decomp_.Decompress(input.chrom, output.position.chrom)) {
        logger->error("Failed to decompress CHROM");
        return false;
    }
    logger->debug("CHROM done, count={}", output.position.chrom.size());

    logger->debug("Decompressing POS...");
    if (!pos_decomp_.Decompress(input.pos, output.position.pos)) {
        logger->error("Failed to decompress POS");
        return false;
    }
    logger->debug("POS done, count={}", output.position.pos.size());

    logger->debug("Decompressing ID...");
    if (!id_decomp_.Decompress(input.id, output.position.id)) {
        logger->error("Failed to decompress ID");
        return false;
    }
    logger->debug("ID done, count={}", output.position.id.size());

    // Decompress sequence fields
    logger->debug("Decompressing REF...");
    if (!ref_decomp_.Decompress(input.ref, output.sequence.ref)) {
        logger->error("Failed to decompress REF");
        return false;
    }
    logger->debug("REF done, count={}", output.sequence.ref.size());

    logger->debug("Decompressing ALT...");
    if (!alt_decomp_.Decompress(input.alt, output.sequence.alt)) {
        logger->error("Failed to decompress ALT");
        return false;
    }
    logger->debug("ALT done, count={}", output.sequence.alt.size());

    // Decompress quality fields
    logger->debug("Decompressing QUAL...");
    if (!qual_decomp_.Decompress(input.qual, output.quality.qual)) {
        logger->error("Failed to decompress QUAL");
        return false;
    }
    logger->debug("QUAL done, count={}", output.quality.qual.size());

    logger->debug("Decompressing FILTER...");
    if (!filter_decomp_.Decompress(input.filter, output.quality.filter)) {
        logger->error("Failed to decompress FILTER");
        return false;
    }
    logger->debug("FILTER done, count={}", output.quality.filter.size());

    // Decompress INFO/END
    if (input.has_end_field) {
        logger->debug("Decompressing END...");
        if (!end_decomp_.Decompress(input.info_end, output.position.pos, output.info.end)) {
            logger->error("Failed to decompress END");
            return false;
        }
        logger->debug("END done, count={}", output.info.end.size());
    }

    // Decompress GT
    logger->debug("Decompressing GT...");
    if (!gt_decomp_.Decompress(input.gt_mask, input.gt_patches, input.gt_phase,
                               output.sample.gt)) {
        logger->error("Failed to decompress GT");
        return false;
    }
    logger->debug("GT done, count={}", output.sample.gt.size());

    // Decompress DP first
    logger->debug("Decompressing DP...");
    if (!dp_decomp_.Decompress(input.dp, output.sample.dp)) {
        logger->error("Failed to decompress DP");
        return false;
    }
    logger->debug("DP done, count={}", output.sample.dp.size());

    // Decompress MIN_DP
    if (input.has_min_dp) {
        logger->debug("Decompressing MIN_DP...");
        if (!min_dp_decomp_.Decompress(input.dp_min_dp_diff, output.sample.dp,
                                       output.sample.min_dp)) {
            logger->error("Failed to decompress MIN_DP");
            return false;
        }
        logger->debug("MIN_DP done, count={}", output.sample.min_dp.size());
    }

    // Decompress PL first (needed for GQ)
    logger->debug("Decompressing PL...");
    if (!pl_decomp_.Decompress(input.pl, output.sample.pl)) {
        logger->error("Failed to decompress PL");
        return false;
    }
    logger->debug("PL done, count={}", output.sample.pl.size());

    // Decompress GQ with PL prediction
    logger->debug("Decompressing GQ...");
    if (!gq_decomp_.Decompress(input.gq, output.sample.pl, output.sample.gq)) {
        logger->error("Failed to decompress GQ");
        return false;
    }
    logger->debug("GQ done, count={}", output.sample.gq.size());

    // Decompress AD
    logger->debug("Decompressing AD...");
    if (!ad_decomp_.Decompress(input.ad, output.sample.ad)) {
        logger->error("Failed to decompress AD");
        return false;
    }
    logger->debug("AD done, count={}", output.sample.ad.size());

    // Decompress unknown fields
    logger->debug("Decompressing unknown INFO fields, count={}...", input.unknown_info.size());
    for (const auto& kv : input.unknown_info) {
        const auto& name = kv.first;
        const auto& field = kv.second;
        logger->debug("  Unknown INFO: {}", name);
        std::vector<std::string> values;
        if (!generic_decomp_.Decompress(field, values)) {
            logger->error("Failed to decompress unknown INFO field: {}", name);
            return false;
        }
        output.unknown_info[name] = std::move(values);
    }
    logger->debug("Unknown INFO done");

    logger->debug("Decompressing unknown FORMAT fields, count={}...", input.unknown_format.size());
    for (const auto& kv : input.unknown_format) {
        const auto& name = kv.first;
        const auto& field = kv.second;
        logger->debug("  Unknown FORMAT: {}", name);
        std::vector<std::string> values;
        if (!generic_decomp_.Decompress(field, values)) {
            logger->error("Failed to decompress unknown FORMAT field: {}", name);
            return false;
        }
        output.unknown_format[name] = std::move(values);
    }
    logger->debug("Unknown FORMAT done");

    logger->debug("GVCFBlockDecompressor::Decompress - complete");
    return true;
}

bool GVCFBlockDecompressor::DecompressPositionFields(const CompressedGVCFBlock& input,
                                                     GVCFBlock& output) {
    output.variant_count = input.variant_count;

    if (!chrom_decomp_.Decompress(input.chrom, output.position.chrom)) {
        return false;
    }

    if (!pos_decomp_.Decompress(input.pos, output.position.pos)) {
        return false;
    }

    if (!id_decomp_.Decompress(input.id, output.position.id)) {
        return false;
    }

    return true;
}

bool GVCFBlockDecompressor::DecompressSampleFields(const CompressedGVCFBlock& input,
                                                   GVCFBlock& output) {
    output.sample_count = input.sample_count;

    if (!gt_decomp_.Decompress(input.gt_mask, input.gt_patches, input.gt_phase,
                               output.sample.gt)) {
        return false;
    }

    if (!dp_decomp_.Decompress(input.dp, output.sample.dp)) {
        return false;
    }

    if (input.has_min_dp) {
        if (!min_dp_decomp_.Decompress(input.dp_min_dp_diff, output.sample.dp,
                                       output.sample.min_dp)) {
            return false;
        }
    }

    if (!pl_decomp_.Decompress(input.pl, output.sample.pl)) {
        return false;
    }

    if (!gq_decomp_.Decompress(input.gq, output.sample.pl, output.sample.gq)) {
        return false;
    }

    if (!ad_decomp_.Decompress(input.ad, output.sample.ad)) {
        return false;
    }

    return true;
}

bool GVCFBlockDecompressor::Validate(const GVCFBlock& original,
                                     const GVCFBlock& decompressed) {
    if (original.variant_count != decompressed.variant_count) {
        return false;
    }

    // Validate position fields
    if (original.position.chrom != decompressed.position.chrom) {
        return false;
    }
    if (original.position.pos != decompressed.position.pos) {
        return false;
    }
    if (original.position.id != decompressed.position.id) {
        return false;
    }

    // Validate sequence fields
    if (original.sequence.ref != decompressed.sequence.ref) {
        return false;
    }

    // Validate GT
    if (original.sample.gt.size() != decompressed.sample.gt.size()) {
        return false;
    }
    for (size_t i = 0; i < original.sample.gt.size(); ++i) {
        if (!(original.sample.gt[i] == decompressed.sample.gt[i])) {
            return false;
        }
    }

    // Validate numeric fields
    if (original.sample.dp != decompressed.sample.dp) {
        return false;
    }
    if (original.sample.gq != decompressed.sample.gq) {
        return false;
    }
    if (original.sample.min_dp != decompressed.sample.min_dp) {
        return false;
    }

    return true;
}

} // namespace gvcf
