/**
 * @file gvcf_field_compress.cpp
 * @brief Implementation of gVCF field-level compression
 */

#include "gvcf_field_compress.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <sstream>

namespace gvcf {

// ============================================================================
// GenotypeData Implementation
// ============================================================================

std::string GenotypeData::ToString() const {
    std::string result;
    if (allele1 < 0) {
        result = ".";
    } else {
        result = std::to_string(allele1);
    }

    result += phased ? '|' : '/';

    if (allele2 < 0) {
        result += ".";
    } else {
        result += std::to_string(allele2);
    }

    return result;
}

GenotypeData GenotypeData::FromString(const std::string& str) {
    GenotypeData gt;

    if (str.empty() || str == "." || str == "./." || str == ".|.") {
        return gt; // Missing
    }

    size_t sep_pos = str.find('/');
    gt.phased = false;
    if (sep_pos == std::string::npos) {
        sep_pos = str.find('|');
        gt.phased = (sep_pos != std::string::npos);
    }

    if (sep_pos == std::string::npos) {
        // Haploid or malformed
        if (str != ".") {
            gt.allele1 = static_cast<int8_t>(std::stoi(str));
        }
        return gt;
    }

    std::string a1_str = str.substr(0, sep_pos);
    std::string a2_str = str.substr(sep_pos + 1);

    gt.allele1 = (a1_str == ".") ? -1 : static_cast<int8_t>(std::stoi(a1_str));
    gt.allele2 = (a2_str == ".") ? -1 : static_cast<int8_t>(std::stoi(a2_str));

    return gt;
}

// ============================================================================
// FieldCompressor Base Implementation
// ============================================================================

bool FieldCompressor::ApplyBackendCompression(const std::vector<uint8_t>& input,
                                              CompressedField& output) {
    if (!backend_) {
        output.data = input;
        return true;
    }

    std::vector<uint8_t> compressed;
    if (!backend_->Compress(input, compressed)) {
        return false;
    }

    // Use compressed only if smaller (accounting for 1-byte flag)
    if (compressed.size() + 1 < input.size()) {
        output.data.clear();
        output.data.reserve(compressed.size() + 1);
        output.data.push_back(1);  // Flag: data is compressed
        output.data.insert(output.data.end(), compressed.begin(), compressed.end());
    } else {
        output.data.clear();
        output.data.reserve(input.size() + 1);
        output.data.push_back(0);  // Flag: data is NOT compressed
        output.data.insert(output.data.end(), input.begin(), input.end());
    }

    return true;
}

// ============================================================================
// ChromCompressor Implementation
// ============================================================================

bool ChromCompressor::Compress(const std::vector<std::string>& data,
                               CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::RLE;

    // RLE encode
    RLEResult rle_result;
    if (!RLEEncoder::Compress(data, rle_result)) {
        return false;
    }

    // Serialize
    std::vector<uint8_t> serialized;
    if (!RLEEncoder::Serialize(rle_result, serialized)) {
        return false;
    }

    // Apply backend compression
    if (!ApplyBackendCompression(serialized, output)) {
        return false;
    }

    output.compression_ratio = static_cast<float>(output.data.size()) /
        (data.size() * 5); // Estimate original size

    return true;
}

GVCFFieldAnalysis ChromCompressor::Analyze(const std::vector<std::string>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "CHROM";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    uint32_t max_freq = 0;
    for (const auto& kv : freq) {
        if (kv.second > max_freq) max_freq = kv.second;
    }
    result.dominant_ratio = static_cast<float>(max_freq) / result.total_count;

    result.recommended_method = FieldCompressionMethod::RLE;

    return result;
}

// ============================================================================
// PosCompressor Implementation
// ============================================================================

bool PosCompressor::Compress(const std::vector<uint64_t>& data,
                             CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::DELTA;

    // Delta encode
    DeltaResult delta_result;
    if (!DeltaEncoder::Compress(data, delta_result)) {
        return false;
    }

    // Serialize
    std::vector<uint8_t> serialized;
    if (!DeltaEncoder::Serialize(delta_result, serialized)) {
        return false;
    }

    // Apply backend compression
    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis PosCompressor::Analyze(const std::vector<uint64_t>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "POS";
    result.total_count = static_cast<uint32_t>(data.size());

    if (!data.empty()) {
        result.min_value = static_cast<int64_t>(data.front());
        result.max_value = static_cast<int64_t>(data.back());

        // Calculate average delta
        double sum_delta = 0;
        for (size_t i = 1; i < data.size(); ++i) {
            sum_delta += static_cast<double>(data[i]) - data[i-1];
        }
        result.avg_value = sum_delta / (data.size() - 1);
    }

    result.recommended_method = FieldCompressionMethod::DELTA;

    return result;
}

// ============================================================================
// EndCompressor Implementation
// ============================================================================

bool EndCompressor::Compress(const std::vector<int64_t>& end_values,
                             const std::vector<uint64_t>& pos_values,
                             CompressedField& output) {
    output.Clear();
    if (end_values.empty()) {
        return true;
    }

    if (end_values.size() != pos_values.size()) {
        return false;
    }

    output.original_count = static_cast<uint32_t>(end_values.size());
    output.method = FieldCompressionMethod::MASK;  // Using inference-based compression

    // Find exceptions: records where END != next_POS - 1
    // For continuous gVCF, 99.97% can be inferred
    std::vector<uint32_t> exception_indices;
    std::vector<int64_t> exception_values;  // Store END-POS diff for exceptions

    for (size_t i = 0; i < end_values.size(); ++i) {
        int64_t end_val = end_values[i];

        // For records without END field (stored as -1), treat as END = POS
        if (end_val < 0) {
            end_val = static_cast<int64_t>(pos_values[i]);
        }

        // Check if END can be inferred from next record's POS
        bool is_exception = true;
        if (i + 1 < end_values.size()) {
            // Can infer if END + 1 == next_POS
            int64_t next_pos = static_cast<int64_t>(pos_values[i + 1]);
            if (end_val + 1 == next_pos) {
                is_exception = false;
            }
        }
        // Last record is always an exception (no next POS to infer from)

        if (is_exception) {
            exception_indices.push_back(static_cast<uint32_t>(i));
            // Store END - POS difference (usually small positive number)
            exception_values.push_back(end_val - static_cast<int64_t>(pos_values[i]));
        }
    }

    // Serialize: [exception_count] [indices...] [values...]
    std::vector<uint8_t> serialized;

    // Exception count
    VarIntUtil::WriteVarUint(exception_indices.size(), serialized);

    // Exception indices (delta encoded for efficiency)
    if (!exception_indices.empty()) {
        // Delta encode indices
        std::vector<uint32_t> delta_indices;
        delta_indices.reserve(exception_indices.size());
        uint32_t prev = 0;
        for (uint32_t idx : exception_indices) {
            delta_indices.push_back(idx - prev);
            prev = idx;
        }
        for (uint32_t delta : delta_indices) {
            VarIntUtil::WriteVarUint(delta, serialized);
        }

        // Exception values (END - POS differences, usually small)
        for (int64_t val : exception_values) {
            VarIntUtil::WriteVarInt(val, serialized);
        }
    }

    return ApplyBackendCompression(serialized, output);
}

bool EndCompressor::CompressStandalone(const std::vector<int64_t>& data,
                                       CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::DELTA;

    DeltaResult delta_result;
    if (!DeltaEncoder::CompressSigned(data, delta_result)) {
        return false;
    }

    std::vector<uint8_t> serialized;
    if (!DeltaEncoder::Serialize(delta_result, serialized)) {
        return false;
    }

    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis EndCompressor::Analyze(const std::vector<int64_t>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "END";
    result.total_count = static_cast<uint32_t>(data.size());

    if (!data.empty()) {
        auto minmax = std::minmax_element(data.begin(), data.end()); auto min_it = minmax.first; auto max_it = minmax.second;
        result.min_value = *min_it;
        result.max_value = *max_it;
        result.avg_value = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }

    result.recommended_method = FieldCompressionMethod::DELTA;

    return result;
}

// ============================================================================
// IdCompressor Implementation
// ============================================================================

bool IdCompressor::Compress(const std::vector<std::string>& data,
                            CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());

    // Check dominant ratio for "."
    uint32_t dot_count = std::count(data.begin(), data.end(), ".");
    float dot_ratio = static_cast<float>(dot_count) / data.size();

    if (dot_ratio >= config_.mask_threshold) {
        // Use mask encoding with "." as dominant
        output.method = FieldCompressionMethod::MASK;

        MaskResult mask_result;
        if (!MaskEncoder::CompressWithDominant(data, ".", mask_result)) {
            return false;
        }

        output.dominant_ratio = mask_result.dominant_ratio;

        std::vector<uint8_t> serialized;
        if (!MaskEncoder::Serialize(mask_result, serialized)) {
            return false;
        }

        return ApplyBackendCompression(serialized, output);
    } else {
        // Use dictionary encoding
        output.method = FieldCompressionMethod::DICTIONARY;

        DictResult dict_result;
        if (!DictEncoder::Compress(data, dict_result)) {
            return false;
        }

        std::vector<uint8_t> serialized;
        if (!DictEncoder::Serialize(dict_result, serialized)) {
            return false;
        }

        return ApplyBackendCompression(serialized, output);
    }
}

GVCFFieldAnalysis IdCompressor::Analyze(const std::vector<std::string>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "ID";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    uint32_t dot_count = 0;
    if (freq.count(".")) {
        dot_count = freq["."];
    }
    result.dominant_ratio = static_cast<float>(dot_count) / result.total_count;

    result.recommended_method = (result.dominant_ratio >= 0.7f) ?
        FieldCompressionMethod::MASK : FieldCompressionMethod::DICTIONARY;

    return result;
}

// ============================================================================
// RefCompressor Implementation
// ============================================================================

bool RefCompressor::Compress(const std::vector<std::string>& data,
                             CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::MASK;  // 2-bit encoding for single bases

    // Count single-base vs multi-base REFs
    // In gVCF, 99.9%+ are single bases (A, C, G, T)
    std::vector<uint32_t> exception_indices;
    std::vector<std::string> exception_values;

    // 2-bit encoding: A=00, C=01, G=10, T=11
    auto base_to_bits = [](char c) -> int {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    };

    // First pass: identify exceptions
    for (size_t i = 0; i < data.size(); ++i) {
        const std::string& ref = data[i];
        if (ref.size() != 1 || base_to_bits(ref[0]) < 0) {
            exception_indices.push_back(static_cast<uint32_t>(i));
            exception_values.push_back(ref);
        }
    }

    // Serialize
    std::vector<uint8_t> serialized;

    // Exception count
    VarIntUtil::WriteVarUint(exception_indices.size(), serialized);

    // Exception indices (delta encoded)
    if (!exception_indices.empty()) {
        uint32_t prev = 0;
        for (uint32_t idx : exception_indices) {
            VarIntUtil::WriteVarUint(idx - prev, serialized);
            prev = idx;
        }

        // Exception values (length-prefixed strings)
        for (const auto& val : exception_values) {
            VarIntUtil::WriteVarUint(val.size(), serialized);
            serialized.insert(serialized.end(), val.begin(), val.end());
        }
    }

    // Pack single bases into 2-bit encoding (4 bases per byte)
    // Only for non-exception positions
    std::vector<uint8_t> packed_bases;
    packed_bases.reserve((data.size() + 3) / 4);

    uint8_t current_byte = 0;
    int bit_pos = 0;

    for (size_t i = 0; i < data.size(); ++i) {
        // Skip exceptions
        bool is_exception = std::binary_search(exception_indices.begin(),
                                                exception_indices.end(),
                                                static_cast<uint32_t>(i));
        if (is_exception) continue;

        int bits = base_to_bits(data[i][0]);
        current_byte |= (bits << bit_pos);
        bit_pos += 2;

        if (bit_pos == 8) {
            packed_bases.push_back(current_byte);
            current_byte = 0;
            bit_pos = 0;
        }
    }

    // Flush remaining bits
    if (bit_pos > 0) {
        packed_bases.push_back(current_byte);
    }

    // Append packed bases
    VarIntUtil::WriteVarUint(packed_bases.size(), serialized);
    serialized.insert(serialized.end(), packed_bases.begin(), packed_bases.end());

    return ApplyBackendCompression(serialized, output);
}

FieldCompressionMethod RefCompressor::SelectMethod(const std::vector<std::string>& data) {
    if (data.empty()) {
        return FieldCompressionMethod::DICTIONARY;
    }

    // Count unique values
    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    float unique_ratio = static_cast<float>(freq.size()) / data.size();

    // Check for run patterns
    uint32_t run_count = 1;
    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] != data[i-1]) {
            run_count++;
        }
    }

    float run_ratio = static_cast<float>(run_count) / data.size();

    // If low unique ratio and good run patterns, use RLE
    if (unique_ratio < 0.3f && run_ratio < 0.3f) {
        return FieldCompressionMethod::RLE;
    }

    return FieldCompressionMethod::DICTIONARY;
}

GVCFFieldAnalysis RefCompressor::Analyze(const std::vector<std::string>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "REF";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    uint32_t max_freq = 0;
    for (const auto& kv : freq) { const auto& count = kv.second;
        if (count > max_freq) max_freq = count;
    }
    result.dominant_ratio = static_cast<float>(max_freq) / result.total_count;

    result.recommended_method = SelectMethod(data);

    return result;
}

// ============================================================================
// AltCompressor Implementation
// ============================================================================

bool AltCompressor::Compress(const std::vector<std::vector<std::string>>& data,
                             CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::MASK;

    // Flatten ALT to first allele only for primary compression
    std::vector<std::string> first_alts;
    first_alts.reserve(data.size());
    for (const auto& alts : data) {
        first_alts.push_back(alts.empty() ? "" : alts[0]);
    }

    // Collect multi-allelic records (records with more than 1 ALT)
    std::vector<std::pair<uint32_t, std::vector<std::string>>> extra_alts;
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i].size() > 1) {
            std::vector<std::string> extras(data[i].begin() + 1, data[i].end());
            extra_alts.emplace_back(static_cast<uint32_t>(i), std::move(extras));
        }
    }

    // Check for <NON_REF> dominance (typical gVCF)
    uint32_t non_ref_count = std::count(first_alts.begin(), first_alts.end(),
                                        config_.default_alt);
    float non_ref_ratio = static_cast<float>(non_ref_count) / first_alts.size();

    output.dominant_ratio = non_ref_ratio;

    MaskResult mask_result;
    if (non_ref_ratio >= config_.mask_threshold) {
        if (!MaskEncoder::CompressWithDominant(first_alts, config_.default_alt,
                                               mask_result)) {
            return false;
        }
    } else {
        if (!MaskEncoder::Compress(first_alts, mask_result)) {
            return false;
        }
    }

    std::vector<uint8_t> serialized;
    if (!MaskEncoder::Serialize(mask_result, serialized)) {
        return false;
    }

    // Append multi-allelic data
    // Format: [has_extra: 1 byte] [if has_extra: extra_count + records]
    if (extra_alts.empty()) {
        serialized.push_back(0);  // No extra ALTs
    } else {
        serialized.push_back(1);  // Has extra ALTs
        VarIntUtil::WriteVarUint(extra_alts.size(), serialized);
        for (const auto& record : extra_alts) {
            VarIntUtil::WriteVarUint(record.first, serialized);  // Record index
            VarIntUtil::WriteVarUint(record.second.size(), serialized);  // Extra ALT count
            for (const auto& alt : record.second) {
                VarIntUtil::WriteVarUint(alt.size(), serialized);  // ALT length
                serialized.insert(serialized.end(), alt.begin(), alt.end());  // ALT data
            }
        }
    }

    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis AltCompressor::Analyze(const std::vector<std::vector<std::string>>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "ALT";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& alts : data) {
        for (const auto& alt : alts) {
            freq[alt]++;
        }
    }

    result.unique_count = static_cast<uint32_t>(freq.size());

    // Check <NON_REF> frequency
    uint32_t non_ref_count = freq.count("<NON_REF>") ? freq["<NON_REF>"] : 0;
    result.dominant_ratio = static_cast<float>(non_ref_count) / result.total_count;

    result.recommended_method = (result.dominant_ratio >= 0.7f) ?
        FieldCompressionMethod::MASK : FieldCompressionMethod::DICTIONARY;

    return result;
}

// ============================================================================
// QualCompressor Implementation
// ============================================================================

bool QualCompressor::Compress(const std::vector<float>& data,
                              CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::NONE;

    // Store as raw float bytes with simple encoding
    std::vector<uint8_t> serialized;
    serialized.reserve(data.size() * sizeof(float) + 4);

    // Write count
    VarIntUtil::WriteVarUint(data.size(), serialized);

    // Write floats
    for (float val : data) {
        uint32_t bits;
        std::memcpy(&bits, &val, sizeof(float));
        serialized.push_back((bits >> 0) & 0xFF);
        serialized.push_back((bits >> 8) & 0xFF);
        serialized.push_back((bits >> 16) & 0xFF);
        serialized.push_back((bits >> 24) & 0xFF);
    }

    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis QualCompressor::Analyze(const std::vector<float>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "QUAL";
    result.total_count = static_cast<uint32_t>(data.size());

    if (!data.empty()) {
        auto minmax = std::minmax_element(data.begin(), data.end()); auto min_it = minmax.first; auto max_it = minmax.second;
        result.min_value = static_cast<int64_t>(*min_it);
        result.max_value = static_cast<int64_t>(*max_it);
        result.avg_value = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }

    result.recommended_method = FieldCompressionMethod::NONE;

    return result;
}

// ============================================================================
// FilterCompressor Implementation
// ============================================================================

bool FilterCompressor::Compress(const std::vector<std::string>& data,
                                CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());

    // Count unique values
    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    // Find dominant value
    std::string dominant;
    uint32_t max_count = 0;
    for (const auto& kv : freq) { const auto& val = kv.first; const auto& count = kv.second;
        if (count > max_count) {
            max_count = count;
            dominant = val;
        }
    }

    float dominant_ratio = static_cast<float>(max_count) / data.size();
    output.dominant_ratio = dominant_ratio;

    std::vector<uint8_t> serialized;

    if (freq.size() == 1) {
        // All same value - use RLE
        output.method = FieldCompressionMethod::RLE;
        RLEResult rle_result;
        RLEEncoder::Compress(data, rle_result);
        RLEEncoder::Serialize(rle_result, serialized);
    } else if (dominant_ratio >= config_.mask_threshold) {
        // Use mask encoding
        output.method = FieldCompressionMethod::MASK;
        MaskResult mask_result;
        MaskEncoder::CompressWithDominant(data, dominant, mask_result);
        MaskEncoder::Serialize(mask_result, serialized);
    } else {
        // Use dictionary
        output.method = FieldCompressionMethod::DICTIONARY;
        DictResult dict_result;
        DictEncoder::Compress(data, dict_result);
        DictEncoder::Serialize(dict_result, serialized);
    }

    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis FilterCompressor::Analyze(const std::vector<std::string>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "FILTER";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    uint32_t max_freq = 0;
    for (const auto& kv : freq) { const auto& count = kv.second;
        if (count > max_freq) max_freq = count;
    }
    result.dominant_ratio = static_cast<float>(max_freq) / result.total_count;

    if (result.unique_count == 1) {
        result.recommended_method = FieldCompressionMethod::RLE;
    } else if (result.dominant_ratio >= 0.7f) {
        result.recommended_method = FieldCompressionMethod::MASK;
    } else {
        result.recommended_method = FieldCompressionMethod::DICTIONARY;
    }

    return result;
}

// ============================================================================
// GTCompressor Implementation
// ============================================================================

bool GTCompressor::Compress(const std::vector<GenotypeData>& data,
                           CompressedField& mask_output,
                           CompressedField& patches_output,
                           CompressedField& phase_output) {
    mask_output.Clear();
    patches_output.Clear();
    phase_output.Clear();

    if (data.empty()) {
        return true;
    }

    uint32_t count = static_cast<uint32_t>(data.size());
    mask_output.original_count = count;
    patches_output.original_count = count;
    phase_output.original_count = count;

    // Convert to string for mask encoding
    std::vector<std::string> gt_strings;
    gt_strings.reserve(data.size());
    for (const auto& gt : data) {
        gt_strings.push_back(gt.ToString());
    }

    // Compress GT strings with mask (dominant = "0/0")
    mask_output.method = FieldCompressionMethod::MASK;

    MaskResult mask_result;
    if (!MaskEncoder::CompressWithDominant(gt_strings, config_.default_gt, mask_result)) {
        return false;
    }

    mask_output.dominant_ratio = mask_result.dominant_ratio;

    // Serialize mask
    std::vector<uint8_t> mask_serialized;
    if (!MaskEncoder::Serialize(mask_result, mask_serialized)) {
        return false;
    }

    if (!ApplyBackendCompression(mask_serialized, mask_output)) {
        return false;
    }

    // Extract and compress phase information
    std::vector<bool> phases;
    phases.reserve(data.size());
    for (const auto& gt : data) {
        phases.push_back(gt.phased);
    }

    phase_output.method = FieldCompressionMethod::RLE;

    std::vector<uint8_t> phase_bitmap = BitmapUtil::FromBools(phases);
    RLEByteResult phase_rle;
    RLEByteEncoder::Compress(phase_bitmap, phase_rle);

    std::vector<uint8_t> phase_serialized;
    VarIntUtil::WriteVarUint(count, phase_serialized); // Original count
    std::vector<uint8_t> rle_buf;
    RLEByteEncoder::Serialize(phase_rle, rle_buf);
    phase_serialized.insert(phase_serialized.end(), rle_buf.begin(), rle_buf.end());

    if (!ApplyBackendCompression(phase_serialized, phase_output)) {
        return false;
    }

    return true;
}

bool GTCompressor::CompressSingle(const std::vector<GenotypeData>& data,
                                  CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::MASK;

    // Convert to string
    std::vector<std::string> gt_strings;
    gt_strings.reserve(data.size());
    for (const auto& gt : data) {
        gt_strings.push_back(gt.ToString());
    }

    MaskResult mask_result;
    if (!MaskEncoder::CompressWithDominant(gt_strings, config_.default_gt, mask_result)) {
        return false;
    }

    output.dominant_ratio = mask_result.dominant_ratio;

    std::vector<uint8_t> serialized;
    if (!MaskEncoder::Serialize(mask_result, serialized)) {
        return false;
    }

    return ApplyBackendCompression(serialized, output);
}

bool GTCompressor::CompressStrings(const std::vector<std::string>& data,
                                   CompressedField& mask_output,
                                   CompressedField& patches_output,
                                   CompressedField& phase_output) {
    std::vector<GenotypeData> gt_data;
    gt_data.reserve(data.size());
    for (const auto& str : data) {
        gt_data.push_back(GenotypeData::FromString(str));
    }
    return Compress(gt_data, mask_output, patches_output, phase_output);
}

GVCFFieldAnalysis GTCompressor::Analyze(const std::vector<GenotypeData>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "GT";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& gt : data) {
        freq[gt.ToString()]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    // Check 0/0 frequency
    uint32_t hom_ref_count = freq.count("0/0") ? freq["0/0"] : 0;
    hom_ref_count += freq.count("0|0") ? freq["0|0"] : 0;
    result.dominant_ratio = static_cast<float>(hom_ref_count) / result.total_count;

    result.recommended_method = FieldCompressionMethod::MASK;

    return result;
}

GVCFFieldAnalysis GTCompressor::AnalyzeStrings(const std::vector<std::string>& data) {
    std::vector<GenotypeData> gt_data;
    gt_data.reserve(data.size());
    for (const auto& str : data) {
        gt_data.push_back(GenotypeData::FromString(str));
    }
    return Analyze(gt_data);
}

// ============================================================================
// DPCompressor Implementation
// ============================================================================

bool DPCompressor::Compress(const std::vector<int32_t>& data,
                            CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());

    // Try RLE first
    RLEIntResult rle_result;
    RLEIntEncoder::Compress(data, rle_result);
    std::vector<uint8_t> rle_buf;
    RLEIntEncoder::Serialize(rle_result, rle_buf);

    // Try dictionary
    DictIntResult dict_result;
    DictIntEncoder::Compress(data, dict_result);
    std::vector<uint8_t> dict_buf;
    DictIntEncoder::Serialize(dict_result, dict_buf);

    // Choose smaller
    if (rle_buf.size() <= dict_buf.size()) {
        output.method = FieldCompressionMethod::RLE;
        return ApplyBackendCompression(rle_buf, output);
    } else {
        output.method = FieldCompressionMethod::DICTIONARY;
        return ApplyBackendCompression(dict_buf, output);
    }
}

GVCFFieldAnalysis DPCompressor::Analyze(const std::vector<int32_t>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "DP";
    result.total_count = static_cast<uint32_t>(data.size());

    if (!data.empty()) {
        auto minmax = std::minmax_element(data.begin(), data.end()); auto min_it = minmax.first; auto max_it = minmax.second;
        result.min_value = *min_it;
        result.max_value = *max_it;
        result.avg_value = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }

    std::unordered_map<int32_t, uint32_t> freq;
    for (int32_t val : data) {
        freq[val]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    result.recommended_method = FieldCompressionMethod::ADAPTIVE;

    return result;
}

// ============================================================================
// MinDPCompressor Implementation
// ============================================================================

bool MinDPCompressor::Compress(const std::vector<int32_t>& min_dp,
                               const std::vector<int32_t>& dp,
                               CompressedField& diff_output,
                               CompressedField& dp_output) {
    diff_output.Clear();
    dp_output.Clear();

    if (min_dp.empty() || dp.empty()) {
        return true;
    }

    if (min_dp.size() != dp.size()) {
        return false;
    }

    uint32_t count = static_cast<uint32_t>(min_dp.size());
    diff_output.original_count = count;
    dp_output.original_count = count;

    // First compress DP
    DPCompressor dp_comp(backend_, config_);
    if (!dp_comp.Compress(dp, dp_output)) {
        return false;
    }

    // Calculate DP - MIN_DP differences
    std::vector<int32_t> diffs;
    diffs.reserve(min_dp.size());
    for (size_t i = 0; i < min_dp.size(); ++i) {
        diffs.push_back(dp[i] - min_dp[i]);
    }

    // Analyze difference distribution
    uint32_t zero_count = std::count(diffs.begin(), diffs.end(), 0);
    float zero_ratio = static_cast<float>(zero_count) / diffs.size();

    if (zero_ratio >= 0.6f) {
        // Many zeros - use mask encoding
        diff_output.method = FieldCompressionMethod::MASK;

        MaskIntResult mask_result;
        MaskIntEncoder::CompressWithDominant(diffs, 0, mask_result);
        diff_output.dominant_ratio = mask_result.dominant_ratio;

        std::vector<uint8_t> serialized;
        MaskIntEncoder::Serialize(mask_result, serialized);
        return ApplyBackendCompression(serialized, diff_output);
    } else {
        // Use RLE or dictionary
        diff_output.method = FieldCompressionMethod::RLE;

        RLEIntResult rle_result;
        RLEIntEncoder::Compress(diffs, rle_result);

        std::vector<uint8_t> serialized;
        RLEIntEncoder::Serialize(rle_result, serialized);
        return ApplyBackendCompression(serialized, diff_output);
    }
}

bool MinDPCompressor::CompressStandalone(const std::vector<int32_t>& data,
                                         CompressedField& output) {
    DPCompressor dp_comp(backend_, config_);
    return dp_comp.Compress(data, output);
}

GVCFFieldAnalysis MinDPCompressor::Analyze(const std::vector<int32_t>& min_dp,
                                           const std::vector<int32_t>& dp) {
    GVCFFieldAnalysis result;
    result.field_name = "MIN_DP";
    result.total_count = static_cast<uint32_t>(min_dp.size());

    // Analyze differences
    std::vector<int32_t> diffs;
    for (size_t i = 0; i < std::min(min_dp.size(), dp.size()); ++i) {
        diffs.push_back(dp[i] - min_dp[i]);
    }

    if (!diffs.empty()) {
        auto minmax = std::minmax_element(diffs.begin(), diffs.end()); auto min_it = minmax.first; auto max_it = minmax.second;
        result.min_value = *min_it;
        result.max_value = *max_it;
        result.avg_value = std::accumulate(diffs.begin(), diffs.end(), 0.0) / diffs.size();
    }

    uint32_t zero_count = std::count(diffs.begin(), diffs.end(), 0);
    result.dominant_ratio = static_cast<float>(zero_count) / result.total_count;

    result.recommended_method = (result.dominant_ratio >= 0.6f) ?
        FieldCompressionMethod::MASK : FieldCompressionMethod::RLE;

    return result;
}

// ============================================================================
// GQCompressor Implementation
// ============================================================================

int32_t GQCompressor::PredictGQFromPL(const std::vector<int32_t>& pl) {
    if (pl.size() < 2) {
        return 0;
    }

    // GQ is typically the second smallest PL value
    std::vector<int32_t> sorted_pl = pl;
    std::partial_sort(sorted_pl.begin(), sorted_pl.begin() + 2, sorted_pl.end());

    return sorted_pl[1] - sorted_pl[0];
}

bool GQCompressor::Compress(const std::vector<int32_t>& gq,
                            const std::vector<std::vector<int32_t>>& pl,
                            CompressedField& output) {
    output.Clear();
    if (gq.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(gq.size());
    output.method = FieldCompressionMethod::MASK;

    // Predict GQ from PL and store exceptions
    std::vector<int32_t> exceptions;
    std::vector<uint32_t> exception_indices;
    std::vector<uint8_t> bitmask((gq.size() + 7) / 8, 0);

    for (size_t i = 0; i < gq.size(); ++i) {
        int32_t predicted = (i < pl.size()) ? PredictGQFromPL(pl[i]) : 0;

        if (gq[i] == predicted) {
            BitmapUtil::SetBit(bitmask, static_cast<uint32_t>(i), true);
        } else {
            exceptions.push_back(gq[i]);
            exception_indices.push_back(static_cast<uint32_t>(i));
        }
    }

    output.dominant_ratio = 1.0f - static_cast<float>(exceptions.size()) / gq.size();

    // Serialize
    std::vector<uint8_t> serialized;

    // Write original count
    VarIntUtil::WriteVarUint(gq.size(), serialized);

    // Write bitmask (RLE compressed)
    RLEByteResult bitmask_rle;
    RLEByteEncoder::Compress(bitmask, bitmask_rle);
    std::vector<uint8_t> bitmask_buf;
    RLEByteEncoder::Serialize(bitmask_rle, bitmask_buf);

    VarIntUtil::WriteVarUint(bitmask_buf.size(), serialized);
    serialized.insert(serialized.end(), bitmask_buf.begin(), bitmask_buf.end());

    // Write exceptions
    VarIntUtil::WriteVarUint(exceptions.size(), serialized);

    // Exception indices (delta encoded)
    if (!exception_indices.empty()) {
        std::vector<uint64_t> indices_u64(exception_indices.begin(), exception_indices.end());
        DeltaResult delta_indices;
        DeltaEncoder::Compress(indices_u64, delta_indices);
        std::vector<uint8_t> idx_buf;
        DeltaEncoder::Serialize(delta_indices, idx_buf);

        VarIntUtil::WriteVarUint(idx_buf.size(), serialized);
        serialized.insert(serialized.end(), idx_buf.begin(), idx_buf.end());
    }

    // Exception values
    for (int32_t val : exceptions) {
        VarIntUtil::WriteVarInt(val, serialized);
    }

    return ApplyBackendCompression(serialized, output);
}

bool GQCompressor::CompressStandalone(const std::vector<int32_t>& data,
                                      CompressedField& output) {
    DPCompressor dp_comp(backend_, config_);
    return dp_comp.Compress(data, output);
}

GVCFFieldAnalysis GQCompressor::Analyze(const std::vector<int32_t>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "GQ";
    result.total_count = static_cast<uint32_t>(data.size());

    if (!data.empty()) {
        auto minmax = std::minmax_element(data.begin(), data.end()); auto min_it = minmax.first; auto max_it = minmax.second;
        result.min_value = *min_it;
        result.max_value = *max_it;
        result.avg_value = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }

    result.recommended_method = FieldCompressionMethod::MASK; // With PL prediction

    return result;
}

// ============================================================================
// PLCompressor Implementation
// ============================================================================

bool PLCompressor::Compress(const std::vector<std::vector<int32_t>>& data,
                            CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::DICTIONARY;

    // Optimized PL encoding:
    // - Standard pattern (99%+): 3 values with PL[0]=0 -> store only PL[1], PL[2]
    // - Exception pattern: use dictionary encoding

    std::vector<uint8_t> is_standard;  // Bitmask: 1=standard, 0=exception
    std::vector<int32_t> standard_pl1;  // PL[1] values for standard patterns
    std::vector<int32_t> standard_pl2;  // PL[2] values for standard patterns

    // Exception patterns dictionary
    std::unordered_map<std::string, uint32_t> exception_to_idx;
    std::vector<std::vector<int32_t>> exception_dict;
    std::vector<uint32_t> exception_indices;

    is_standard.resize((data.size() + 7) / 8, 0);

    for (size_t i = 0; i < data.size(); ++i) {
        const auto& pl = data[i];

        // Check if standard pattern: exactly 3 values and PL[0] == 0
        bool is_std = (pl.size() == 3 && pl[0] == 0);

        if (is_std) {
            // Set bit in bitmask
            is_standard[i / 8] |= (1 << (i % 8));
            standard_pl1.push_back(pl[1]);
            standard_pl2.push_back(pl[2]);
        } else {
            // Exception pattern - use dictionary
            std::string pattern;
            for (int32_t val : pl) {
                pattern += std::to_string(val) + ",";
            }

            auto it = exception_to_idx.find(pattern);
            if (it == exception_to_idx.end()) {
                uint32_t idx = static_cast<uint32_t>(exception_dict.size());
                exception_dict.push_back(pl);
                exception_to_idx[pattern] = idx;
                exception_indices.push_back(idx);
            } else {
                exception_indices.push_back(it->second);
            }
        }
    }

    // Serialize
    std::vector<uint8_t> serialized;

    // Header: version marker for new format
    serialized.push_back(2);  // Version 2: optimized PL encoding

    // Write total count
    VarIntUtil::WriteVarUint(data.size(), serialized);

    // Write standard pattern bitmask (RLE compressed)
    RLEByteResult bitmask_rle;
    RLEByteEncoder::Compress(is_standard, bitmask_rle);
    std::vector<uint8_t> bitmask_buf;
    RLEByteEncoder::Serialize(bitmask_rle, bitmask_buf);
    VarIntUtil::WriteVarUint(bitmask_buf.size(), serialized);
    serialized.insert(serialized.end(), bitmask_buf.begin(), bitmask_buf.end());

    // Write standard PL[1] values (Delta encoded for better compression)
    VarIntUtil::WriteVarUint(standard_pl1.size(), serialized);
    for (int32_t val : standard_pl1) {
        VarIntUtil::WriteVarInt(val, serialized);
    }

    // Write standard PL[2] values
    for (int32_t val : standard_pl2) {
        VarIntUtil::WriteVarInt(val, serialized);
    }

    // Write exception dictionary
    VarIntUtil::WriteVarUint(exception_dict.size(), serialized);
    for (const auto& pl : exception_dict) {
        VarIntUtil::WriteVarUint(pl.size(), serialized);
        for (int32_t val : pl) {
            VarIntUtil::WriteVarInt(val, serialized);
        }
    }

    // Write exception indices (RLE compressed)
    VarIntUtil::WriteVarUint(exception_indices.size(), serialized);
    if (!exception_indices.empty()) {
        std::vector<int32_t> indices_signed(exception_indices.begin(), exception_indices.end());
        RLEIntResult rle_indices;
        RLEIntEncoder::Compress(indices_signed, rle_indices);
        std::vector<uint8_t> idx_buf;
        RLEIntEncoder::Serialize(rle_indices, idx_buf);
        VarIntUtil::WriteVarUint(idx_buf.size(), serialized);
        serialized.insert(serialized.end(), idx_buf.begin(), idx_buf.end());
    }

    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis PLCompressor::Analyze(const std::vector<std::vector<int32_t>>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "PL";
    result.total_count = static_cast<uint32_t>(data.size());

    // Count unique patterns
    std::unordered_map<std::string, uint32_t> patterns;
    for (const auto& pl : data) {
        std::string pattern;
        for (int32_t val : pl) {
            pattern += std::to_string(val) + ",";
        }
        patterns[pattern]++;
    }

    result.unique_count = static_cast<uint32_t>(patterns.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    uint32_t max_freq = 0;
    for (const auto& kv : patterns) { const auto& count = kv.second;
        if (count > max_freq) max_freq = count;
    }
    result.dominant_ratio = static_cast<float>(max_freq) / result.total_count;

    result.recommended_method = FieldCompressionMethod::DICTIONARY;

    return result;
}

// ============================================================================
// ADCompressor Implementation
// ============================================================================

bool ADCompressor::Compress(const std::vector<std::vector<int32_t>>& data,
                            CompressedField& output) {
    // Similar to PL compression
    return PLCompressor(backend_, config_).Compress(data, output);
}

GVCFFieldAnalysis ADCompressor::Analyze(const std::vector<std::vector<int32_t>>& data) {
    GVCFFieldAnalysis result = PLCompressor::Analyze(data);
    result.field_name = "AD";
    return result;
}

// ============================================================================
// GenericFieldCompressor Implementation
// ============================================================================

bool GenericFieldCompressor::Compress(const std::vector<std::string>& data,
                                      CompressedField& output) {
    output.Clear();
    if (data.empty()) {
        return true;
    }

    output.original_count = static_cast<uint32_t>(data.size());
    output.method = FieldCompressionMethod::DICTIONARY;

    DictResult dict_result;
    if (!DictEncoder::Compress(data, dict_result)) {
        return false;
    }

    std::vector<uint8_t> serialized;
    if (!DictEncoder::Serialize(dict_result, serialized)) {
        return false;
    }

    return ApplyBackendCompression(serialized, output);
}

GVCFFieldAnalysis GenericFieldCompressor::Analyze(const std::vector<std::string>& data) {
    GVCFFieldAnalysis result;
    result.field_name = "GENERIC";
    result.total_count = static_cast<uint32_t>(data.size());

    std::unordered_map<std::string, uint32_t> freq;
    for (const auto& val : data) {
        freq[val]++;
    }

    result.unique_count = static_cast<uint32_t>(freq.size());
    result.repetition_ratio = static_cast<float>(result.unique_count) / result.total_count;

    uint32_t max_freq = 0;
    for (const auto& kv : freq) { const auto& count = kv.second;
        if (count > max_freq) max_freq = count;
    }
    result.dominant_ratio = static_cast<float>(max_freq) / result.total_count;

    result.recommended_method = FieldCompressionMethod::DICTIONARY;

    return result;
}

// ============================================================================
// GVCFBlockCompressor Implementation
// ============================================================================

GVCFBlockCompressor::GVCFBlockCompressor(std::shared_ptr<CompressionBackend> backend,
                                         const GVCFBlockConfig& config)
    : backend_(backend)
    , config_(config)
    , chrom_comp_(backend, config)
    , pos_comp_(backend, config)
    , end_comp_(backend, config)
    , id_comp_(backend, config)
    , ref_comp_(backend, config)
    , alt_comp_(backend, config)
    , qual_comp_(backend, config)
    , filter_comp_(backend, config)
    , gt_comp_(backend, config)
    , dp_comp_(backend, config)
    , min_dp_comp_(backend, config)
    , gq_comp_(backend, config)
    , pl_comp_(backend, config)
    , ad_comp_(backend, config)
    , generic_comp_(backend, config) {}

bool GVCFBlockCompressor::Compress(const GVCFBlock& input, CompressedGVCFBlock& output) {
    output.Clear();
    output.variant_count = input.variant_count;
    output.sample_count = input.sample_count;

    // Compress position fields
    if (!chrom_comp_.Compress(input.position.chrom, output.chrom)) {
        return false;
    }

    if (!pos_comp_.Compress(input.position.pos, output.pos)) {
        return false;
    }

    if (!id_comp_.Compress(input.position.id, output.id)) {
        return false;
    }

    // Compress sequence fields
    if (!ref_comp_.Compress(input.sequence.ref, output.ref)) {
        return false;
    }

    if (!alt_comp_.Compress(input.sequence.alt, output.alt)) {
        return false;
    }

    // Compress quality fields
    if (!qual_comp_.Compress(input.quality.qual, output.qual)) {
        return false;
    }

    if (!filter_comp_.Compress(input.quality.filter, output.filter)) {
        return false;
    }

    // Compress INFO/END
    if (!input.info.end.empty()) {
        output.has_end_field = true;
        if (!end_comp_.Compress(input.info.end, input.position.pos, output.info_end)) {
            return false;
        }
    }

    // Compress sample fields
    if (!gt_comp_.Compress(input.sample.gt, output.gt_mask, output.gt_patches, output.gt_phase)) {
        return false;
    }

    // Compress DP and MIN_DP together
    if (!input.sample.min_dp.empty()) {
        output.has_min_dp = true;
        if (!min_dp_comp_.Compress(input.sample.min_dp, input.sample.dp,
                                   output.dp_min_dp_diff, output.dp)) {
            return false;
        }
    } else if (!input.sample.dp.empty()) {
        if (!dp_comp_.Compress(input.sample.dp, output.dp)) {
            return false;
        }
    }

    // Compress GQ with PL prediction
    if (!input.sample.gq.empty()) {
        if (!gq_comp_.Compress(input.sample.gq, input.sample.pl, output.gq)) {
            return false;
        }
    }

    // Compress PL
    if (!input.sample.pl.empty()) {
        if (!pl_comp_.Compress(input.sample.pl, output.pl)) {
            return false;
        }
    }

    // Compress AD
    if (!input.sample.ad.empty()) {
        if (!ad_comp_.Compress(input.sample.ad, output.ad)) {
            return false;
        }
    }

    // Compress unknown fields
    for (const auto& kv : input.unknown_info) {
        const auto& name = kv.first;
        const auto& values = kv.second;
        CompressedField field;
        if (!generic_comp_.Compress(values, field)) {
            return false;
        }
        output.unknown_info[name] = std::move(field);
    }

    for (const auto& kv : input.unknown_format) {
        const auto& name = kv.first;
        const auto& values = kv.second;
        CompressedField field;
        if (!generic_comp_.Compress(values, field)) {
            return false;
        }
        output.unknown_format[name] = std::move(field);
    }

    return true;
}

GVCFBlockAnalysis GVCFBlockCompressor::Analyze(const GVCFBlock& block) {
    GVCFBlockAnalysis result;
    result.variant_count = block.variant_count;
    result.sample_count = block.sample_count;

    result.chrom_analysis = ChromCompressor::Analyze(block.position.chrom);
    result.pos_analysis = PosCompressor::Analyze(block.position.pos);
    result.id_analysis = IdCompressor::Analyze(block.position.id);
    result.ref_analysis = RefCompressor::Analyze(block.sequence.ref);
    result.alt_analysis = AltCompressor::Analyze(block.sequence.alt);
    result.gt_analysis = GTCompressor::Analyze(block.sample.gt);
    result.dp_analysis = DPCompressor::Analyze(block.sample.dp);
    result.gq_analysis = GQCompressor::Analyze(block.sample.gq);

    if (!block.sample.min_dp.empty()) {
        result.min_dp_analysis = MinDPCompressor::Analyze(block.sample.min_dp, block.sample.dp);
    }

    if (!block.info.end.empty()) {
        result.end_analysis = EndCompressor::Analyze(block.info.end);
    }

    // Check if typical gVCF
    result.is_typical_gvcf = (result.gt_analysis.dominant_ratio >= 0.7f) &&
                             (result.alt_analysis.dominant_ratio >= 0.7f);

    return result;
}

} // namespace gvcf
