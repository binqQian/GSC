/**
 * @file gvcf_field_decompress.h
 * @brief gVCF field-level decompression
 *
 * Provides decompression for each gVCF field type.
 * Reference: ref_code/gvcf_refCode/decoding/
 */
#pragma once

#include "gvcf_encoding.h"
#include "gvcf_block.h"
#include "gvcf_field_compress.h"

namespace gvcf {

// 所有字段解压器的基类：先处理后端压缩标记，再交给具体字段恢复。

class FieldDecompressor {
protected:
    std::shared_ptr<CompressionBackend> backend_;

public:
    explicit FieldDecompressor(std::shared_ptr<CompressionBackend> backend)
        : backend_(backend) {}

    virtual ~FieldDecompressor() = default;

    // 读取字段头部的“是否压缩”标记，并在需要时调用后端解压。
    bool ApplyBackendDecompression(const CompressedField& input,
                                   std::vector<uint8_t>& output);
};

// ============================================================================
// CHROM Field Decompressor
// ============================================================================

class ChromDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<std::string>& output);
};

// ============================================================================
// POS Field Decompressor
// ============================================================================

class PosDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<uint64_t>& output);
};

// ============================================================================
// INFO/END Field Decompressor
// ============================================================================

class EndDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    // Decompress END relative to POS
    bool Decompress(const CompressedField& input,
                   const std::vector<uint64_t>& pos_values,
                   std::vector<int64_t>& output);

    // Decompress END standalone
    bool DecompressStandalone(const CompressedField& input,
                             std::vector<int64_t>& output);
};

// ============================================================================
// ID Field Decompressor
// ============================================================================

class IdDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<std::string>& output);
};

// ============================================================================
// REF Field Decompressor
// ============================================================================

class RefDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<std::string>& output);
};

// ============================================================================
// ALT Field Decompressor
// ============================================================================

class AltDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input,
                   std::vector<std::vector<std::string>>& output);
};

// ============================================================================
// QUAL Field Decompressor
// ============================================================================

class QualDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<float>& output);
};

// ============================================================================
// FILTER Field Decompressor
// ============================================================================

class FilterDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<std::string>& output);
};

// GT 解压器：把 mask / patches / phase 三段数据拼回完整基因型。

class GTDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    // Decompress from three compressed fields
    bool Decompress(const CompressedField& mask_input,
                   const CompressedField& patches_input,
                   const CompressedField& phase_input,
                   std::vector<GenotypeData>& output);

    // Decompress from single compressed field
    bool DecompressSingle(const CompressedField& input,
                         std::vector<GenotypeData>& output);

    // Decompress to string representation
    bool DecompressStrings(const CompressedField& mask_input,
                          const CompressedField& patches_input,
                          const CompressedField& phase_input,
                          std::vector<std::string>& output);
};

// ============================================================================
// DP Field Decompressor
// ============================================================================

class DPDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<int32_t>& output);
};

// MIN_DP 解压器：优先借助 DP 恢复，避免重复存整列数据。

class MinDPDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    // Decompress MIN_DP from differences and DP
    bool Decompress(const CompressedField& diff_input,
                   const std::vector<int32_t>& dp,
                   std::vector<int32_t>& output);

    // Decompress MIN_DP standalone
    bool DecompressStandalone(const CompressedField& input,
                             std::vector<int32_t>& output);
};

// GQ 解压器：优先从 PL 推断，再回填异常值。

class GQDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    // Decompress GQ with PL prediction
    bool Decompress(const CompressedField& input,
                   const std::vector<std::vector<int32_t>>& pl,
                   std::vector<int32_t>& output);

    // Decompress GQ standalone
    bool DecompressStandalone(const CompressedField& input,
                             std::vector<int32_t>& output);
};

// PL / AD 解压器负责恢复 sample 级数组字段。

class PLDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input,
                   std::vector<std::vector<int32_t>>& output);
};

// ============================================================================
// AD Field Decompressor
// ============================================================================

class ADDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input,
                   std::vector<std::vector<int32_t>>& output);
};

// ============================================================================
// Generic Field Decompressor
// ============================================================================

class GenericFieldDecompressor : public FieldDecompressor {
public:
    using FieldDecompressor::FieldDecompressor;

    bool Decompress(const CompressedField& input, std::vector<std::string>& output);
};

// Block 解压调度器：按字段依赖顺序恢复一个完整的原始块。

class GVCFBlockDecompressor {
public:
    explicit GVCFBlockDecompressor(std::shared_ptr<CompressionBackend> backend);

    // 解整个块，适合正常解压路径。
    bool Decompress(const CompressedGVCFBlock& input, GVCFBlock& output);

    // 只解部分字段，适合范围查询或轻量检查。
    bool DecompressPositionFields(const CompressedGVCFBlock& input, GVCFBlock& output);
    bool DecompressSampleFields(const CompressedGVCFBlock& input, GVCFBlock& output);

    // Validate decompression (compare with original)
    bool Validate(const GVCFBlock& original, const GVCFBlock& decompressed);

private:
    std::shared_ptr<CompressionBackend> backend_;

    ChromDecompressor chrom_decomp_;
    PosDecompressor pos_decomp_;
    EndDecompressor end_decomp_;
    IdDecompressor id_decomp_;
    RefDecompressor ref_decomp_;
    AltDecompressor alt_decomp_;
    QualDecompressor qual_decomp_;
    FilterDecompressor filter_decomp_;
    GTDecompressor gt_decomp_;
    DPDecompressor dp_decomp_;
    MinDPDecompressor min_dp_decomp_;
    GQDecompressor gq_decomp_;
    PLDecompressor pl_decomp_;
    ADDecompressor ad_decomp_;
    GenericFieldDecompressor generic_decomp_;
};

} // namespace gvcf
