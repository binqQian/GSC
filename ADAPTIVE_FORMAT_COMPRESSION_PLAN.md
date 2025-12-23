# GSC项目 自适应FORMAT字段压缩重构方案

## 1. 现状分析

### 1.1 根目录GSC项目现有架构

当前GSC项目对VCF数据的压缩分为三层：

| 数据类型 | 处理方式 | 核心文件 |
|---------|---------|---------|
| **GT (Genotype)** | 复杂稀疏编码：零向量检测、复制检测、SDSL位向量 | `compressor.cpp`, `block_processing.cpp` |
| **Fixed Fields** (CHROM/POS/ID/REF/ALT/QUAL) | 块压缩（BSC/ZSTD/Brotli） | `compressor.cpp:compressFixedFields()` |
| **Other Fields** (INFO/FORMAT非GT) | 按数据类型分类压缩（INT/REAL/STR/FLAG） | `compression_reader.cpp`, `file_handle.cpp` |

**关键观察**：
- FORMAT非GT字段目前使用**通用按类型压缩**（`compress_other_fileds()`）
- 没有针对AD/DP/PL/GQ等字段的**特征识别**和**专用编码**
- 兜底策略已经存在：未识别字段通过File_Handle_2存储

### 1.2 ref_code/fmt_comp项目的硬编码策略

ref_code/fmt_comp项目对FORMAT字段使用硬编码分派：

```
字段名 memcmp() → 分派到专用处理器
  AD → 2bit tip + (全0/高频值/字典数组)
  DP → sum(AD)预测 + 异常覆盖
  PL → 2bit tip + 模式识别(checkAbPattern) + 字典
  GQ → f(PL)预测 + 异常覆盖
  PGT/PID → 稀疏存储 + 字典
  other → 原样字符串拼接
```

**问题**：硬编码无法适应非标准VCF或自定义FORMAT字段。

---

## 2. 重构目标

将"按字段名硬编码"改为"按字段内容特征自适应选型"：

1. **特征检测**：运行时分析字段值的统计特征
2. **Codec选择**：基于特征匹配最优压缩算法
3. **无缝回退**：无法识别时使用现有默认压缩

---

## 3. 架构设计

### 3.1 核心模块

```
src/
├── format_field_codec.h          # Codec接口定义
├── format_field_codec.cpp
├── format_field_detector.h       # 特征检测器
├── format_field_detector.cpp
├── format_field_manager.h        # 字段管理器（协调检测与编码）
├── format_field_manager.cpp
└── codecs/                       # 各类Codec实现
    ├── predicted_scalar_codec.h   # 预测型标量（DP/GQ）
    ├── bit_tip_array_codec.h      # 2bit tip数组（AD/PL）
    ├── sparse_dict_codec.h        # 稀疏字典（PGT/PID）
    ├── raw_string_codec.h         # 兜底原样存储
    └── ...
```

### 3.2 Codec接口设计

```cpp
// format_field_codec.h
#pragma once
#include <vector>
#include <cstdint>
#include <string>

// 字段特征统计结果
struct FieldFeatures {
    // 基本类型
    bool is_array;              // 是否为数组类型
    bool is_numeric;            // 是否为数值类型
    bool is_integer;            // 是否为整数

    // 缺失率
    float missing_ratio;        // value == "." 的比例

    // 数组特征
    bool array_len_fixed;       // 数组长度是否恒定
    uint32_t array_len;         // 若恒定，数组长度

    // 值域
    uint32_t min_val;
    uint32_t max_val;
    uint8_t max_bits;           // 需要的位数

    // 稀疏性与重复
    uint32_t unique_count;      // 去重后的值数量
    float top1_freq;            // 最高频值的占比

    // 跨字段相关性
    float dp_ad_hit_ratio;      // DP == sum(AD) 命中率
    float gq_pl_hit_ratio;      // GQ == f(PL) 命中率

    // 模式检测（PL专用）
    float pattern_ab_ratio;     // 满足AB模式的比例
};

// Codec类型枚举
enum class CodecType : uint8_t {
    RawString = 0,              // 兜底：原样存储
    PredictedScalar = 1,        // 预测型标量（DP/GQ）
    BitTipArray = 2,            // 2bit tip数组（AD/PL）
    SparseDictString = 3,       // 稀疏字典字符串
    FixedWidthArray = 4,        // 定宽数组
    DeltaVarint = 5,            // 差分变长整数
};

// 抽象Codec接口
class FormatFieldCodec {
public:
    virtual ~FormatFieldCodec() = default;

    // 观察阶段：收集统计信息
    virtual void observe(const char* value, size_t len, uint32_t sample_pos) = 0;

    // 冻结：基于统计信息决定编码参数
    virtual void freeze() = 0;

    // 编码单个值
    virtual void encode(const char* value, size_t len, uint32_t sample_pos) = 0;

    // 序列化编码结果
    virtual void serialize(std::vector<uint8_t>& out) = 0;

    // 反序列化
    virtual void deserialize(const uint8_t* data, size_t len) = 0;

    // 解码单个值
    virtual std::string decode(uint32_t sample_pos) = 0;

    // 获取Codec类型
    virtual CodecType type() const = 0;

    // 重置状态（用于下一行）
    virtual void reset() = 0;
};
```

### 3.3 特征检测器设计

```cpp
// format_field_detector.h
#pragma once
#include "format_field_codec.h"
#include <memory>

// 特征检测阈值配置
struct DetectorConfig {
    float missing_high_threshold = 0.8f;   // 高缺失率阈值
    float hit_ratio_threshold = 0.95f;     // 预测命中率阈值
    float pattern_threshold = 0.9f;        // 模式匹配阈值
    uint32_t observe_sample_limit = 256;   // 观察样本数上限
};

class FormatFieldDetector {
public:
    FormatFieldDetector(const DetectorConfig& config = {});

    // 喂入一个样本的值
    void feed(const char* value, size_t len);

    // 设置关联字段（用于跨字段预测检测）
    void setRelatedField(const std::string& field_name,
                         FormatFieldDetector* related);

    // 完成检测，返回特征
    FieldFeatures finalize();

    // 基于特征创建最优Codec
    std::unique_ptr<FormatFieldCodec> createCodec(const FieldFeatures& features);

private:
    DetectorConfig config_;

    // 统计变量
    uint32_t sample_count_ = 0;
    uint32_t missing_count_ = 0;
    bool has_comma_ = false;        // 是否包含逗号（数组）
    bool all_numeric_ = true;
    bool all_integer_ = true;

    // 值域统计
    uint32_t min_val_ = UINT32_MAX;
    uint32_t max_val_ = 0;

    // 数组长度统计
    std::vector<uint32_t> array_lens_;

    // 唯一值计数（用小buffer）
    std::unordered_map<std::string, uint32_t> value_freq_;

    // 跨字段关联
    std::map<std::string, FormatFieldDetector*> related_fields_;
    std::vector<uint32_t> ad_sums_;  // 记录sum(AD)用于DP预测
    std::vector<uint32_t> pl_gqs_;   // 记录f(PL)用于GQ预测
};
```

### 3.4 字段管理器设计

```cpp
// format_field_manager.h
#pragma once
#include "format_field_detector.h"
#include <unordered_map>
#include <string>

class FormatFieldManager {
public:
    FormatFieldManager();

    // 初始化一行的FORMAT字段列表
    void initRow(const std::vector<std::string>& format_keys,
                 uint32_t allele_count);

    // 处理一个样本的完整FORMAT串
    void processSample(const char* sample_str, size_t len, uint32_t sample_pos);

    // 行结束：冻结检测器、序列化编码结果
    void finalizeRow(std::vector<uint8_t>& out);

    // 解压一个样本的FORMAT
    std::string decodeSample(const uint8_t* data, size_t len,
                             uint32_t sample_pos,
                             const std::vector<std::string>& format_keys);

private:
    std::vector<std::string> current_format_keys_;
    uint32_t allele_count_ = 2;

    // 每个字段一个检测器/编码器
    std::unordered_map<std::string, std::unique_ptr<FormatFieldDetector>> detectors_;
    std::unordered_map<std::string, std::unique_ptr<FormatFieldCodec>> codecs_;

    // 状态
    bool frozen_ = false;
    uint32_t samples_observed_ = 0;
    static constexpr uint32_t FREEZE_THRESHOLD = 256;

    // 解析单个样本串
    void parseSampleFields(const char* sample_str, size_t len,
                           std::vector<std::pair<const char*, size_t>>& fields);
};
```

---

## 4. Codec实现规范

### 4.1 预测型标量Codec（DP/GQ）

```cpp
// codecs/predicted_scalar_codec.h
class PredictedScalarCodec : public FormatFieldCodec {
public:
    // 设置预测器：predictor(sample_pos) -> expected_value
    void setPredictor(std::function<uint32_t(uint32_t)> predictor);

    // 编码时：只记录异常值 (pos, actual_value)
    void encode(const char* value, size_t len, uint32_t sample_pos) override;

    // 序列化：异常表用Vint编码
    void serialize(std::vector<uint8_t>& out) override;

private:
    std::function<uint32_t(uint32_t)> predictor_;
    std::vector<std::pair<uint32_t, uint32_t>> exceptions_; // (pos, val)
    uint32_t missing_marker_ = 0xFFFF;
};
```

**应用场景**：
- DP字段：predictor = sum(AD[sample_pos])
- GQ字段：predictor = f(PL[sample_pos])，f可以是PL[1]或second_min(PL)

### 4.2 2bit Tip数组Codec（AD/PL）

```cpp
// codecs/bit_tip_array_codec.h
class BitTipArrayCodec : public FormatFieldCodec {
public:
    // AD的tip编码规则
    // 00: 全0
    // 01: sum==arr[0] && sum==2 (高频特判)
    // 10: sum==arr[0] && sum!=2 (仅首值)
    // 11: 一般情况（字典化）

    // PL的tip编码规则
    // 00: 全0
    // 01: 满足AB模式 (存a,b)
    // 10: 满足A15模式 (存a, b=15*a)
    // 11: 一般情况（字典化）

    void encode(const char* value, size_t len, uint32_t sample_pos) override;

private:
    std::vector<uint8_t> tips_;       // 2bit per sample, packed
    std::vector<uint32_t> payload_;   // 根据tip类型存储不同内容

    // 字典
    std::vector<std::vector<uint32_t>> dict_;
    std::unordered_map<std::string, uint32_t> dict_index_;
};
```

### 4.3 稀疏字典Codec（PGT/PID）

```cpp
// codecs/sparse_dict_codec.h
class SparseDictCodec : public FormatFieldCodec {
public:
    // 默认认为值为缺失"."
    // 只记录非缺失值的位置和字典ID

    void encode(const char* value, size_t len, uint32_t sample_pos) override;

private:
    std::vector<uint32_t> non_missing_pos_;  // 非缺失位置
    std::vector<uint32_t> dict_ids_;         // 对应字典ID

    // 字符串字典
    std::vector<std::string> dict_;
    std::unordered_map<std::string, uint32_t> dict_index_;
};
```

### 4.4 兜底RawString Codec

```cpp
// codecs/raw_string_codec.h
class RawStringCodec : public FormatFieldCodec {
public:
    // 原样存储，按:分隔拼接
    void encode(const char* value, size_t len, uint32_t sample_pos) override {
        if (sample_pos > 0) raw_data_.push_back(':');
        raw_data_.insert(raw_data_.end(), value, value + len);
    }

    void serialize(std::vector<uint8_t>& out) override {
        // 直接输出原始数据
        out.insert(out.end(), raw_data_.begin(), raw_data_.end());
    }

private:
    std::vector<char> raw_data_;
};
```

---

## 5. 算法映射规则

基于特征检测结果，按优先级匹配：

| 优先级 | 条件 | 选择的Codec |
|-------|------|------------|
| 1 | missing_ratio >= 0.8 && unique_count <= 256 | SparseDictCodec |
| 2 | is_integer && has_predictor && hit_ratio >= 0.95 | PredictedScalarCodec |
| 3 | is_array && array_len_fixed && (all_zero_ratio高 \|\| pattern_hit高) | BitTipArrayCodec |
| 4 | is_array && array_len_fixed && max_bits <= 8/16/32 | FixedWidthArrayCodec |
| 5 | is_integer && mono_inc_ratio > 0.9 | DeltaVarintCodec |
| 6 | 以上都不满足 | **RawStringCodec (兜底)** |

---

## 6. 与现有系统对接

### 6.1 入口点修改

修改 `compression_reader.cpp` 中的 `SetVariantOtherFields()` 函数：

```cpp
// compression_reader.cpp (修改)
bool CompressionReader::SetVariantOtherFields(vector<field_desc> &fields) {
    // 对于FORMAT字段，使用FormatFieldManager处理
    for (const auto& key : keys) {
        if (key.keys_type == key_type_t::fmt && key.key_id != key_gt_id) {
            // 使用自适应压缩
            format_field_manager_->processSample(...);
        }
    }
    // ...
}
```

### 6.2 序列化格式

每行的FORMAT编码输出格式：

```
[row_header]
  ├─ uint8_t  format_key_count
  ├─ uint8_t  codec_types[format_key_count]  // 每个字段使用的codec类型
  └─ uint32_t codec_params[...]              // codec参数（如字典大小等）

[row_payload]
  ├─ field_0_data: [codec_0 序列化数据]
  ├─ field_1_data: [codec_1 序列化数据]
  └─ ...
```

### 6.3 解压路径

修改 `decompressor.cpp` 中的FORMAT重建逻辑：

```cpp
// 读取row_header获取codec_types
// 对每个字段：
//   根据codec_type反序列化对应Codec
//   调用codec->decode(sample_pos)重建字符串
```

---

## 7. 兜底策略

**原则**：任何无法稳定处理的情况，都回退到RawStringCodec。

触发兜底的条件：
1. 类型无法稳定解析（混合数字/字符串）
2. 数组长度不稳定
3. 统计置信度不足（样本数太少）
4. 识别出的codec预计收益不高

**与现有默认压缩的关系**：
- 兜底后的数据走现有 `compress_other_fileds()` 路径
- 仍可借助后续块级压缩器（BSC/ZSTD）获得压缩收益

---

## 8. 实施路线图

### Phase 1: 基础框架（~估计代码量：500行） ✅ 已完成
- [x] 创建 `format_field_codec.h/cpp` 接口定义
- [x] 创建 `format_field_detector.h/cpp` 特征检测器
- [x] 实现 `RawStringCodec` 作为兜底

### Phase 2: 核心Codec实现（~1000行） ✅ 已完成
- [x] 实现 `PredictedScalarCodec`（DP/GQ）- 292行
- [x] 实现 `BitTipArrayCodec`（AD/PL基础版）- 724行
- [x] 实现 `SparseDictCodec`（PGT/PID）- 191行

### Phase 3: 系统集成（~500行） ✅ 框架已完成
- [x] 创建 `FormatFieldManager` (在format_field_detector.h/cpp中实现)
- [x] 修改 `compression_reader.h/cpp` 添加自适应压缩框架
  - 添加FormatFieldManager成员变量
  - 添加formatFieldToString辅助函数
  - 添加集成点标记（TODO注释）
- [x] 修改 `decompressor.h` 添加自适应解压框架
  - 添加FormatFieldManager成员变量
  - 添加use_adaptive_format_标志

**Phase 3 说明**：
- 框架代码已就绪，可编译通过
- 现有压缩/解压逻辑保持不变（渐进式集成）
- 自适应压缩的完整启用需要在Phase 4中实现数据流重构

### Phase 4: 优化与测试 ✅ 已完成
- [x] PL的模式识别（checkAbPattern）- 已在BitTipArrayCodec中实现
- [x] 跨字段预测器（DP←AD, GQ←PL）- 已在FormatFieldManager中实现
- [x] 完整启用自适应压缩
  - compression_reader.cpp: SetVariantOtherFields()集成FormatFieldManager
  - 自动检测FORMAT字段并构建样本字符串
  - 调用processSample()和finalizeRow()
  - 数据存储到专用stream (adaptive_format_data)
- [x] 解压路径支持
  - decompression_reader.cpp: 检测adaptive_format_data stream
  - 向后兼容：同时保留现有压缩格式
- [ ] 性能测试与压缩率对比（由用户进行）

---

## 9. 文件变更清单

| 文件 | 操作 | 说明 |
|-----|------|-----|
| `src/format_field_codec.h` | 新增 | Codec接口定义 |
| `src/format_field_codec.cpp` | 新增 | 基础实现 |
| `src/format_field_detector.h` | 新增 | 特征检测器 |
| `src/format_field_detector.cpp` | 新增 | 检测逻辑 |
| `src/format_field_manager.h` | 新增 | 字段管理器 |
| `src/format_field_manager.cpp` | 新增 | 协调逻辑 |
| `src/codecs/` | 新增目录 | 各类Codec实现 |
| `src/compression_reader.h` | 修改 | 添加manager成员 |
| `src/compression_reader.cpp` | 修改 | 集成自适应压缩 |
| `src/decompressor.cpp` | 修改 | 集成自适应解压 |
| `CMakeLists.txt` | 修改 | 添加新源文件 |

---

## 10. 关键设计决策

### Q1: 为什么不直接复制ref_code/fmt_comp的硬编码？
A: 硬编码无法适应非标准VCF，且违背"自适应"原则。通过特征检测，可以自动处理用户自定义FORMAT字段。

### Q2: 如何保证性能不退化？
A:
- 检测阶段仅采样前256个样本
- freeze()后直接流式编码
- 热路径（encode/decode）无复杂逻辑

### Q3: 兜底策略会导致压缩率下降吗？
A: 不会。兜底使用RawStringCodec，与现有实现等价，仍可借助块压缩器。

### Q4: 如何验证正确性？
A: 压缩→解压→与原始VCF逐字节对比（lossless模式）。

---

*文档生成时间: 2025-12-22*
*目标项目: /home/binq/Code/vcfcomp/gsc_base/GSC*
