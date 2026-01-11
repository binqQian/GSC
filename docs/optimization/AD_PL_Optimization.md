# AD/PL 字段压缩优化方案

> 本文档分析 VCF FORMAT/AD 和 FORMAT/PL 字段的数据特征、当前编码现状、BSC vs Brotli 压缩差异原因，并提出进一步优化方案。

---

## 1. 字段特征分析

### 1.1 测试数据集

| 数据集               | 路径                                                | 样本数 | 位点数 | 特点                   |
| -------------------- | --------------------------------------------------- | ------ | ------ | ---------------------- |
| **LC** (数据集A)     | `/home/binq/data/1000GPi/lc_bams.first50000.vcf.gz` | 2,560  | 46,601 | 低覆盖度 WGS，模式多样 |
| **Subset** (数据集B) | `/home/binq/data/output_subset_5000_20000.vcf.gz`   | 5,000  | 20,000 | 高覆盖度子集，模式规整 |

### 1.2 AD 字段特征

#### 1.2.1 数值分布

| 指标       | LC     | Subset |
| ---------- | ------ | ------ |
| 零值比例   | 41.69% | 56.47% |
| 均值       | 7.91   | 4.47   |
| 标准差     | 20.78  | 7.22   |
| 二等位比例 | 94.31% | 87.18% |

#### 1.2.2 当前 Tip 编码分布 (codec_id=4)

| Tip | 含义            | LC         | Subset     | 存储        |
| --- | --------------- | ---------- | ---------- | ----------- |
| 00  | 全零 sum(AD)==0 | 0.00%      | 6.95%      | 无          |
| 01  | ref-only (n,0)  | **76.54%** | **86.62%** | delta_sum   |
| 10  | alt-only (0,n)  | 1.55%      | 0.17%      | delta_sum   |
| 11  | 一般情况        | **19.34%** | **0.72%**  | tag+payload |

**关键差异**：
- LC 有 19.34% 需要复杂编码（Tip 11），Subset 仅 0.72%
- LC 的 AD 模式数：159,388；Subset 仅 7,681（差 20 倍）
- LC 杂合比例 20.24%；Subset 仅 0.84%

#### 1.2.3 高频 AD 模式

```
LC Top 5:     (6,0): 7.75%  (5,0): 7.75%  (7,0): 7.27%  (4,0): 7.09%  (8,0): 6.47%
Subset Top 5: (4,0): 5.45%  (5,0): 5.41%  (3,0): 5.35%  (6,0): 5.24%  (0,0): 5.24%
```

绝大多数是 `(n,0)` 形式，即 **纯合位点只有 REF 深度**。

### 1.3 PL 字段特征

#### 1.3.1 数值分布

| 指标     | LC        | Subset    |
| -------- | --------- | --------- |
| 值域     | [0, 5369] | [0, 7614] |
| 均值     | 131.62    | 132.58    |
| 中位数   | 24.0      | 21.0      |
| 零值比例 | 31.33%    | 36.15%    |

#### 1.3.2 当前 Type 分布 (codec_id=4)

| Type       | 条件         | LC         | Subset     | 存储               |
| ---------- | ------------ | ---------- | ---------- | ------------------ |
| 0 (Tip 00) | 全零 [0,0,0] | 0.00%      | 7.54%      | 无                 |
| 2 (Tip 01) | b=15*a       | 0.14%      | **29.82%** | delta_a            |
| 1 (Tip 10) | ratio≈11     | **90.10%** | 61.93%     | delta_a + residual |
| 3 (Tip 11) | 兜底         | **9.76%**  | **0.71%**  | tag+payload        |

**关键差异**：
- LC 的 Type 3 占 9.76%，几乎全部是 `PL[0]!=0` 导致（已通过 min-pos 归一化解决）
- Subset 有 29.82% 满足 `b=15*a`（GATK 特征），只需存储 `a`

#### 1.3.3 Type 1 的 b/a 比值分布

| 比值         | LC        | Subset    |
| ------------ | --------- | --------- |
| ratio=10     | 11.9%     | 14.8%     |
| **ratio=11** | **37.0%** | **25.7%** |
| **ratio=12** | **44.6%** | **25.9%** |
| ratio=13     | 1.6%      | 6.0%      |

**发现**：74-93% 的 b/a 集中在 10-12，使用 `residual = b - 11*a` 可将大部分 b 压缩到小残差。

### 1.4 跨字段关联

| 关系                 | LC     | Subset      |
| -------------------- | ------ | ----------- |
| DP == sum(AD)        | 86.04% | **100.00%** |
| GQ == second_min(PL) | 94.13% | 98.99%      |

**结论**：异常表策略对 Subset 完美适用，对 LC 有 14% 异常需存储。

---

## 2. BSC vs Brotli 压缩差异分析

### 2.1 实验现象

| 数据集           | BSC 压缩率 | Brotli 压缩率 | 更优算法 |
| ---------------- | ---------- | ------------- | -------- |
| LC (数据集A)     | **更好**   | 较差          | BSC      |
| Subset (数据集B) | 较差       | **更好**      | Brotli   |

**本仓库 toy 数据快速复现（仅用于 sanity check）**：

```bash
GSC_LOG_LEVEL=debug ./build/gsc compress --in toy/final_subset.vcf.gz --out tmp/adpl_review_bsc.gsc --compressor bsc
GSC_LOG_LEVEL=debug ./build/gsc compress --in toy/final_subset.vcf.gz --out tmp/adpl_review_brotli.gsc --compressor brotli
```

在 `toy/final_subset.vcf.gz`（500 samples, 1300 variants）上：Brotli 总体更小（约 923KB vs 1008KB），且 PL 仍占主要比例（约 71-72%）。

### 2.2 算法原理对比

#### 2.2.1 BSC (Block-Sorting Compressor)

基于 **BWT (Burrows-Wheeler Transform)** + **MTF** + 熵编码：

```
原理：
1. BWT 对输入块进行全局字符排序，将相似上下文聚集
2. MTF (Move-To-Front) 将重复字符转换为小整数
3. 熵编码压缩最终流

优势：
- 擅长捕获 **长范围相关性**（数千字节跨度）
- 对 **高度重复但分散** 的模式效果好
- 全局视角，能发现跨越大距离的相似性

劣势：
- 块排序开销较大
- 对已经很规整的数据，开销可能超过收益
```

#### 2.2.2 Brotli

基于 **LZ77** + **Huffman/ANS** 编码：

```
原理：
1. LZ77 在滑动窗口内查找重复字符串
2. 用 (距离, 长度) 对替换重复
3. Huffman/ANS 编码最终符号流

优势：
- 对 **短范围重复** 效果好（窗口内匹配）
- 开销低，适合高度规整的数据
- 通用性强

劣势：
- 窗口大小有限（默认 32KB-16MB）
- 对长范围相关性捕获不如 BWT
```

### 2.3 差异原因深度分析

#### 2.3.1 LC 数据集为何 BSC 更优

**数据特征**：
- **高局部熵**：19.34% 的 AD 需要复杂编码，9.76% 的 PL 无法模式匹配
- **模式数量多**：AD 有 159,388 种唯一模式
- **长范围相关性强**：相邻位点 sum(AD) 相关性 0.990

**BSC 优势体现**：

```
样本1: [8,0] [8,0] [9,0] [8,1] [7,0] ...  (位点1-5)
样本2: [7,0] [7,0] [8,0] [7,1] [6,0] ...
样本3: [9,0] [9,0] [10,0] [9,1] [8,0] ...
       ↓
BWT 排序后：
[6,0] [7,0] [7,0] [7,0] [7,1] [8,0] [8,0] [8,0] [8,1] [9,0] ...
       ↓
MTF 转换：相同值变为 0，新值移到前端
0, 0, 0, 1, 0, 0, 0, 1, 0, ...
       ↓
熵编码：大量 0/1，压缩率高
```

即使模式分散在不同位点，BWT 能将它们聚集在一起。

**Brotli 劣势**：

```
原始流：[8,0] [7,0] [9,0] ... [8,0] [7,0] [9,0] ...
        ↑                     ↑
        如果距离 > 窗口大小，无法匹配

LC 数据的 2560 样本 × 多种模式 → 相同模式可能相距很远
LZ77 窗口无法覆盖
```

#### 2.3.2 Subset 数据集为何 Brotli 更优

**数据特征**：
- **极低熵**：仅 0.72% 的 AD 需要复杂编码，0.71% 的 PL 无法模式匹配
- **模式数量少**：AD 仅 7,681 种（LC 的 1/20）
- **高度规整**：99.28% 纯合，29.82% 满足 b=15a

**Brotli 优势体现**：

```
编码后的字节流：
[tip=01][delta=0] [tip=01][delta=0] [tip=01][delta=1] [tip=01][delta=0] ...
                  ↑ 大量连续重复 ↑

LZ77 滑动窗口：
"[tip=01][delta=0]" 重复 → (距离=N, 长度=K)

极短的模式就能达到很高压缩率
```

**BSC 劣势**：

```
Subset 数据已经很规整，模式少
BWT 的块排序、MTF 转换等开销
反而成为负担，压缩率不如直接 LZ77
```

### 2.4 理论分析：熵与冗余

#### 2.4.1 数据熵估算

| 数据集 | AD 模式熵 (log2)              | PL 模式熵 (log2)                | 总体熵 |
| ------ | ----------------------------- | ------------------------------- | ------ |
| LC     | log2(159,388) ≈ **17.3 bits** | log2(1,058,654) ≈ **20.0 bits** | **高** |
| Subset | log2(7,681) ≈ **12.9 bits**   | log2(316,746) ≈ **18.3 bits**   | **中** |

#### 2.4.2 冗余类型

| 冗余类型     | LC     | Subset   | 适合算法      |
| ------------ | ------ | -------- | ------------- |
| 短范围重复   | 中等   | **极高** | LZ77 (Brotli) |
| 长范围相关   | **高** | 中等     | BWT (BSC)     |
| 符号频率偏斜 | 高     | 极高     | Huffman/ANS   |

### 2.5 结论

| 数据特征          | 推荐后端   | 原因                       |
| ----------------- | ---------- | -------------------------- |
| 高熵 + 长范围相关 | **BSC**    | BWT 全局排序捕获分散相似性 |
| 低熵 + 短范围重复 | **Brotli** | LZ77 窗口足够，开销低      |
| 混合特征          | 需实测     | 取决于具体数据分布         |

**建议**：实现自适应后端选择，根据预处理后数据的熵估算选择 BSC 或 Brotli。

---

## 3. 当前编码现状评估

### 3.1 AD 编码 (codec_id=4)

**当前策略（与代码实现对齐）**：

```cpp
// 代码位置（以当前 mainline 为准）:
// - 编码: src/compressor.cpp（FMT AD special codec_id=4）
// - 解码: src/decompression_reader.cpp（FMT AD codec_id=4）
//
// stream layout（codec_id=4, codec marker 后）:
//   [block_size_vint]              // getenv("GSC_AD_BASELINE_BLOCK"), 默认 16384
//   for each record:
//     if (rec_idx % block_size == 0): [baseline_sum[s] vint] * n_samples
//     [tip_bytes][payload_len_vint][payload...]
//
// Tips（2-bit）:
//   00: all_zero                      -> no payload; baseline_sum=0
//   01: ref_only (AD[1..]=0)          -> zigzag(delta_sum)               (sum vs baseline_sum)
//   10: alt_only (biallelic, AD[0]=0) -> zigzag(delta_sum)               (sum vs baseline_sum)
//   11: tag-based:
//         tag=0: raw u32 bits (vint)
//         tag=1: biallelic both non-zero -> [zigzag(delta_sum)][alt]      (lossless)
//         tag>=2: dict_id = tag - 2      -> AD dictionary item (lossless)
```

**效果评估**：

| 指标              | LC     | Subset | 评价             |
| ----------------- | ------ | ------ | ---------------- |
| Tip 00+01+10 覆盖 | 78.09% | 93.74% | 良好             |
| Tip 11 占比       | 19.34% | 0.72%  | LC 较高          |
| 压缩后占比        | 25.13% | 24.70% | **仍是主要开销** |

**问题**：
1. Tip 11 并不等价于“完整存储 AD 数组”。在当前实现中，二等位杂合（ref>0 && alt>0）会走 `tag=1` 的 `sum+alt` 快速路径，依然是无损且通常比 raw 更轻量。
2. 仍然可能偏重的部分主要来自：多等位（per_sample>2）、含缺失/特殊值、或无法装入小字典 item 的样本（回退到 tag=0 raw）。

#### 3.1.1 已实现优化：AD codec_id=5（Δsum 子流拆分）

**目标**：把最常见的 ref-only / alt-only 的 `zigzag(Δsum)` 从“混合 payload”里拆出来，降低字节交错熵，让后端压缩更容易利用重复性。

```cpp
// 代码位置:
// - 编码: src/compressor.cpp（FMT AD special codec_id=5）
// - 解码: src/decompression_reader.cpp（FMT AD codec_id=5）
//
// 启用:
//   GSC_AD_CODEC_ID=5
//
// stream layout（codec_id=5, codec marker 后）:
//   [block_size_vint]
//   for each record:
//     if (rec_idx % block_size == 0): [baseline_sum[s] vint] * n_samples
//     [tip_bytes][payload_len_vint][payload...]
//
// payload（record）:
//   [delta_len_vint][delta_stream...][fallback_stream...]
//   - delta_stream: 仅包含 Tip 01 / Tip 10 的 zigzag(Δsum)
//   - fallback_stream: Tip 11 的 tag+payload（raw / sum+alt / dict），以及含缺失/特殊值的 raw
```

### 3.2 PL 编码 (codec_id=4)

**当前策略（与代码实现对齐）**：

```cpp
// 代码位置（以当前 mainline 为准）:
// - 编码: src/compressor.cpp（FMT PL special codec_id=4）
// - 解码: src/decompression_reader.cpp（FMT PL codec_id=4）
//
// stream layout（codec_id=4, codec marker 后）:
//   [block_size_vint]              // getenv("GSC_PL_BLOCK_SIZE"), 默认 16384
//   for each record:
//     if (rec_idx % block_size == 0): [baseline_a[s] vint] * n_samples   // 仅 len==3 时有意义
//     [tip_bytes]
//     [perm_flag][perm_present?][perm_vals?]                             // len==3 && PL[0]!=0 时的 min-pos 归一化
//     [payload_len_vint][payload...]
//
// Tips（2-bit）:
//   00: all_zero
//   01: Type2 (b=15*a)                   -> zigzag(delta_a)
//   10: Type1 (len==3 且 |Δa|,|res|<=63)  -> zigzag(delta_a) + zigzag(residual)   // residual = b - 11*a
//   11: tag-based:
//         tag=0: raw u32 bits (vint)
//         tag=1: pattern1 (a,b) fallback (超阈值时回退)
//         tag>=2: dict_id = tag - 2
```

**效果评估**：

| 指标          | LC        | Subset | 评价         |
| ------------- | --------- | ------ | ------------ |
| Tip 10 覆盖率 | **94.1%** | 88.9%  | 优秀         |
| Tip 11 兜底   | 5.9%      | ~11%   | 可接受       |
| 压缩后占比    | 68.37%    | 71.96% | **主要瓶颈** |

**问题**：
1. PL 占压缩后 68-72%，是最大优化目标
2. Tip 10 当前是“字节对齐的两次 varint”，在进入后端压缩前确实是 2B/样本；但这两字节本身熵很低，后端（BSC/Brotli）通常还能继续压缩，因此“预处理层节省 25%”不等价于“最终文件节省 25%”，需要以真实压缩结果为准。

#### 3.2.1 已实现优化：PL codec_id=6（Type2 / Type1 子流拆分）

**动机**：codec_id=4 把不同类型样本的 varint 交错写入同一条 payload；codec_id=6 将 Type2/Tip10 的 (Δa,residual) 拆成独立子流，使后端压缩更容易吃到重复性。

```cpp
// 代码位置:
// - 编码: src/compressor.cpp（FMT PL special codec_id=6）
// - 解码: src/decompression_reader.cpp（FMT PL codec_id=6）
//
// 启用:
//   GSC_PL_CODEC_ID=6
//
// stream layout（codec_id=6, codec marker 后）:
//   [block_size_vint]
//   for each record:
//     if (rec_idx % block_size == 0): [baseline_a[s] vint] * n_samples
//     [tip_bytes]
//     [perm_flag][perm_present?][perm_vals?]
//     [payload_len_vint][payload...]
//
// payload（record）:
//   [len_type2][len_t10_delta][len_t10_residual]
//   [type2_stream...][t10_delta_stream...][t10_residual_stream...][fallback_stream...]
//   - Type2 stream: Tip 01 的 zigzag(Δa)
//   - Tip10 streams: Tip 10 的 zigzag(Δa) 与 zigzag(residual=b-11*a) 分开存
//   - fallback_stream: Tip 11 的 tag+payload（raw / pattern1(a,b) / dict 等）
```

---

## 4. 优化方案

### 4.1 PL 优化（已实现）：子流拆分 (codec_id=6)

**结论（实测）**：相比 codec_id=4，codec_id=6 在 LC 数据集（BSC）上能显著降低 PL 字节数，并且总大小也更小；这是目前优先推荐的无损优化路径。

**启用方式**：

```bash
GSC_PL_CODEC_ID=6 GSC_LOG_LEVEL=debug ./build/gsc compress --in $VCF --out $OUT --compressor bsc
```

**实现要点**：将 Type2 的 `zigzag(Δa)`、Tip10 的 `zigzag(Δa)`、`zigzag(residual)` 拆成独立子流（见 3.2.1），把“高频小 varint”聚到一起，降低交错熵。

**LC（46601×2560, BSC）对比**：

| 配置 | PL Compressed | TOTAL Compressed |
| ---- | ------------- | ---------------- |
| PL4 + AD4 | 152.12 MB | 225.95 MB |
| **PL6 + AD4** | **146.14 MB** | **219.98 MB** |

#### 4.1.1 经验教训：bit-pack (codec_id=5) 在 LC 上会变差

`codec_id=5`（12-bit packed Tip10）在进入后端压缩前确实更紧凑，但在 LC 这类高熵数据上，bit-level 打包会破坏字节级重复结构，使 BSC/Brotli 更难压，实测整体会变差；因此不建议作为主路径（保留为对照/实验）。

### 4.2 AD 优化（已实现）：Δsum 子流拆分 (codec_id=5)；候选：biallelic 杂合 `alt` 小值打包（无损）

#### 4.2.1 已实现：codec_id=5（Δsum 子流拆分）

**启用方式**：

```bash
GSC_AD_CODEC_ID=5 GSC_LOG_LEVEL=debug ./build/gsc compress --in $VCF --out $OUT --compressor bsc
```

**LC（46601×2560, BSC）对比（在 PL6 的前提下）**：

| 配置 | AD Compressed | TOTAL Compressed |
| ---- | ------------- | ---------------- |
| PL6 + AD4 | 56.55 MB | 219.98 MB |
| **PL6 + AD5** | **53.42 MB** | **216.85 MB** |

#### 4.2.2 现状澄清（与当前实现一致）

当前 codec_id=4 对二等位杂合（ref>0 && alt>0）已经有无损快速路径：`tip=11 + tag=1: [Δsum][alt]`，并不需要“完整存储 AD 数组”。因此 AD 的进一步优化空间通常小于 PL，需要先用真实数据统计 `alt` 的取值分布。

#### 4.2.3 不建议：比例量化（会引入有损）

早期草案中的 `ratio_code = ref * 15 / sum` 属于量化近似，无法无损重构 `(ref,alt)`（存在取整误差），**不适用于无损压缩路径**；除非明确引入“有损模式”并给出误差界与验证策略。

#### 4.2.4 可行的无损方向：小值打包 + 逃逸（需先做分布统计）

如果 `alt` 大量集中在小值区间（例如 [0,14] / [0,31]），可以考虑把 `tag=1` 路径中的 `alt` 做小值打包：
- `alt <= 14`：用 4bit 存储 `alt`
- `alt == 15`：作为 escape，后面跟一个 varint(alt)

这种做法对 `alt` 小值密集的数据理论上能把 `alt` 的平均开销从 ~1B/样本压到 ~0.5B/样本，但收益完全取决于分布，需要先用脚本统计再决定是否值得实现。

#### 4.2.5 预期效果

在没有 `alt` 分布统计前，不建议给出固定的“AD 总节省”百分比；建议以 `GSC_LOG_LEVEL=debug` 的 Field Compression Statistics 为最终评估口径。

### 4.3 DP 优化：sum(AD) 异常表

#### 4.3.1 当前状态

- Subset: `DP == sum(AD)` 100%，可完全推断
- LC: `DP == sum(AD)` 86.04%，14% 需要存储

#### 4.3.2 优化方案

```cpp
// 存储策略：
// 1. 假设 DP == sum(AD)
// 2. 仅存储异常项：[sample_id, delta_dp]

struct DPException {
    uint32_t sample_id;
    int16_t delta;  // DP - sum(AD)，通常是小值
};

// 编码：[exception_count][exceptions...]
// LC 14% 异常 → 存储 14% × 样本数 × 约 3B = 远小于全量存储
```

### 4.4 自适应后端选择

#### 4.4.1 策略

```cpp
enum class BackendType { BSC, BROTLI, ZSTD };

BackendType selectBackend(const PreprocessedData& data) {
    // 估算熵
    double entropy = estimateEntropy(data);
    double short_range_redundancy = measureShortRangeRedundancy(data);
    double long_range_correlation = measureLongRangeCorrelation(data);

    // 决策逻辑
    if (entropy > 15.0 && long_range_correlation > 0.8) {
        return BackendType::BSC;  // 高熵 + 长范围相关 → BSC
    } else if (entropy < 13.0 || short_range_redundancy > 0.7) {
        return BackendType::BROTLI;  // 低熵或短范围重复 → Brotli
    } else {
        // 中间情况：小规模测试两者，选更优
        return testAndSelect(data, {BackendType::BSC, BackendType::BROTLI});
    }
}
```

#### 4.4.2 熵估算方法

```cpp
double estimateEntropy(const vector<uint8_t>& data) {
    // 使用样本估算
    array<uint32_t, 256> freq{};
    for (uint8_t b : data) freq[b]++;

    double entropy = 0.0;
    for (uint32_t f : freq) {
        if (f > 0) {
            double p = (double)f / data.size();
            entropy -= p * log2(p);
        }
    }
    return entropy;  // bits per byte
}
```

---

## 5. 实现优先级

| 优先级 | 优化项            | 预期收益          | 实现难度 | 目标字段 |
| ------ | ----------------- | ----------------- | -------- | -------- |
| **P0** | PL 子流拆分（codec_id=6） | 已实测：LC/BSC 下 PL 与 TOTAL 明显下降 | 中 | PL |
| **P1** | AD Δsum 子流拆分（codec_id=5） | 已实测：LC/BSC 下 AD 与 TOTAL 进一步下降 | 低-中 | AD |
| **P2** | AD alt 小值打包（无损） | 数据依赖（需先统计） | 中 | AD |
| **P3** | DP 异常表         | DP -90%，总 -1%   | 低       | DP       |
| **P4** | 自适应后端选择    | 5-15%（数据依赖） | 高       | 全局     |

---

## 6. 实测效果对照

> 说明：以下为实际运行 `GSC_LOG_LEVEL=debug` 的 Field Compression Statistics（最终压缩后字节数），不同后端/数据集会有差异。

### 6.1 LC 数据集（46601 variants × 2560 samples）

后端：BSC（`--compressor bsc`）

| 配置 | PL Compressed | AD Compressed | TOTAL Compressed |
| ---- | ------------- | ------------- | ---------------- |
| PL4 + AD4 | 152.12 MB | 56.55 MB | 225.95 MB |
| PL6 + AD4 | 146.14 MB | 56.55 MB | 219.98 MB |
| **PL6 + AD5** | **146.14 MB** | **53.42 MB** | **216.85 MB** |

### 6.2 Toy Subset（`toy/final_subset.vcf.gz`, 1300 variants × 500 samples）

| 后端 | 配置 | PL Compressed | AD Compressed | TOTAL Compressed |
| ---- | ---- | ------------- | ------------- | ---------------- |
| BSC | PL4 + AD4 | 726 KB | 160 KB | 1008 KB |
| BSC | PL6 + AD4 | 647 KB | 160 KB | 929 KB |
| BSC | PL6 + AD5 | 647 KB | 157 KB | 927 KB |
| Brotli | PL4 + AD4 | 658 KB | 153 KB | 923 KB |
| Brotli | PL6 + AD4 | 597 KB | 153 KB | 862 KB |

---

## 7. 测试验证

### 7.1 无损验证

```bash
# 压缩
./build/gsc compress --in toy/final_subset.vcf.gz --out /tmp/test.gsc --compressor brotli

# 解压
./build/gsc decompress --in /tmp/test.gsc --out /tmp/test.vcf

# MD5 校验
zcat toy/final_subset.vcf.gz | md5sum
md5sum /tmp/test.vcf
# 必须一致
```

### 7.2 压缩率验证

```bash
# 对比（LC 推荐先看 BSC）
# Baseline: PL4 + AD4
GSC_LOG_LEVEL=debug GSC_PL_CODEC_ID=4 GSC_AD_CODEC_ID=4 ./build/gsc compress \
    --in /home/binq/data/1000GPi/lc_bams.first50000.vcf.gz \
    --out /tmp/lc_pl4_ad4.gsc \
    --compressor bsc

# Optimized: PL6 + AD5
GSC_LOG_LEVEL=debug GSC_PL_CODEC_ID=6 GSC_AD_CODEC_ID=5 ./build/gsc compress \
    --in /home/binq/data/1000GPi/lc_bams.first50000.vcf.gz \
    --out /tmp/lc_pl6_ad5.gsc \
    --compressor bsc

# 查看 Field Compression Statistics
```

### 7.3 不同后端对比

```bash
# BSC
./build/gsc compress --in $VCF --out /tmp/test_bsc.gsc --compressor bsc
ls -la /tmp/test_bsc.gsc

# Brotli
./build/gsc compress --in $VCF --out /tmp/test_brotli.gsc --compressor brotli
ls -la /tmp/test_brotli.gsc
```

---

## 8. 总结

### 8.1 BSC vs Brotli 结论

| 数据特征                               | 推荐后端   |
| -------------------------------------- | ---------- |
| 低覆盖度、高杂合、模式多样 (如 LC)     | **BSC**    |
| 高覆盖度、高纯合、模式规整 (如 Subset) | **Brotli** |

根本原因：**数据熵和冗余类型不同**。BSC 的 BWT 擅长长范围相关，Brotli 的 LZ77 擅长短范围重复。

### 8.2 优化路线图

1. **已完成 (P0)**：PL 子流拆分（`GSC_PL_CODEC_ID=6`），对 LC/BSC 实测有效
2. **已完成 (P1)**：AD Δsum 子流拆分（`GSC_AD_CODEC_ID=5`），在 PL6 基础上进一步降低 AD/TOTAL
3. **后续 (P2)**：统计 AD 二等位杂合 `alt` 分布，评估 `alt` 小值打包（无损）
4. **后续 (P4)**：自适应后端选择（数据依赖）

---

*文档版本: 2026-01-11 (codec6+ad5 update)*
