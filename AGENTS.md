# GSC 项目上下文（以源码为准 / AI Agent 知识库）

> 本文档面向“需要改代码的人”（含 AI Agent）。所有描述以 `/src` 与实际实现为准；若与旧文档/论文不一致，以源码为准。

---

## 1. 项目速览

### 1.1 GSC 是什么

**GSC (Genotype Sparse Compression)** 是一个面向大型群体 VCF/BCF 的基因型（GT）为核心的压缩/解压与查询工具，支持：

- **多样本 VCF/BCF → `.gsc`**（主线，支持 lossy/lossless 两种“打包模式”）
- **单样本 gVCF → gVCF 专用压缩格式**（独立模块，文件头 `GVCF_FILE_MAGIC`，支持范围查询）

> 说明：项目里“lossy/lossless”更多是在描述“是否把 FILTER/INFO/FORMAT 等 *其他字段* 一并打包进 `.gsc`”。GT 与 VCF 固定列字段（CHROM/POS/ID/REF/ALT/QUAL）在两种模式都会参与。

### 1.2 CLI 入口（与 `src/main.cpp` 一致）

```bash
gsc compress           # 多样本 VCF/BCF 压缩
gsc decompress         # 多样本 .gsc 解压/查询（可输出 VCF/BCF/BED/PGEN/BGEN）
gsc gvcf               # 单样本 gVCF 压缩（独立格式）
gsc gvcf-decompress    # gVCF 解压
gsc gvcf-query         # gVCF 范围查询
```

关键参数（摘要）：
- `-M/--mode_lossly`：选择 lossy 模式（默认是 lossless）。**解压时必须与文件模式匹配**（见 `Decompressor::decompressProcess()`）。
- `--compressor bsc|zstd|brotli`：选择后端压缩（GT 与多个子流会用到）。
- `--max-block-rows`：每个 GT row_block 的最大变体数（默认 `10000`）。
- `--max-block-cols`：每个 GT column_block 的最大单倍体数（默认 `10000`，0 或 ≥ 总单倍体表示不切列）。

---

## 2. 目录结构与关键文件（与实际目录一致）

```
GSC/
├── src/
│   ├── main.cpp                         # CLI + 任务路由
│   ├── gsc_params.h                     # 全局参数结构 GSC_Params + 枚举
│   ├── defs.h                           # 常量与宏（no_variants_in_buf 等）
│   ├── logger.cpp                       # spdlog 初始化（配合 include/logger.h）
│   │
│   ├── compressor.h / compressor.cpp    # 多样本压缩主流程
│   ├── decompressor.h / decompressor.cpp# 多样本解压/查询/格式转换
│   ├── compression_reader.*             # 读取 VCF/BCF，构建 GT bitplane + 固定字段
│   ├── decompression_reader.*           # 读取 `.gsc`，解码 fixed fields 与 GT index
│   ├── block_processing.*               # GT 置换/差分/XOR/稀疏化 + copy/zero 标注
│   │
│   ├── compression_strategy.*           # 压缩后端抽象（BSC/Zstd/Brotli）
│   ├── bsc.*                            # libbsc 封装（CBSCWrapper）
│   ├── zstd_compress.*                  # zstd 封装（支持 ptr 解压避免拷贝）
│   │
│   ├── bit_memory.*                     # CBitMemory：位级写入/读取
│   ├── vint_code.*                      # vint_code：数组变长编码（0 作为终止/零编码）
│   ├── utils.*                          # append/read/Barrier/merge_sort 等工具
│   ├── queues.h                         # 多线程队列（GtBlockQueue 等）
│   ├── file_handle.*                    # File_Handle_2：lossless other fields 的多流容器
│   ├── variant.h                        # variant_desc_t / fixed_field_block / fixed_field_chunk 等
│   ├── field_stats.h                    # 可选：字段级压缩统计
│   │
│   ├── fmt_compress/                    # lossless：FORMAT 特殊字段压缩（AD/DP/PL/GQ/PGT/PID…）
│   │   ├── fmt_field_processor.h
│   │   ├── fmt_dictionaries.h
│   │   ├── fmt_utils.h
│   │   └── vint_codec.h                 # EBML 风格 VInt（注意：与 vint_code 不同）
│   │
│   ├── gvcf/                            # gVCF 独立模块（单样本）
│   │   ├── gvcf_compressor.*            # gVCF 压缩/解压/query
│   │   ├── gvcf_block.*                 # gVCF 块结构
│   │   ├── gvcf_encoding.*              # RLE/Delta/Mask/Dict 等编码
│   │   └── gvcf_field_*.*               # 字段级压缩/解压
│   │
│   ├── parallel_vcf_reader.*            # VCF.gz 并行读/解析（可选）
│   └── parallel_vcf_writer.*            # VCF/BCF 并行写（单线程写出 + 多线程准备）
│
├── include/
│   ├── htslib/                          # 解析 VCF/BCF
│   ├── sdsl/                            # rrr_vector / rank_support 等
│   ├── cpp-mmf/                         # mmap 读取（decompress 侧）
│   ├── libbsc.h
│   └── logger.h                         # LogManager（spdlog wrapper）
│
├── lib/                                 # 预编译静态库（hts/sdsl/bsc/zstd）
└── docs/                                # 文档（本文件）
```

> 纠错：`include/` 里并没有 `zstd-1.5.2/` 这样的源码目录；zstd 是以静态库形式放在 `lib/libzstd.a`，并在 `src/zstd_compress.*` 封装调用。

---

## 3. 核心数据结构（关键类型以 `src/variant.h` 为准）

### 3.1 变体描述（固定列字段的逻辑表示）

`variant_desc_t`（`src/variant.h`）承载 CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO/FORMAT 等文本字段，压缩侧会把其中一部分转为二进制块写入；解压侧会重建这些字段用于输出。

### 3.2 Fixed fields 的块/Chunk（用于范围查询与输出）

`fixed_field_block`（`src/variant.h`）：
- `chrom/id/alt/qual/pos/ref`：各字段对应的压缩 payload
- `gt_block`：**GT 索引/稀疏矩阵编码后的 payload**（不是原始 GT）
- `no_variants`：该 row_block 内的变体数

`fixed_field_chunk`（`src/variant.h`）：
- `row_blocks[]`：多个 `fixed_field_block`（按 row_block 切分）
- `row_meta[]`：每个 row_block 的位置范围信息（`first_pos/last_pos`）
- `gt_row_blocks[]`：**tiled 模式**下每个 row_block 汇总的 GT index（包含所有 column_block 的段）

### 3.3 Column tiling 元数据

项目支持按“单倍体维度”切列以控制内存（`max_block_cols`）。在 on-disk 侧会保存 `n_col_blocks` 与每个 column_block 的 `(start_haplotype, n_haplotypes)` 信息（详见 §6）。

### 3.4 Lossless other fields 的 Key 与 Field 表示

`key_desc`（`src/variant.h`）：
- `keys_type`: `flt/info/fmt`
- `type`: `BCF_HT_*`（FLAG/INT/REAL/STR）
- `name`: 字段名（例如 `DP`、`AD`、`PL`…）

`field_desc`（`src/variant.h`）：
- `present`: 本条记录是否存在该字段
- `data/data_size`: 序列化后 payload（注意 move/copy 语义与内存释放）

---

## 4. 多样本 GT 压缩：数据建模与位级编码（关键纠错点）

### 4.1 GT 并不是“2bit/基因型字符串”，而是“2bit/单倍体 allele”

在 `CompressionReader::addVariant()` 中，GT 先被解析为 allele index（每个样本有 `ploidy` 个 allele），然后用 **两条 bitplane** 表示（因此每个变体输出 `2 * total_haplotypes` 位）：

- 第一条向量（MSB plane）：
  - allele ∈ {0,1} → 0
  - 其他（如 allele=2 或 missing）→ 1
- 第二条向量（LSB plane）：
  - allele ∈ {1,2} → 1
  - allele=0 或 missing → 0

因此解码端在 `Decompressor::initialLut()` 中使用两个字节（MSB/LSB）恢复每个 haplotype 的字符：

| MSB | LSB | 解码字符       |
| --- | --- | -------------- |
| 0   | 0   | `0`            |
| 0   | 1   | `1`            |
| 1   | 0   | `.`（missing） |
| 1   | 1   | `2`            |

> 这也是为什么 `.gsc` 的“向量数”是 `2 * 变体数`：每个变体贡献两条向量。

### 4.2 多等位位点处理（源码行为必须在文档中明确）

`CompressionReader::ProcessFixedVariants()` 对 `n_allele > 2` 的记录并不是通用支持：
- `n_allele <= 2`：常规二等位位点
- 特殊情况：当 `n_allele == 3` 且第三个 ALT 为 `<M>` 或 `<N>` 时，会将位点拆分/改写并继续处理
- 其他更复杂的多等位情况，压缩逻辑并不保证覆盖

这会影响：
- 输出变体数与原始 VCF 记录数不一定 1:1
- lossless 模式会记录 `actual_variants[]` 用于解压时匹配（见 `compressor.cpp` 写入 part2 参数）

---

## 5. 多样本压缩主流程（从源码抽象出的真实架构）

### 5.1 总体数据流

```
VCF/BCF 输入
  ↓ (htslib / 可选 ParallelVCFReader 解析 VCF.gz)
CompressionReader
  ├─ fixed fields: CHROM/POS/ID/REF/ALT/QUAL（写入 fixed_field_*）
  ├─ GT: bitplane (MSB/LSB) → 位流缓冲（CBitMemory）
  └─ (lossless) FILTER/INFO/FORMAT → field_desc[] → CBuffer 分流
  ↓
GtBlockQueue(block_id, col_block_id, data, num_rows, v_vcf_data)
  ↓ (多个 GT 处理线程)
BlockProcess::ProcessSquareBlock()
  ├─ 单倍体置换（clustering，popcnt/bit_cost）
  ├─ XOR 差分编码（按 byte 内相邻 bit）
  ├─ 稀疏化：把每条向量转成“1 的位置列表”（delta + 0 终止）
  └─ Copy/Zero 标注（窗口复制 depth = max_replication_depth）
  ↓
fixed_field_chunk 聚合（row_block + gt_row_blocks）
  ↓
compressFixedFieldsChunk() → writeTempChunkRB()（写临时文件）
  ↓
writeCompressFlie() 拼装 `.gsc` 并追加 SDSL rrr_vector（zeros/copies）
```

### 5.2 Row block 与 Column tiling（控制维度与内存）

- **row_block（行块）**：按变体数切分，阈值 `params.max_block_rows`（默认 10000）
  - GT 向量条数 = `2 * row_block_variants`
- **column_block（列块）**：按单倍体数切分，阈值 `params.max_block_cols`（默认 10000）
  - `n_col_blocks = ceil(total_haplotypes / max_block_cols)`
  - 每块的 `vec_len = ceil(col_block_size / 8)`（按 byte 存 8 个 haplotype 位）

设计要点：
- column tiling 启用时，**只有最后一个 `col_block_id`** 会携带 `v_vcf_data_compress`（变体描述）；其他列块的该 vector 为空（见 `CompressionReader::addVariant()` 的 push 逻辑）。下游逻辑必须显式处理这一点。

### 5.3 BlockProcess 的“稀疏矩阵”编码（源码级语义）

BlockProcess 输出的核心 payload 是 `samples_indexes`（`vector<uint8_t>`），它来自：
- 先得到某种“变换域”的 0/1 向量（包含置换与 XOR 差分）
- 对每条 **非零且非 copy** 的向量：
  - 按从小到大枚举所有 `1` 的 bit 位置（0-based）
  - 写入 **delta 编码**（当前 1 的位置 - 上一个 1 的位置），但编码端用 1-based 临时变量确保 delta ≥ 1
  - 行结束追加一个 `0` 作为终止符
- 全部整数序列用 `vint_code::EncodeArray()` 编码为 bytes（注意其“0 → '\\0'”的约定）

解码端（`Decompressor::decoded_vector_row()`）会：
- `vint_code::DecodeArray()` 得到 delta 序列
- 用 `prev=-1` 的方式恢复 0-based 的 1 位置
- 把 bit 置回到 `cur_decomp_data[byte]` 中

> 重要：`0` 在该编码里是“行终止符/零值编码”，编码侧必须保证不会产生合法的 delta=0。

### 5.4 XOR 差分编码是“byte 内相邻 bit 的 XOR”（与 `initialXORLut` 对齐）

解压端 `Decompressor::initialXORLut()` 建了一个 `map_t256[256]`：
- 输入：一个字节（其 bit 表示“相邻 bit 的 XOR 差分”）
- 输出：恢复后的原始字节

这意味着压缩侧在稀疏化前（或在“变换域”构造中）使用了“byte 内差分 XOR”，以降低 1 的密度并提升稀疏编码效果。

### 5.5 Copy/Zero 标注与 copy 映射压缩（核心性能点）

BlockProcess 会对每条向量：
- 全零向量：`zeros_only[i]=true`
- 在窗口 `max_replication_depth` 内遇到完全相同的向量：标记为 copy，并记录其 origin

全局汇总后（`Compressor::compressReplicatedRow()`）：
- 建 `zeros_only_bit_vector[2]`、`copy_bit_vector[2]`（按 vec_id parity 分离）
- 建 `unique` 位图并用 `rank_support_v5` 计算映射
- 将 copy 的“原始 unique id 差值”用 `used_bits_cp` 位打包写入 `bm_comp_copy_orgl_id`（`CBitMemory`）
- 最终把 zeros/copy 两组位图以 `sdsl::rrr_vector<>` 形式序列化写到文件尾（`sdsl_offset`）

解码端据此可以快速判定：
- 某条向量是全零 / copy / unique
- copy 时用 bit-packed 差值映射找到 origin unique 向量
- unique 向量可选缓存（见 §7.3）

---

## 6. 多样本 `.gsc` 文件格式（真实的 on-disk 结构）

### 6.1 顶层文件头（注意：没有全局 Magic）

多样本 `.gsc` 的起始不是 `GSCF` magic，而是：

1. `bool mode_type`：是否包含 lossless 的 other fields（`true`=lossless，`false`=lossy）
2. `uint64_t other_fields_offset`：other fields 区域起始偏移
3. `uint64_t sdsl_offset`：SDSL rrr_vector 区域起始偏移

解压端会严格校验 CLI 的 `-M/--mode_lossly` 与 `mode_type` 是否匹配（`Decompressor::decompressProcess()`）。

### 6.2 主 Archive 区（从 offset=17 开始，到 `other_fields_offset` 结束）

这部分由 `Compressor::writeCompressFlie()` 写入，包含（顺序按代码）：

- `chunks_streams_size` 与 `chunks_streams[]`（每个 chunk 的实际位置 + 在 archive 内的 offset）
- `ploidy`、`max_block_rows`、`max_block_cols`
- `vec_len`、`no_vec`、`no_copy`
- `used_bits_cp`、`bm_comp_cp_size`、`bm_comp_copy_orgl_id bytes`
- `n_samples`
- `chunks_min_pos[]`（每 chunk 的最小 pos）
- `where_chrom[]`（染色体边界映射）
- GT column tiling 元数据：`n_col_blocks` + `(start_hap, block_size)` * `n_col_blocks`
- 置换表：
  - legacy：`map<row_block_id, vint_bytes>`
  - tiled：`map<(row_block_id,col_block_id), vint_bytes>`
- meta：压缩后的 VCF header 文本与样本名列表（每段以 `[uint32 size][bytes]` 存）
- fixed fields chunks：来自临时文件拼接而来（每个 chunk 前还有一层“payload_size”包裹）

> 平台注意：archive 中出现 `size_t`（例如 `chunks_streams` 的 offset 字段），因此当前格式默认假定 64-bit little-endian 平台。若要跨平台/跨 ABI，需引入显式宽度与版本号。

### 6.3 fixed fields chunk 内部格式：`GSCF`/`GSC_FIXED_FIELDS_RB_VERSION_*`

在每个 fixed fields chunk 内部，存在一个“行块目录格式”，以 `GSC_FIXED_FIELDS_RB_MAGIC = 0x46435347 ("GSCF")` 识别（见 `variant.h` 与 `Compressor::writeTempChunkRB()`）：

- header：`magic/version/total_variants/row_block_count/flags`
- directory：每个 row_block 的 `variant_count/first_pos/last_pos` + 各字段的 `(off,size)`；v2 还包含 per-row_block 的 `gt_off/gt_size`
- data region：按 directory 定义的顺序紧跟各字段 payload

解压端 `DecompressionReader::readFixedFields()` 会自动识别：
- 新格式：读入目录，支持 `DecoderByRange()`（无需解码整个 chunk 即可按 pos 裁剪）
- 旧格式：以 “no_variants 开头的直排块” 方式读取（兼容历史文件）

### 6.4 other fields 区（lossless 才存在）

当 `mode_type==true` 时，`other_fields_offset .. sdsl_offset` 区间保存的是 `File_Handle_2` 管理的多流容器（`src/file_handle.*`）：
- 多个 stream（按 key/type 划分 size/data 子流）
- 尾部包含 footer（stream 目录与偏移表），`deserialize()` 通过从文件尾回读 footer_size 来加载目录

此外，lossless 模式还会在 `part2_params` 流中写入：
- `actual_variants[]`
- `keys[]`（含字段名，用于 FORMAT 特殊 codec）
- `params.backend`（以及 FORMAT dictionaries blob：AD/PL/PID）

### 6.5 SDSL 区（文件尾，从 `sdsl_offset` 开始）

尾部序列化了四个 `sdsl::rrr_vector<>`：
- `zeros` 两个（parity=0/1）
- `copies` 两个（parity=0/1）

解压端会先从 `sdsl_offset` 读取这些结构，再 mmap/读取主 archive 区。

---

## 7. 多样本解压/查询与格式转换（核心实现点）

### 7.1 解压模式与一致性校验

`Decompressor::decompressProcess()` 的关键行为：
- 读取 `.gsc` 头部后，若文件是 lossless（`file_mode_type=true`），**禁止**用户用 `-M/--mode_lossly` 解压；反之亦然
- lossless 模式若指定 `--samples`，会进入“样本查询”路径：**只输出 fixed fields + GT**，不解码 other fields（日志会提示）

### 7.2 Range query 的加速路径：`fixed_fields_rb_dir + DecoderByRange`

当：
- 用户给了 `--range`
- 文件包含新 fixed-fields 行块目录（`has_fixed_fields_rb_dir=true`）
- 且不是 legacy path

解压会走 `DecompressionReader::DecoderByRange()`：只解码覆盖范围内的 row_block/记录，降低 I/O 与 CPU。

### 7.3 Unique 向量缓存（`MB_memory`）

`Decompressor::decoded_vector_row()` 支持缓存 unique 向量：
- `done_unique: unordered_map<uint64_t, uint8_t*>` 存储已解码的 unique row
- `stored_unique` 维护一个简单的淘汰顺序
- `max_stored_unique` 来自 `params.max_MB_memory`（以 MB 估算）或默认不限制

该缓存能显著加速 copy-heavy 或重复访问场景，但也会带来：
- 堆内存管理复杂度（必须确保清理时释放 `uint8_t*`）
- 对不同输出路径（VCF/bed/pgen/bgen）的释放点要一致（源码中多处 `done_unique.clear()` 清理）

### 7.4 输出格式与快速转换

多样本解压支持输出：
- VCF/BCF（默认）
- PLINK1：BED/BIM/FAM（`BedFormatDecompress*`）
- PLINK2：PGEN/PVAR/PSAM（实验性；`initGscToPgenLUT()`/`finalizePgen()` 会 patch header）
- BGEN v1.2（实验性；`initGscToBgenLUT()`/`finalizeBgen()`）

这些路径通常会利用预计算 LUT（例如 `gt_lookup_table[256][256][8]`）减少 per-sample 字符串拼接成本。

---

## 8. Lossless：FILTER/INFO/FORMAT 的真实压缩策略（文档缺失补全）

### 8.1 Key 收集与字段顺序图（field_order_graph）

在 lossless 模式下，`CompressionReader::InitVarinats()` 会：
- 从 header 中收集 FILTER/INFO/FORMAT 的 key 列表（`keys[]`）
- 为每个 key 分配 `CBuffer`（`max_buffer_size` 控制 flush）

每条 record 在解析时会记录该 record 中字段的出现顺序，并累积到 `field_order_graph`，最终 `topo_sort()` 得到一个全局 order；压缩结束前调用 `CompressionReader::UpdateKeys(keys)` 用这个 order 修正 `actual_field_id`，并写入 `part2_params`，确保解压侧的解释一致。

### 8.2 FORMAT 特殊字段 codec（重要优化点）

源码对部分 FORMAT 字段做了“语义级”压缩，而不是纯字节流压缩（位置主要在 `CompressionReader::GetVariantFromRec()`，并复用 `fmt_compress/*` 工具）：

- `AD`：
  - 识别常见模式（例如总和为 0、只有第一个分量非零等）
  - 使用 2-bit tip + 可选字典 ID（`FmtDictionaries::getADItemId`）
- `DP`：
  - 可由 `sum(AD)` 预测，存储 exceptions/原始值
- `MIN_DP`：
  - 可由 `DP` 预测，存储 exceptions 或 raw fallback
- `PL`：
  - 识别多个 pattern（含 `b==15*a` 特例），否则落入字典（`getPLItemId`）
- `GQ`：
  - 由 `PL` 的二小值/特定 pattern 预测，仅存 exceptions
- `PGT`：
  - 稀疏存储（只记录非 `.` 的位置与短字符串）
- `PID`：
  - 稀疏 + 字典（`getPIDItemId`）

同时，压缩结束会把 AD/PL/PID 的 dictionaries blob（+ 版本号）写入 `part2_params`，以便解压复原。

---

## 9. gVCF 模块（单样本独立格式）

gVCF 模块与多样本 `.gsc` **不是同一个文件格式**。

### 9.1 文件头与版本

`src/gvcf/gvcf_compressor.h`：
- `GVCF_FILE_MAGIC = 0x47564346 ("GVCF")`
- `GVCF_FILE_VERSION = 4`（v4 具备 block index 用于 range query）

### 9.2 编码组件（`gvcf_encoding.*`）

该模块实现了若干可序列化的编码器：
- RLE（字符串/整数/字节）
- Delta（POS/END）
- Mask（dominant value + bitmask + patches）
- Dictionary（字符串/整数）

并在 `GVCFCompressor/GVCFDecompressor/GVCFQueryer` 中组织为 block 级压缩与索引。

---

## 10. 开发者注意事项（最容易踩坑的点）

1. **解压模式必须匹配文件模式**：lossless 文件不能用 `-M` 解压；lossy 文件必须带 `-M`（代码会直接报错退出）。
2. **GT 编码的 4 态含义不可改动**：`00/01/10/11 → 0/1/./2` 的假设贯穿 LUT、稀疏化与多种输出格式。
3. **vint_code 的 `0` 是约定的终止/零编码**：任何会产生 delta=0 的改动都会破坏解码（尤其是稀疏矩阵行编码）。
4. **XOR 差分是按 byte 内相邻 bit**：压缩侧任何关于“差分域”的改动都必须同步更新 `initialXORLut()` 的逆变换逻辑。
5. **column tiling 下只有最后一个列块携带 variant_desc**：压缩线程里不能假设每个 `GtBlockQueue` pop 都有 `v_vcf_data_io`。
6. **GT row_block payload 尾部 marker**：固定字段块里 GT index 压缩后会 `push_back(marker)`；解压端会 `pop_back()` 再按 marker 选 codec。修改 codec/marker 规则会导致不可解码。
7. **fixed-fields chunk 有版本**：`GSC_FIXED_FIELDS_RB_VERSION_V1/V2` 已在用；变更目录布局必须 bump 版本并保持兼容读取。
8. **`.gsc` 主 archive 含 `size_t`**：若未来要跨平台/跨语言读取，必须把所有 `size_t` 替换为固定宽度并引入全局 magic/version。
9. **临时文件与路径**：
   - lossless 解压会在 `tmp/` 下创建 `gsc_part2_*`，压缩也会创建 `*.temp/.com_tmp_gsc` 等文件；改动时注意异常路径清理。
10. **并发与队列结束信号**：任何新增线程/队列都必须保证 `Complete()` 一定被调用，否则会死锁（尤其是异常/early return）。

---

## 11. 构建与调试（最小可用指引）

```bash
cmake -S . -B build
cmake --build build -j

# toy 示例
./build/gsc compress -i toy/toy.vcf -o toy/test.gsc
./build/gsc decompress -i toy/test.gsc -o toy/restored.vcf
```

日志级别：
```bash
export GSC_LOG_LEVEL=debug   # trace/debug/info/warn/error/critical/off
```

---

## 更新日志

| 日期       | 更新内容                                                                                                                                      |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| 2026-01-29 | 对旧版 `PROJECT_CONTEXT.md` 做“以源码为准”的重构：修正目录结构、GT 位级语义、`.gsc` 文件格式与 lossless FORMAT 特殊 codec，并补充开发者易错点 |

