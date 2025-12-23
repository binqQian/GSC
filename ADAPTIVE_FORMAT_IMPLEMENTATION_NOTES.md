# 自适应 FORMAT 压缩：实现现状/约束/下一步（给下次AI的速读文档）

> 目标：让后续 AI/开发者无需重新通读代码，就能理解“自适应 FORMAT”当前做到了什么、为什么这样做、如何验证、还缺什么。

---

## 1. 背景与目标

### 1.1 背景（legacy 路径）

GSC 的 lossless 压缩大体分层：

- GT：专用稀疏/复制检测 + 位向量编码（GT 主体仍是 GSC 的核心）
- fixed fields：块压缩（BSC/ZSTD/Brotli）
- other fields（INFO + FORMAT 非 GT）：按类型（INT/REAL/STR/FLAG）分 stream，写入后由 `Compressor::compress_other_fileds()` 走后端压缩（BSC/ZSTD/Brotli）

### 1.2 自适应 FORMAT 的目标（“primary” 的最终形态）

1) 正确性：lossless round-trip，`compress -> decompress` 输出 VCF 与 legacy(off) 一致（或可定义的 canonicalization 一致）。  
2) 可控启用：默认不破坏原有行为；支持 shadow/primary 两种运行模式。  
3) 压缩率：primary 模式下整体 `.gsc` 不能显著劣化；尤其不能出现“把 typed binary + 压缩”变成“raw bytes/字符串”导致变大。  
4) 兼容性：新旧 `.gsc` 文件都能解压（旧文件无 adaptive stream，新文件可选择使用 adaptive 或 legacy）。

---

## 2. 预期行为（CLI 与语义）

### 2.1 压缩端

`./gsc compress --adaptive-format off|shadow|primary ...`
`./gsc compress --adaptive-format-compressor auto|follow|raw|bsc|zstd|brotli ...`

- `off`：完全关闭 adaptive（不写 `adaptive_format_data`），与历史行为等价。
- `shadow`：同时写 legacy non-GT FORMAT（other_fields 路径）+ `adaptive_format_data`（用于回归验证；体积可能更大）。
- `primary`：写 `adaptive_format_data`，并且**仅对“adaptive 里实际输出的 tag”抑制 legacy non-GT FORMAT**；其余 tag 保留 legacy（避免 RawString 放大）。

### 2.2 解压端

`./gsc decompress --use-adaptive-format ...`

- 不加 `--use-adaptive-format`：保持旧行为（只用 legacy other_fields 恢复 FORMAT）。
- 加 `--use-adaptive-format`：若输入文件存在 `adaptive_format_data` stream，则优先使用 adaptive 恢复 non-GT FORMAT；对于 adaptive 未覆盖的 tag，仍从 legacy 恢复（保证 primary 的“选择性覆盖”可用）。

备注：当前 adaptive 解压主要面向 VCF 输出；BCF 输出仍优先 legacy typed 更新。

---

## 3. 最近提交做了什么、主要问题是什么（按时间顺序）

### 3.1 `8c0f011`（phase 1&2）

新增/实现：
- `FormatFieldCodec` 接口 + 多个 codec：`BitTipArrayCodec`（AD/PL），`PredictedScalarCodec`（DP/GQ 预测型），`SparseDictCodec`（PGT/PID），`RawString` 兜底。
- `FormatFieldDetector`：采样统计特征并选 codec（`selectCodecType()`）。

潜在问题（当时/默认实现）：
- `selectCodecType()` 会返回未真正实现的 codec 类型（最终回退 RawString），容易误判“应该有收益”。
- 字典/序列化细节（ZigZag、vector_end）若处理不严，会影响 lossless。

### 3.2 `8e68e76`（phase 3&4）

新增/集成：
- `CompressionReader` 侧把 adaptive 数据写入新 stream：`adaptive_format_data`（按 variant 写 row）。
- `DecompressionReader` 侧检测 `adaptive_format_data` stream 的存在。
- `Decompressor` 增加框架字段（use flag、manager）。

当时存在的关键问题（导致压缩率/正确性不佳）：
- `adaptive_format_data` **裸写入、没有任何压缩**；而 legacy other_fields 会走压缩后端 → primary 会明显变大。
- per-record FORMAT tag 集合/顺序没有严格遵循输入 record；可能“凭空带上” header 里有但 record 未出现的 FORMAT tag。
- 解码契约/预测器绑定存在漏洞：codec header/payload 划分、预测器生命周期、sample_pos 对齐等容易出错。
- `BitTipArrayCodec` 字典的 ZigZag/element_bytes 选择存在截断风险（例如 135 被错误还原成 7）。

### 3.3 `35bdd91`（FORMAT字段自适应压缩0）

主要内容：
- 让 `toy/toy.vcf` 上 `primary + --use-adaptive-format` 的解压输出与 legacy/shadow 对齐（端到端可用）。
- 引入 `--adaptive-format off|shadow|primary` 与 `--use-adaptive-format`，并补齐 adaptive stream 的写入/读取/对齐框架。
- 生成了分阶段 patch 文件（`tmp/patches/01..03-*.patch`）用于 review/回归。

遗留问题：
- `adaptive_format_data` 仍是裸写，会让 `primary` 体积劣化（尤其是大样本/多 tag 的真实 VCF）。
- 对 `RawString` 退化字段（例如 GQ/RGQ/MIN_DP/SB 等）若也走 adaptive，会把 typed+压缩变成文本+元数据，压缩率很差。

### 3.4 `7690f56`（FORMAT字段自适应压缩1）

主要内容（针对“primary 压缩率差”与正确性补洞）：
- `adaptive_format_data` part 增加 framing，并按 part 做 zstd 压缩写入；解压端透明解压，兼容旧 raw part。
- primary 模式下支持“选择性覆盖”：adaptive row 可省略 `RawString` 字段，并且 legacy 抑制只针对 adaptive row 实际包含的 tag。
- 解压端按 `keys[]` 顺序逐 tag 混合 legacy/adaptive 更新，避免 FORMAT tag 顺序变化或值错位。
- 修复 `BitTipArrayCodec` 字典 ZigZag 截断（`135 -> 7`）的 element_bytes 选择错误，恢复多等位 PL 的正确性。

---

## 4. 最近完成的改造（本次工作区 patch 的核心内容）

### 4.1 primary 压缩率：把 adaptive stream 变成“可压缩”的 part

文件：
- `src/compression_reader.cpp`
- `src/decompression_reader.cpp`

做法：
- `adaptive_format_data` 的每个 part 写入 **带 framing 的 payload**：
  - magic：`AFD1`
  - method：`0=raw`，`1=zstd`，`2=bsc`，`3=brotli`
  - raw_size：`uint32 little-endian`
  - payload：raw 或 zstd 压缩后的 bytes
- 解压端识别 magic；若为 zstd 则解压后再喂给 row parser；否则兼容旧 raw part。

效果（toy/final_subset.vcf 示例）：
- `shadow`：双写（预期最大）
- `primary`：体积显著下降（见“验证命令/指标”）

### 4.2 primary “选择性覆盖”：避免 RawString 放大

核心点：
- `FormatFieldManager::finalizeRow(out, omit_rawstring_fields)`：primary 模式下可选择**不把 RawString-coded 的 tag 写入 adaptive row**。
- 压缩端 primary 模式下的 legacy 抑制策略调整为：
  - 仅对 **adaptive row 实际包含的 tag** 才 suppress legacy non-GT FORMAT
  - adaptive 未覆盖的 tag 仍保留 legacy typed 压缩（例如 GQ、RGQ、MIN_DP、SB 等）

相关文件：
- `src/format_field_detector.h`
- `src/format_field_detector.cpp`
- `src/compression_reader.cpp`

### 4.3 解压端“混合来源”保持正确与 FORMAT 顺序

之前的风险：
- primary 模式下只抑制部分 legacy tag，会出现“部分 tag 来自 adaptive、部分来自 legacy”；如果用 `bcf_update_format_*` 的顺序不一致，会导致输出 FORMAT tag 顺序变化或值错位。

当前做法：
- 在 `appendVCFToRec()` 内按 `keys[]` 顺序遍历 FORMAT：
  - 若 tag 在 adaptive row 中：从 adaptive 解码并 typed 更新该 tag
  - 否则：按 legacy `_fields` typed 更新该 tag

文件：
- `src/decompressor.cpp`

### 4.4 `BitTipArrayCodec` 字典 ZigZag 截断修复（关键正确性）

问题：
- 字典序列化用 ZigZag（`zz=(v<<1)^(v>>63)`），但 element_bytes 的选择仅基于 `abs(v)`，会低估 ZigZag 后的最大值（例如 `135 -> zz=270` 需要 2 bytes）。

修复：
- element_bytes 选择帮助函数改为按 `max_zz`（近似 `2*abs(v)+1`）估算。

文件：
- `src/codecs/bit_tip_array_codec.cpp`

---

## 5. 当前限制/已知欠缺（下一步待做）

### 5.1 性能

- `--use-adaptive-format` 已实现 “decodeRow once” 优化：每个 variant 只反序列化一次 codecs，避免每个 sample 重复反序列化。
  - 参考：`FormatFieldManager::prepareRowDecoder()` + `Decompressor::appendVCFToRec()` 的按行缓存。
  - 在 `toy/final_subset.vcf` 上，adaptive 解压耗时从 ~30s 级降到 ~1.5s（机器相关），接近 legacy 路径。

### 5.2 codec 覆盖面

- 目前收益主要来自：AD/PL(BitTipArray) + DP<-AD(PredictedScalar)。
- GQ<-PL predictor 仍未完成（若实现可进一步减少 GQ 的存储量）。
- 多等位 PL 的 AB-pattern 目前为了安全只在二等位（len==3）启用；多等位的模式重建需要完整实现 VCF spec 的 PL 索引规则。

### 5.3 primary 策略可继续优化

- 当前 primary 采用“只覆盖非 RawString 的 tag”是一个保守策略，确保压缩率不劣化。
- 下一步可增加策略开关：例如允许对某些 int 标量（GQ/RGQ/MIN_DP）实现更好的 codec 后再纳入 adaptive 覆盖。

---

## 6. 验证命令（建议保留给 CI/回归脚本）

### 6.1 基础一致性（toy/toy.vcf）

```bash
make -j
./gsc compress --adaptive-format off    --in toy/toy.vcf --out tmp/toy_off.gsc
./gsc compress --adaptive-format shadow --in toy/toy.vcf --out tmp/toy_shadow.gsc
./gsc compress --adaptive-format primary --in toy/toy.vcf --out tmp/toy_primary.gsc

./gsc decompress --in tmp/toy_off.gsc --out tmp/toy_off.vcf
./gsc decompress --in tmp/toy_shadow.gsc --out tmp/toy_shadow.vcf
./gsc decompress --use-adaptive-format --in tmp/toy_primary.gsc --out tmp/toy_primary.vcf

diff -q tmp/toy_off.vcf tmp/toy_shadow.vcf
diff -q tmp/toy_off.vcf tmp/toy_primary.vcf
```

### 6.2 你的目标用例（toy/final_subset.vcf）

```bash
./gsc compress --adaptive-format off --in toy/final_subset.vcf --out toy/test_off.gsc
./gsc compress --adaptive-format primary --in toy/final_subset.vcf --out toy/test_fs_1.gsc
ls -lh toy/test_off.gsc toy/test_fs_1.gsc

./gsc decompress --in toy/test_off.gsc --out toy/test_off.vcf
./gsc decompress --use-adaptive-format --in toy/test_fs_1.gsc --out toy/test_fs_1.vcf
diff -q toy/test_off.vcf toy/test_fs_1.vcf
```

（示例指标，环境相关）：`off ~1.8M`，`primary ~2.1M`，`shadow ~4.6M`。

---

## 7. 关键实现位置（下次AI定位入口）

- adaptive 写入入口：`src/compression_reader.cpp` -> `CompressionReader::SetVariantOtherFields()`
- adaptive stream 读取：`src/decompression_reader.cpp` -> `GetNextAdaptiveFormatRow()` + part framing 解压
- adaptive 应用到输出：`src/decompressor.cpp` -> `appendVCFToRec()`（按 key 顺序混合 legacy/adaptive）
- codec 框架：`src/format_field_detector.*`（含 `FormatFieldManager` 实现）
- AD/PL/DP 编码核心：`src/codecs/bit_tip_array_codec.cpp`、`src/codecs/predicted_scalar_codec.*`
