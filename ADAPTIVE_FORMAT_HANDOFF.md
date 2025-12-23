# 自适应 FORMAT（adaptive_format_data）Handoff：目标/语义/现状/最近变更/待做

> 目的：让下次 AI/开发者不需要重新通读代码即可理解：自适应 FORMAT 做到哪一步、为什么 primary 可能变大/变慢、以及下一步该做什么。

## 1) 目标与边界

- 目标：lossless round-trip（`compress -> decompress` 输出 VCF 与 legacy(off) 一致），并在 `primary` 模式下尽量不劣化总体 `.gsc` 体积/解压性能。
- 默认不破坏旧行为：压缩端默认 `off`；解压端默认不用 adaptive（需显式 `--use-adaptive-format`）。
- 自适应 FORMAT 目前主要服务 VCF 输出；BCF 输出仍优先 legacy typed 路径（避免更新策略/顺序差异）。

## 2) CLI 新增参数与语义（相对 legacy 的变化）

### 压缩端

`./gsc compress --adaptive-format off|shadow|primary ...`

- `off`：不写 `adaptive_format_data`（完全等价历史行为）。
- `shadow`：同时写 legacy non-GT FORMAT（other_fields 路径）+ `adaptive_format_data`（用于回归验证；体积通常最大）。
- `primary`：写 `adaptive_format_data`，并**只对 adaptive row 实际包含的 tag 抑制 legacy non-GT FORMAT**；adaptive 未覆盖的 tag 仍走 legacy。
  - `primary` 额外策略：若某 tag 的 codec 退化成 `RawString`，该 tag 会被 **omit** 出 adaptive row（避免“typed+后端压缩”退化成“raw string + 元数据”导致体积暴涨）。
- `--adaptive-format-compressor auto|follow|raw|bsc|zstd|brotli`：选择 `adaptive_format_data` 的 part 压缩后端。
  - `auto`（默认）：在 `zstd/bsc/brotli` 中择优（挑最小）。
  - `follow`：跟随全局 `--compressor`。
  - `raw`：不压缩（仅 framing）。

### 解压端

`./gsc decompress --use-adaptive-format ...`

- 不加该参数：保持旧行为（只从 legacy other_fields 恢复 FORMAT）。
- 加该参数：若输入存在 `adaptive_format_data` stream，则优先用 adaptive 恢复 non-GT FORMAT；adaptive 未覆盖的 tag 仍从 legacy 恢复（保证 primary 的“选择性覆盖”可用）。

## 3) 为什么 `--adaptive-format primary` 可能“压缩率差”

常见原因（按影响从大到小）：

1. `shadow`/误配置导致“双写”：legacy + adaptive 都保留会直接变大；`primary` 才会 suppress legacy。
2. adaptive 的 payload 原先只做 zstd（而 legacy other_fields 默认走 BSC/Brotli 且可跨 chunk 聚合压缩），在某些数据分布下 zstd 对 adaptive bytes 的效果不如 legacy → primary 可能仍略大。
3. codec 覆盖面/收益评估还保守：目前主要收益来自 AD/PL(BitTipArray) + DP<-AD(PredictedScalar)；若真实数据里这些字段在 legacy 已非常可压，adaptive 的 row 级编码未必稳赢。

已做的应对（见 5.3）：adaptive part 现在会在 zstd 与 bsc 间择优压缩，从而显著缩小 primary 的体积劣化。

## 4) 数据格式与关键代码位置

### 4.1 adaptive_format_data 的 part framing

- stream：`adaptive_format_data`（写入/读取在 `File_Handle_2` 内）
- part framing（新）：`AFD1` + `method` + `raw_size(u32le)` + `payload`
  - `method=0` raw；`method=1` zstd；`method=2` bsc；`method=3` brotli
  - payload 解压后为 row stream bytes（见下）
- 兼容性：老文件可能是“无 framing 的 raw part”；读端会自动兼容。

关键文件：
- 写入：`src/compression_reader.cpp`（`encodeAdaptiveFormatPart()` + `adaptive_format_stream_id_` flush）
- 读取：`src/decompression_reader.cpp`（`refillAdaptiveFormatBuffer()`）

### 4.2 row stream 格式

每个 row：
- `[row_size: u32le][row_payload: row_size bytes]`

row_payload：
- `field_count(vint)`
- 循环 `field_count` 次：
  - `name_len(vint) + name(bytes)`
  - `codec_size(vint) + codec_bytes`

codec_bytes 内含 `CodecParams` + codec payload（由各 codec 自己序列化）。

关键文件：
- 生成 row：`src/format_field_detector.cpp`（`gsc::FormatFieldManager::finalizeRow()`）
- 解析 row：`src/format_field_detector.cpp`（`gsc::FormatFieldManager::prepareRowDecoder()`）

### 4.3 解压端应用策略（混合来源 + 保序）

- 入口：`src/decompressor.cpp` -> `Decompressor::appendVCFToRec()`
- 按 `keys[]` 顺序遍历 FORMAT，保证输出 tag 顺序稳定：
  - 若 tag 在 adaptive row 中：从 adaptive 解码并 `bcf_update_format_*`
  - 否则：走 legacy `_fields` 的 typed 更新
- sample subset：legacy non-GT FORMAT 不能 subset；若 `-s/--samples` 开启 subset，只能依赖 adaptive 覆盖的 tag 才能正确 subset。

## 5) 最近两次 commit 评审 + 当前工作区的“下一阶段变更”

### 5.1 `35bdd91`（FORMAT字段自适应压缩0）

合理性：
- 完成了 adaptive 框架打通（codec/detector/manager + stream 写入/读取 + `--adaptive-format`/`--use-adaptive-format`）。
- 端到端（toy）可用，为后续优化留了接口。

潜在问题：
- adaptive stream 裸写或压缩不足，`shadow/primary` 容易显著变大（尤其大量 RawString 字段）。
- 解压端容易出现 per-sample 反复反序列化/字符串拆分开销（可到 30s 级）。

### 5.2 `7690f56`（FORMAT字段自适应压缩1）

合理性：
- 引入 `AFD1` framing 并 zstd 压缩 part；primary 做“选择性 suppress legacy”（仅对 adaptive 实际包含的 tag）。
- primary omit RawString-coded 字段，避免最常见的体积暴涨路径。
- 修复 `BitTipArrayCodec` 字典 ZigZag 截断（正确性关键洞）。

潜在问题/注意点：
- framing 内 `raw_size` 主要用于校验/健壮性；若遇到损坏文件，当前读端仍偏“容错继续”（会 fallback raw），未来可考虑更严格的错误处理。
- codec 选择策略仍偏启发式，可能对某些非典型字段误选（需收益评估/阈值进一步收紧）。

### 5.3 当前工作区（下一阶段）已做的增量改进

1) **保持 legacy 默认行为**：`adaptive_format_mode` 默认改为 `off`（避免默认 shadow 导致文件变大/行为变化）。  
2) **primary 体积优化**：`AFD1` part 现在会在 `zstd` 与 `bsc` 之间择优压缩（`method=1/2`），显著改善 primary 体积劣化。  
3) **解压性能优化**：
   - `BitTipArrayCodec` 反序列化后构建 per-sample payload index cache，避免 `decode()` 的 O(sample_pos) 扫描。
   - 解压端对 `PredictedScalar`/`BitTipArray` 增加 typed fast path（直接 `int32`/固定数组写入），减少字符串解析与 `strtol/strtod` 的热点。

关键文件：
- `src/gsc_params.h`
- `src/compression_reader.cpp`
- `src/decompression_reader.cpp`
- `src/codecs/bit_tip_array_codec.{h,cpp}`
- `src/codecs/predicted_scalar_codec.{h,cpp}`
- `src/decompressor.cpp`

## 6) 已完成修复（汇总）

- `adaptive_format_data` part framing：`AFD1` + method + raw_size + payload（兼容旧 raw part）。
- primary “选择性覆盖”：只 suppress legacy 中 adaptive 实际发出的 tag；RawString-coded tag 默认不进入 adaptive row。
- `BitTipArrayCodec` ZigZag element_bytes 截断修复（lossless）。
- 解压端 row-level codec 反序列化缓存（prepareRowDecoder once）+ typed fast path（减少 per-sample 重复工作）。
- `--adaptive-format` / `--use-adaptive-format` CLI 语义落地（可控启用）。

## 7) 待做任务（建议优先级）

P0（直接影响体积/性能/正确性）：
- codec 选择“收益评估/阈值收紧”：避免对非典型数组字段误用 BitTipArray，必要时在 primary 中回退 legacy。
- 多等位 PL 的 AB-pattern 完整重建（按 VCF spec 的 PL 索引规则），并在正确前提下扩展 GQ<-PL predictor。
- framing 的健壮性：对 raw_size/解压失败策略做更严格的校验与报错（避免静默输出损坏）。

P1（工程化/可维护性）：
- 增加 `gsc inspect` 或 debug 统计：输出各 stream/part 的 size 与压缩比，便于解释“为什么变大/变慢”。
- 把 toy 回归命令固化到 `test_scripts/`（至少 size + diff + time 三件套）。

P2（更大结构调整）：
- 将“按 row 存 FORMAT”升级为“按 tag/按 chunk 存 FORMAT”（减少重复元数据，增强跨 row 压缩），才有机会系统性超过 legacy。

## 8) 验证命令（建议保留）

```bash
make -j

./gsc compress --adaptive-format off     -i toy/final_subset.vcf -o toy/test_off.gsc
./gsc compress --adaptive-format primary -i toy/final_subset.vcf -o toy/test_primary.gsc
ls -lh toy/test_off.gsc toy/test_primary.gsc

./gsc decompress --in toy/test_off.gsc --out tmp/off.vcf
./gsc decompress --use-adaptive-format --in toy/test_primary.gsc --out tmp/primary.vcf
diff -q tmp/off.vcf tmp/primary.vcf
```
