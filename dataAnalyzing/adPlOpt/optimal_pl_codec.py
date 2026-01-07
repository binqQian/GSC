#!/usr/bin/env python3
"""
最优 PL 编码方案分析与对比

方案设计：
1. a 值差分编码：利用位点间 a 值高度相关（66-76% 差分在 [-3,3]）
2. b = ratio * a + residual：利用 b/a 比值规律（74-93% 在 10-12）
3. 组合编码：tip(2bit) + payload

编码格式：
  Tip 00: a=0, b=0 (全零，无存储)
  Tip 01: delta_a 小，ratio 匹配 → [delta_a_zigzag][residual_zigzag] (最紧凑)
  Tip 10: delta_a 小，ratio 不匹配 → [delta_a_zigzag][b_vint]
  Tip 11: delta_a 大或其他 → [a_vint][b_vint] (兜底)
"""

import sys
import gzip
from collections import defaultdict

def zigzag_encode(v):
    """Zigzag 编码有符号整数"""
    return (v << 1) ^ (v >> 31)

def varint_len(v):
    """Varint 编码字节数"""
    if v < 0: v = zigzag_encode(v)
    if v < 128: return 1
    if v < 16384: return 2
    if v < 2097152: return 3
    return 4

def analyze_optimal_codec(vcf_path, max_variants=20000, max_samples=500):
    """分析最优 PL 编码方案"""

    # 配置参数
    RATIOS = [10, 11, 12]  # 支持的比值
    DELTA_A_THRESHOLD = 15  # delta_a 小的阈值
    RESIDUAL_THRESHOLD = 7  # residual 小的阈值

    # 统计
    total_samples = 0
    prev_a = {}  # sample_id -> prev_a

    # 各方案字节数
    bytes_baseline = 0     # 当前编码：a + b
    bytes_ratio = 0        # 方案1：比值编码
    bytes_delta_ratio = 0  # 方案2：差分 + 比值编码
    bytes_optimal = 0      # 方案3：最优组合编码

    # 分类统计
    tip_counts = defaultdict(int)

    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        sample_count = 0
        variant_count = 0
        pl_idx = -1

        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                sample_count = min(len(parts) - 9, max_samples)
                continue

            if variant_count >= max_variants:
                break
            variant_count += 1

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            fmt = parts[8].split(':')
            try:
                pl_idx = fmt.index('PL')
            except ValueError:
                continue

            for s_idx in range(min(len(parts) - 9, sample_count)):
                sample_data = parts[9 + s_idx].split(':')
                if pl_idx >= len(sample_data):
                    continue

                pl_str = sample_data[pl_idx]
                if pl_str == '.' or pl_str == './.':
                    continue

                try:
                    pl_vals = [int(x) if x != '.' else -1 for x in pl_str.split(',')]
                except:
                    continue

                if len(pl_vals) != 3:
                    continue

                # 只分析 Type 1: [0, a, b] 且非 Type 0/Type 2
                if pl_vals[0] != 0 or pl_vals[1] < 0 or pl_vals[2] < 0:
                    continue

                a, b = pl_vals[1], pl_vals[2]

                # 跳过 Type 0 和 Type 2
                if a == 0 and b == 0:
                    continue
                if a > 0 and b == 15 * a:
                    continue

                total_samples += 1

                # 计算差分
                delta_a = a - prev_a.get(s_idx, a)  # 首次使用自身作为基线
                prev_a[s_idx] = a

                # 计算最佳 ratio 和 residual
                best_ratio = None
                best_residual = None
                if a > 0:
                    for r in RATIOS:
                        res = b - r * a
                        if abs(res) <= RESIDUAL_THRESHOLD:
                            best_ratio = r
                            best_residual = res
                            break

                # ========== 方案对比 ==========

                # 基线：a + b
                bytes_baseline += varint_len(a) + varint_len(b)

                # 方案1：比值编码 (a + ratio_code + residual or a + b)
                if best_ratio is not None:
                    bytes_ratio += varint_len(a) + 1  # a + 1字节(ratio+residual)
                else:
                    bytes_ratio += varint_len(a) + varint_len(b)

                # 方案2：差分 + 比值 (delta_a + ratio_code + residual or delta_a + b)
                if abs(delta_a) <= DELTA_A_THRESHOLD:
                    if best_ratio is not None:
                        bytes_delta_ratio += varint_len(zigzag_encode(delta_a)) + 1
                    else:
                        bytes_delta_ratio += varint_len(zigzag_encode(delta_a)) + varint_len(b)
                else:
                    bytes_delta_ratio += varint_len(a) + varint_len(b)

                # 方案3：最优组合编码 (2bit tip + payload)
                small_delta = abs(delta_a) <= DELTA_A_THRESHOLD

                if small_delta and best_ratio is not None:
                    # Tip 01: delta_a + residual (最紧凑)
                    tip_counts['01'] += 1
                    # delta_a: zigzag, 大多数 1 字节
                    # residual: zigzag, -7~7 = 1 字节
                    bytes_optimal += varint_len(zigzag_encode(delta_a)) + varint_len(zigzag_encode(best_residual))
                elif small_delta:
                    # Tip 10: delta_a + b
                    tip_counts['10'] += 1
                    bytes_optimal += varint_len(zigzag_encode(delta_a)) + varint_len(b)
                else:
                    # Tip 11: a + b (兜底)
                    tip_counts['11'] += 1
                    bytes_optimal += varint_len(a) + varint_len(b)

    # 输出结果
    print(f"="*70)
    print(f"文件: {vcf_path}")
    print(f"分析位点数: {variant_count}, 样本数: {sample_count}")
    print(f"Type 1 样本总数: {total_samples}")
    print(f"="*70)

    print(f"\n[Tip 分布 - 最优方案]")
    for tip in ['00', '01', '10', '11']:
        cnt = tip_counts[tip]
        pct = 100 * cnt / total_samples if total_samples > 0 else 0
        print(f"  Tip {tip}: {cnt:,} ({pct:.1f}%)")

    print(f"\n[编码字节数对比]")
    print(f"  基线 (a + b):              {bytes_baseline:,} bytes")
    print(f"  方案1 (比值编码):          {bytes_ratio:,} bytes, 节省 {100*(1-bytes_ratio/bytes_baseline):.1f}%")
    print(f"  方案2 (差分+比值):         {bytes_delta_ratio:,} bytes, 节省 {100*(1-bytes_delta_ratio/bytes_baseline):.1f}%")
    print(f"  方案3 (最优组合):          {bytes_optimal:,} bytes, 节省 {100*(1-bytes_optimal/bytes_baseline):.1f}%")

    # Tip 开销估算
    tip_bytes = (total_samples * 2 + 7) // 8
    total_optimal = bytes_optimal + tip_bytes
    print(f"\n  方案3 + Tip 开销:          {total_optimal:,} bytes, 净节省 {100*(1-total_optimal/bytes_baseline):.1f}%")

    # 每样本平均字节数
    print(f"\n[每样本平均字节数]")
    print(f"  基线:     {bytes_baseline/total_samples:.2f} bytes/sample")
    print(f"  最优方案: {total_optimal/total_samples:.2f} bytes/sample")

    print()
    return {
        'baseline': bytes_baseline,
        'optimal': total_optimal,
        'saving_pct': 100*(1-total_optimal/bytes_baseline)
    }

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf.gz> [max_variants] [max_samples]")
        sys.exit(1)

    vcf = sys.argv[1]
    max_var = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    max_samp = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    analyze_optimal_codec(vcf, max_var, max_samp)
