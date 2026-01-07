#!/usr/bin/env python3
"""
最优 PL 编码方案 - 无 Tip 版本

核心思想：
1. a 值差分编码：delta_a = a - prev_a
2. b 值用 (a, ratio) 推导：residual = b - ratio * a
3. 存储：[delta_a_zigzag][residual_zigzag]

无需 Tip，因为解码时可以反推：
  a = prev_a + delta_a
  b = ratio * a + residual

问题：ratio 不固定怎么办？
方案：使用 **自适应 ratio** 或 **联合残差编码**
"""

import sys
import gzip
from collections import defaultdict

def zigzag_encode(v):
    if v >= 0:
        return v * 2
    else:
        return (-v) * 2 - 1

def varint_len(v):
    if v < 0: v = zigzag_encode(v)
    if v < 128: return 1
    if v < 16384: return 2
    if v < 2097152: return 3
    return 4

def analyze_best_codec(vcf_path, max_variants=20000, max_samples=500):
    """分析最优 PL 编码方案"""

    # 统计
    total_samples = 0
    prev_a = {}

    # 各方案字节数
    bytes_baseline = 0          # 基线：a + b
    bytes_fixed_ratio = {}      # 固定 ratio 方案
    bytes_adaptive_ratio = 0    # 自适应 ratio 方案
    bytes_joint_delta = 0       # 联合差分方案

    # 初始化 ratio 统计
    for r in [10, 11, 12]:
        bytes_fixed_ratio[r] = 0

    # 联合差分统计
    prev_b = {}
    delta_b_dist = defaultdict(int)

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

                if pl_vals[0] != 0 or pl_vals[1] < 0 or pl_vals[2] < 0:
                    continue

                a, b = pl_vals[1], pl_vals[2]

                # 跳过 Type 0 和 Type 2
                if a == 0 and b == 0:
                    continue
                if a > 0 and b == 15 * a:
                    continue

                total_samples += 1

                # 基线
                bytes_baseline += varint_len(a) + varint_len(b)

                # delta_a
                delta_a = a - prev_a.get(s_idx, a)
                prev_a[s_idx] = a

                # delta_b (联合差分)
                delta_b = b - prev_b.get(s_idx, b)
                prev_b[s_idx] = b

                if abs(delta_b) <= 31:
                    delta_b_dist['<=31'] += 1
                elif abs(delta_b) <= 127:
                    delta_b_dist['32-127'] += 1
                else:
                    delta_b_dist['>127'] += 1

                # 方案1：固定 ratio（无 ratio 标记）
                for r in [10, 11, 12]:
                    residual = b - r * a if a > 0 else b
                    bytes_fixed_ratio[r] += varint_len(zigzag_encode(delta_a)) + varint_len(zigzag_encode(residual))

                # 方案2：自适应 ratio（需要 2bit ratio 标记）
                if a > 0:
                    best_r = 11
                    best_res = abs(b - 11 * a)
                    for r in [10, 12]:
                        res = abs(b - r * a)
                        if res < best_res:
                            best_r = r
                            best_res = res
                    residual = b - best_r * a
                    # 2 bit ratio + delta_a + residual
                    bytes_adaptive_ratio += 0.25 + varint_len(zigzag_encode(delta_a)) + varint_len(zigzag_encode(residual))
                else:
                    bytes_adaptive_ratio += 0.25 + varint_len(zigzag_encode(delta_a)) + varint_len(b)

                # 方案3：联合差分（delta_a + delta_b）
                bytes_joint_delta += varint_len(zigzag_encode(delta_a)) + varint_len(zigzag_encode(delta_b))

    print(f"="*70)
    print(f"文件: {vcf_path}")
    print(f"Type 1 样本数: {total_samples:,}")
    print(f"="*70)

    print(f"\n[delta_b 分布]")
    for k, v in sorted(delta_b_dist.items()):
        print(f"  {k}: {v:,} ({100*v/total_samples:.1f}%)")

    print(f"\n[编码字节数对比]")
    print(f"  基线 (a + b):              {bytes_baseline:,} bytes ({bytes_baseline/total_samples:.2f} B/sample)")

    for r in [10, 11, 12]:
        b = bytes_fixed_ratio[r]
        print(f"  固定 ratio={r}:            {b:,} bytes ({b/total_samples:.2f} B/sample), 节省 {100*(1-b/bytes_baseline):.1f}%")

    print(f"  自适应 ratio:              {int(bytes_adaptive_ratio):,} bytes ({bytes_adaptive_ratio/total_samples:.2f} B/sample), 节省 {100*(1-bytes_adaptive_ratio/bytes_baseline):.1f}%")

    print(f"  联合差分 (Δa + Δb):        {bytes_joint_delta:,} bytes ({bytes_joint_delta/total_samples:.2f} B/sample), 节省 {100*(1-bytes_joint_delta/bytes_baseline):.1f}%")

    # 最优方案
    best_bytes = min(bytes_fixed_ratio[11], bytes_adaptive_ratio, bytes_joint_delta)
    print(f"\n  ★ 最优方案节省: {100*(1-best_bytes/bytes_baseline):.1f}%")
    print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf.gz>")
        sys.exit(1)

    vcf = sys.argv[1]
    max_var = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    max_samp = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    analyze_best_codec(vcf, max_var, max_samp)
