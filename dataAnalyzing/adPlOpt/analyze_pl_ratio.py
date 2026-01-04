#!/usr/bin/env python3
"""分析 PL 字段的 b/a 比值分布，找出最优编码策略"""

import sys
import gzip
from collections import defaultdict

def analyze_pl(vcf_path, max_variants=20000, max_samples=500):
    """分析 PL 的 (a, b) 值分布"""

    # 统计
    type1_count = 0
    a_values = []
    b_values = []
    ratio_distribution = defaultdict(int)
    residual_distribution = defaultdict(int)

    # 字节统计
    bytes_current = 0  # 当前编码
    bytes_ratio12 = 0  # 比值编码
    bytes_delta_a = 0  # a 差分编码

    # 差分统计
    prev_a = {}  # sample_id -> prev_a
    a_deltas = []

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

                # 检查 Type 1: [0, a, b]
                if pl_vals[0] == 0 and pl_vals[1] >= 0 and pl_vals[2] >= 0:
                    a, b = pl_vals[1], pl_vals[2]

                    # 跳过 Type 0 (全零) 和 Type 2 (b=15*a)
                    if a == 0 and b == 0:
                        continue
                    if a > 0 and b == 15 * a:
                        continue

                    type1_count += 1
                    a_values.append(a)
                    b_values.append(b)

                    # 计算比值
                    if a > 0:
                        ratio = b / a
                        ratio_int = round(ratio)
                        ratio_distribution[ratio_int] += 1

                        # 计算不同比值的残差
                        for r in [10, 11, 12, 13, 14, 15]:
                            res = b - r * a
                            if abs(res) <= 7:
                                residual_distribution[(r, res)] += 1

                    # 当前编码字节数
                    def varint_len(v):
                        if v < 128: return 1
                        if v < 16384: return 2
                        if v < 2097152: return 3
                        return 4

                    bytes_current += varint_len(a) + varint_len(b)

                    # a 差分统计
                    if s_idx in prev_a:
                        delta = a - prev_a[s_idx]
                        a_deltas.append(delta)
                    prev_a[s_idx] = a

    print(f"="*70)
    print(f"文件: {vcf_path}")
    print(f"分析位点数: {variant_count}, 样本数: {sample_count}")
    print(f"Type 1 样本数: {type1_count}")
    print(f"="*70)

    if not a_values:
        print("没有找到 Type 1 数据")
        return

    # a 值分布
    print(f"\n[a 值分布]")
    a_buckets = defaultdict(int)
    for a in a_values:
        if a < 10: a_buckets['<10'] += 1
        elif a < 30: a_buckets['10-30'] += 1
        elif a < 50: a_buckets['30-50'] += 1
        elif a < 128: a_buckets['50-128'] += 1
        else: a_buckets['>=128'] += 1
    for k, v in sorted(a_buckets.items()):
        print(f"  {k}: {v} ({100*v/len(a_values):.1f}%)")

    # b 值分布
    print(f"\n[b 值分布]")
    b_buckets = defaultdict(int)
    for b in b_values:
        if b < 50: b_buckets['<50'] += 1
        elif b < 128: b_buckets['50-128'] += 1
        elif b < 256: b_buckets['128-256'] += 1
        elif b < 512: b_buckets['256-512'] += 1
        else: b_buckets['>=512'] += 1
    for k, v in sorted(b_buckets.items()):
        print(f"  {k}: {v} ({100*v/len(b_values):.1f}%)")

    # b/a 比值分布
    print(f"\n[b/a 比值分布 (整数)]")
    total = sum(ratio_distribution.values())
    for ratio in sorted(ratio_distribution.keys()):
        cnt = ratio_distribution[ratio]
        if cnt / total > 0.01:  # 只显示 >1% 的
            print(f"  ratio={ratio}: {cnt} ({100*cnt/total:.1f}%)")

    # 残差分布
    print(f"\n[最优比值及残差分布]")
    best_ratio = 12
    best_coverage = 0
    for r in [10, 11, 12, 13, 14]:
        coverage = sum(residual_distribution[(r, res)] for res in range(-7, 8))
        pct = 100 * coverage / type1_count if type1_count > 0 else 0
        print(f"  ratio={r}, |residual|<=7: {coverage} ({pct:.1f}%)")
        if coverage > best_coverage:
            best_coverage = coverage
            best_ratio = r

    print(f"\n  最优比值: {best_ratio}, 覆盖率: {100*best_coverage/type1_count:.1f}%")

    # a 差分分布
    print(f"\n[a 值位点间差分分布]")
    if a_deltas:
        delta_buckets = defaultdict(int)
        for d in a_deltas:
            if abs(d) <= 3: delta_buckets['|d|<=3'] += 1
            elif abs(d) <= 7: delta_buckets['4<=|d|<=7'] += 1
            elif abs(d) <= 15: delta_buckets['8<=|d|<=15'] += 1
            else: delta_buckets['|d|>15'] += 1
        for k, v in sorted(delta_buckets.items()):
            print(f"  {k}: {v} ({100*v/len(a_deltas):.1f}%)")

    # 字节数对比估算
    print(f"\n[编码字节数估算]")
    print(f"  当前编码 (a + b): {bytes_current} bytes")

    # 比值编码估算
    bytes_ratio = 0
    for i, (a, b) in enumerate(zip(a_values, b_values)):
        if a > 0:
            res = b - best_ratio * a
            if abs(res) <= 7:
                # a + 4bit残差
                bytes_ratio += (1 if a < 128 else 2) + 0.5
            else:
                # a + b 原始
                bytes_ratio += (1 if a < 128 else 2) + (1 if b < 128 else 2)
        else:
            bytes_ratio += 1 + (1 if b < 128 else 2)
    print(f"  比值编码 (ratio={best_ratio}): {int(bytes_ratio)} bytes, 节省 {100*(1-bytes_ratio/bytes_current):.1f}%")

    # a 差分编码估算
    if a_deltas:
        bytes_delta = len(a_values)  # 假设每个 a 差分用 1 字节
        bytes_delta += sum(1 if b < 128 else 2 for b in b_values)
        print(f"  a 差分编码: {bytes_delta} bytes, 节省 {100*(1-bytes_delta/bytes_current):.1f}%")

    print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf.gz> [max_variants] [max_samples]")
        sys.exit(1)

    vcf = sys.argv[1]
    max_var = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    max_samp = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    analyze_pl(vcf, max_var, max_samp)
