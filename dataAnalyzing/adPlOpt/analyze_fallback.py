#!/usr/bin/env python3
"""
分析 PL codec_id=4 兜底路径的特征，寻找进一步优化空间
"""

import sys
import gzip
from collections import defaultdict

def zigzag_encode(v):
    return (v << 1) ^ (v >> 31) if v >= 0 else ((-v) << 1) - 1

def analyze_fallback(vcf_path, max_variants=30000, max_samples=500):
    """分析兜底路径的特征"""

    THRESHOLD = 63

    # 兜底原因统计
    fallback_reasons = defaultdict(int)
    fallback_delta_a = []
    fallback_residual = []
    fallback_a_values = []
    fallback_b_values = []

    # 位打包潜力
    bitpack_5_5 = 0  # |delta_a| <= 15 && |residual| <= 15
    bitpack_6_6 = 0  # |delta_a| <= 31 && |residual| <= 31
    bitpack_7_7 = 0  # |delta_a| <= 63 && |residual| <= 63

    total_type1 = 0
    prev_a = {}

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

                p0, p1, p2 = pl_vals
                if p0 < 0 or p1 < 0 or p2 < 0:
                    continue

                min_v = min(p0, p1, p2)
                if min_v != 0:
                    continue

                min_pos = 0
                if p1 == 0: min_pos = 1
                elif p2 == 0: min_pos = 2

                if min_pos == 0:
                    a, b = p1, p2
                elif min_pos == 1:
                    a, b = p0, p2
                else:
                    a, b = p1, p0

                if a == 0 and b == 0:
                    continue
                if a > 0 and b == 15 * a:
                    prev_a[s_idx] = a
                    continue
                if a <= 0:
                    continue

                total_type1 += 1
                delta_a = a - prev_a.get(s_idx, a)
                prev_a[s_idx] = a
                residual = b - 11 * a

                # 位打包潜力
                if abs(delta_a) <= 15 and abs(residual) <= 15:
                    bitpack_5_5 += 1
                if abs(delta_a) <= 31 and abs(residual) <= 31:
                    bitpack_6_6 += 1
                if abs(delta_a) <= 63 and abs(residual) <= 63:
                    bitpack_7_7 += 1

                # 兜底分析
                if abs(delta_a) > THRESHOLD or abs(residual) > THRESHOLD:
                    fallback_delta_a.append(delta_a)
                    fallback_residual.append(residual)
                    fallback_a_values.append(a)
                    fallback_b_values.append(b)

                    if abs(delta_a) > THRESHOLD and abs(residual) > THRESHOLD:
                        fallback_reasons['both'] += 1
                    elif abs(delta_a) > THRESHOLD:
                        fallback_reasons['delta_a'] += 1
                    else:
                        fallback_reasons['residual'] += 1

    print(f"=" * 70)
    print(f"文件: {vcf_path}")
    print(f"Type 1 样本数: {total_type1:,}")
    print(f"=" * 70)

    fallback_total = len(fallback_delta_a)
    print(f"\n[兜底路径分析] ({fallback_total:,} 样本, {100*fallback_total/total_type1:.1f}%)")

    print(f"\n  兜底原因分布:")
    for reason in ['delta_a', 'residual', 'both']:
        cnt = fallback_reasons[reason]
        print(f"    {reason}: {cnt:,} ({100*cnt/fallback_total:.1f}%)")

    print(f"\n  兜底样本的 delta_a 分布:")
    buckets = defaultdict(int)
    for d in fallback_delta_a:
        ad = abs(d)
        if ad <= 63: buckets['<=63'] += 1
        elif ad <= 127: buckets['64-127'] += 1
        elif ad <= 255: buckets['128-255'] += 1
        else: buckets['>255'] += 1
    for k in ['<=63', '64-127', '128-255', '>255']:
        if k in buckets:
            print(f"    |Δa| {k}: {buckets[k]:,} ({100*buckets[k]/fallback_total:.1f}%)")

    print(f"\n  兜底样本的 residual 分布:")
    buckets = defaultdict(int)
    for r in fallback_residual:
        ar = abs(r)
        if ar <= 63: buckets['<=63'] += 1
        elif ar <= 127: buckets['64-127'] += 1
        elif ar <= 255: buckets['128-255'] += 1
        else: buckets['>255'] += 1
    for k in ['<=63', '64-127', '128-255', '>255']:
        if k in buckets:
            print(f"    |res| {k}: {buckets[k]:,} ({100*buckets[k]/fallback_total:.1f}%)")

    print(f"\n  兜底样本的 a 值分布:")
    buckets = defaultdict(int)
    for a in fallback_a_values:
        if a < 10: buckets['<10'] += 1
        elif a < 30: buckets['10-30'] += 1
        elif a < 50: buckets['30-50'] += 1
        elif a < 100: buckets['50-100'] += 1
        else: buckets['>=100'] += 1
    for k in ['<10', '10-30', '30-50', '50-100', '>=100']:
        if k in buckets:
            print(f"    a {k}: {buckets[k]:,} ({100*buckets[k]/fallback_total:.1f}%)")

    print(f"\n[位打包潜力]")
    print(f"  5+5 bit (|v|<=15): {bitpack_5_5:,} ({100*bitpack_5_5/total_type1:.1f}%)")
    print(f"  6+6 bit (|v|<=31): {bitpack_6_6:,} ({100*bitpack_6_6/total_type1:.1f}%)")
    print(f"  7+7 bit (|v|<=63): {bitpack_7_7:,} ({100*bitpack_7_7/total_type1:.1f}%)")

    # 字节数估算
    print(f"\n[位打包 vs 当前编码]")
    current_bytes = bitpack_7_7 * 2 + (total_type1 - bitpack_7_7) * 4  # 2B tip10, 4B fallback
    bitpack_bytes_5 = bitpack_5_5 * 1.25 + (total_type1 - bitpack_5_5) * 4
    bitpack_bytes_6 = bitpack_6_6 * 1.5 + (total_type1 - bitpack_6_6) * 4
    print(f"  当前 (7+7 bit varint): {current_bytes/1e6:.2f} MB")
    print(f"  位打包 5+5 bit:        {bitpack_bytes_5/1e6:.2f} MB, 节省 {100*(1-bitpack_bytes_5/current_bytes):.1f}%")
    print(f"  位打包 6+6 bit:        {bitpack_bytes_6/1e6:.2f} MB, 节省 {100*(1-bitpack_bytes_6/current_bytes):.1f}%")

    # 进一步优化建议
    print(f"\n[优化建议]")
    if bitpack_5_5 / total_type1 > 0.7:
        print(f"  ★ 建议使用 5+5 bit 位打包 (覆盖 {100*bitpack_5_5/total_type1:.0f}%)")
    elif bitpack_6_6 / total_type1 > 0.8:
        print(f"  ★ 建议使用 6+6 bit 位打包 (覆盖 {100*bitpack_6_6/total_type1:.0f}%)")

    print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf.gz>")
        sys.exit(1)

    vcf = sys.argv[1]
    max_var = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
    max_samp = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    analyze_fallback(vcf, max_var, max_samp)
