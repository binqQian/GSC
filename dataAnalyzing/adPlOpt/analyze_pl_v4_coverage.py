#!/usr/bin/env python3
"""
分析 PL codec_id=4 的覆盖情况，找出进一步优化空间
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

def analyze_coverage(vcf_path, max_variants=20000, max_samples=500):
    """分析各种编码路径的覆盖情况"""

    # 当前实现的阈值
    THRESHOLD_63 = 63
    THRESHOLD_127 = 127

    # 统计
    total_type1 = 0
    total_type2 = 0
    total_type0 = 0
    total_type3 = 0

    # Tip 覆盖统计
    tip_10_63 = 0      # |delta_a|<=63 && |residual|<=63
    tip_10_127 = 0     # |delta_a|<=127 && |residual|<=127
    tip_11_tag1 = 0    # Type 1 超阈值

    # 残差分布（多种 ratio）
    residual_by_ratio = defaultdict(list)

    # delta_a 分布
    delta_a_dist = defaultdict(int)

    # b 值差分
    prev_b = {}
    delta_b_small = 0  # |delta_b| <= 31

    # 字节数统计
    bytes_current = 0       # 当前实现 (阈值 63)
    bytes_threshold_127 = 0 # 阈值 127
    bytes_multi_ratio = 0   # 多 ratio 支持
    bytes_delta_b = 0       # delta_b 编码

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

                # 归一化处理
                p0, p1, p2 = pl_vals
                if p0 < 0 or p1 < 0 or p2 < 0:
                    continue

                # min-pos normalization
                min_v = min(p0, p1, p2)
                if min_v != 0:
                    continue  # 跳过无效数据

                min_pos = 0
                if p1 == 0: min_pos = 1
                elif p2 == 0: min_pos = 2

                # 归一化为 [0, a, b]
                if min_pos == 0:
                    a, b = p1, p2
                elif min_pos == 1:
                    a, b = p0, p2
                else:
                    a, b = p1, p0

                # 分类
                if a == 0 and b == 0:
                    total_type0 += 1
                    continue
                elif a > 0 and b == 15 * a:
                    total_type2 += 1
                    # Type 2 使用 delta_a 编码
                    delta_a = a - prev_a.get(s_idx, a)
                    prev_a[s_idx] = a
                    continue
                elif a > 0:
                    total_type1 += 1
                else:
                    total_type3 += 1
                    continue

                # Type 1 分析
                delta_a = a - prev_a.get(s_idx, a)
                prev_a[s_idx] = a

                # delta_b
                delta_b = b - prev_b.get(s_idx, b)
                prev_b[s_idx] = b
                if abs(delta_b) <= 31:
                    delta_b_small += 1

                # 多种 ratio 的残差
                for r in [10, 11, 12]:
                    res = b - r * a
                    residual_by_ratio[r].append(res)

                # 当前实现 (阈值 63)
                residual_11 = b - 11 * a
                if abs(delta_a) <= THRESHOLD_63 and abs(residual_11) <= THRESHOLD_63:
                    tip_10_63 += 1
                    bytes_current += 2  # delta_a + residual
                else:
                    tip_11_tag1 += 1
                    bytes_current += varint_len(a) + varint_len(b) + 1  # tag + a + b

                # 阈值 127
                if abs(delta_a) <= THRESHOLD_127 and abs(residual_11) <= THRESHOLD_127:
                    tip_10_127 += 1
                    bytes_threshold_127 += 2
                else:
                    bytes_threshold_127 += varint_len(a) + varint_len(b) + 1

                # 多 ratio 支持
                best_res = abs(residual_11)
                best_ratio = 11
                for r in [10, 12]:
                    res = abs(b - r * a)
                    if res < best_res:
                        best_res = res
                        best_ratio = r
                if abs(delta_a) <= THRESHOLD_127 and best_res <= THRESHOLD_127:
                    bytes_multi_ratio += 2.25  # delta_a + residual + 2bit ratio
                else:
                    bytes_multi_ratio += varint_len(a) + varint_len(b) + 1

                # delta_b 编码
                if abs(delta_a) <= THRESHOLD_127 and abs(delta_b) <= THRESHOLD_127:
                    bytes_delta_b += 2  # delta_a + delta_b
                else:
                    bytes_delta_b += varint_len(a) + varint_len(b) + 1

                # delta_a 分布
                if abs(delta_a) <= 3:
                    delta_a_dist['<=3'] += 1
                elif abs(delta_a) <= 7:
                    delta_a_dist['4-7'] += 1
                elif abs(delta_a) <= 15:
                    delta_a_dist['8-15'] += 1
                elif abs(delta_a) <= 31:
                    delta_a_dist['16-31'] += 1
                elif abs(delta_a) <= 63:
                    delta_a_dist['32-63'] += 1
                elif abs(delta_a) <= 127:
                    delta_a_dist['64-127'] += 1
                else:
                    delta_a_dist['>127'] += 1

    # 输出结果
    print(f"=" * 70)
    print(f"文件: {vcf_path}")
    print(f"分析位点数: {variant_count}, 样本数: {sample_count}")
    print(f"=" * 70)

    total = total_type0 + total_type1 + total_type2 + total_type3
    print(f"\n[PL 类型分布]")
    print(f"  Type 0 (全零):     {total_type0:,} ({100*total_type0/total:.1f}%)")
    print(f"  Type 1 ([0,a,b]):  {total_type1:,} ({100*total_type1/total:.1f}%)")
    print(f"  Type 2 (b=15*a):   {total_type2:,} ({100*total_type2/total:.1f}%)")
    print(f"  Type 3 (其他):     {total_type3:,} ({100*total_type3/total:.1f}%)")

    print(f"\n[Type 1 编码覆盖率]")
    print(f"  Tip 10 (阈值 63):  {tip_10_63:,} ({100*tip_10_63/total_type1:.1f}%)")
    print(f"  Tip 10 (阈值 127): {tip_10_127:,} ({100*tip_10_127/total_type1:.1f}%)")
    print(f"  Tip 11 兜底:       {tip_11_tag1:,} ({100*tip_11_tag1/total_type1:.1f}%)")

    print(f"\n[delta_a 分布]")
    for k in ['<=3', '4-7', '8-15', '16-31', '32-63', '64-127', '>127']:
        if k in delta_a_dist:
            print(f"  |Δa| {k}: {delta_a_dist[k]:,} ({100*delta_a_dist[k]/total_type1:.1f}%)")

    print(f"\n[残差分布 (各 ratio)]")
    for r in [10, 11, 12]:
        residuals = residual_by_ratio[r]
        small = sum(1 for x in residuals if abs(x) <= 63)
        medium = sum(1 for x in residuals if abs(x) <= 127)
        print(f"  ratio={r}: |res|<=63: {100*small/len(residuals):.1f}%, |res|<=127: {100*medium/len(residuals):.1f}%")

    print(f"\n[delta_b 潜力]")
    print(f"  |Δb| <= 31: {delta_b_small:,} ({100*delta_b_small/total_type1:.1f}%)")

    # 简化计算
    baseline_bytes = total_type1 * 2.7  # 平均 2.7 字节

    print(f"\n[Type 1 编码字节数对比]")
    print(f"  基线 (a + b):           ~{baseline_bytes/1e6:.2f} MB ({baseline_bytes/total_type1:.2f} B/sample)")
    print(f"  当前实现 (阈值 63):     ~{bytes_current/1e6:.2f} MB ({bytes_current/total_type1:.2f} B/sample), 节省 {100*(1-bytes_current/baseline_bytes):.1f}%")
    print(f"  阈值 127:               ~{bytes_threshold_127/1e6:.2f} MB ({bytes_threshold_127/total_type1:.2f} B/sample), 节省 {100*(1-bytes_threshold_127/baseline_bytes):.1f}%")
    print(f"  多 ratio (10/11/12):    ~{bytes_multi_ratio/1e6:.2f} MB ({bytes_multi_ratio/total_type1:.2f} B/sample), 节省 {100*(1-bytes_multi_ratio/baseline_bytes):.1f}%")
    print(f"  delta_b 编码:           ~{bytes_delta_b/1e6:.2f} MB ({bytes_delta_b/total_type1:.2f} B/sample), 节省 {100*(1-bytes_delta_b/baseline_bytes):.1f}%")

    # 额外优化潜力
    print(f"\n[额外优化潜力]")
    extra_63_to_127 = tip_10_127 - tip_10_63
    print(f"  阈值 63→127: +{extra_63_to_127:,} 样本 (+{100*extra_63_to_127/total_type1:.1f}%)")
    print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf.gz>")
        sys.exit(1)

    vcf = sys.argv[1]
    max_var = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    max_samp = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    analyze_coverage(vcf, max_var, max_samp)
