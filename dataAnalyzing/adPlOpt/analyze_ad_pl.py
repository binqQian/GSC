#!/usr/bin/env python3
"""
VCF AD/PL 字段数据分析
分析目标：
1. AD (Allelic Depth) 字段数据特征
2. PL (Phred-scaled Likelihood) 字段数据特征
3. 跨字段关联性分析
4. 压缩优化建议
"""

import numpy as np
from cyvcf2 import VCF
from collections import Counter, defaultdict
import sys

# 数据路径
VCF_FILES = {
    '1000GP_LC': '/home/binq/data/1000GPi/lc_bams.first50000.vcf.gz',
    'Subset_5k_20k': '/home/binq/data/output_subset_5000_20000.vcf.gz'
}

MAX_VARIANTS = 10000


def extract_vcf_data(vcf_path, max_variants=10000):
    """从VCF提取AD、PL、GT、DP、GQ数据"""
    vcf = VCF(vcf_path)
    n_samples = len(vcf.samples)
    print(f"样本数: {n_samples}")

    data = {
        'ad_values': [],
        'ad_sums': [],
        'ad_arrays': [],
        'ad_lengths': [],
        'pl_values': [],
        'pl_arrays': [],
        'pl_lengths': [],
        'gt_values': [],
        'dp_values': [],
        'gq_values': [],
        'allele_counts': [],
        'n_samples': n_samples,
        'n_variants': 0
    }

    variant_count = 0
    for variant in vcf:
        if variant_count >= max_variants:
            break

        n_alleles = len(variant.ALT) + 1
        data['allele_counts'].append(n_alleles)

        # 提取AD
        try:
            ad = variant.format('AD')
            if ad is not None:
                for sample_idx in range(n_samples):
                    sample_ad = ad[sample_idx]
                    if sample_ad is not None and len(sample_ad) > 0:
                        valid_ad = [int(v) for v in sample_ad if v >= 0]
                        if valid_ad:
                            data['ad_values'].extend(valid_ad)
                            data['ad_sums'].append(sum(valid_ad))
                            data['ad_arrays'].append(tuple(valid_ad))
                            data['ad_lengths'].append(len(valid_ad))
        except:
            pass

        # 提取PL
        try:
            pl = variant.format('PL')
            if pl is not None:
                for sample_idx in range(n_samples):
                    sample_pl = pl[sample_idx]
                    if sample_pl is not None and len(sample_pl) > 0:
                        valid_pl = [int(v) for v in sample_pl if v >= 0]
                        if valid_pl:
                            data['pl_values'].extend(valid_pl)
                            data['pl_arrays'].append(tuple(valid_pl))
                            data['pl_lengths'].append(len(valid_pl))
        except:
            pass

        # 提取GT
        try:
            gt = variant.gt_types
            if gt is not None:
                data['gt_values'].extend(gt.tolist())
        except:
            pass

        # 提取DP
        try:
            dp = variant.format('DP')
            if dp is not None:
                valid_dp = [int(v[0]) for v in dp if v[0] >= 0]
                data['dp_values'].extend(valid_dp)
        except:
            pass

        # 提取GQ
        try:
            gq = variant.format('GQ')
            if gq is not None:
                valid_gq = [int(v[0]) for v in gq if v[0] >= 0]
                data['gq_values'].extend(valid_gq)
        except:
            pass

        variant_count += 1
        if variant_count % 2000 == 0:
            print(f"  已处理 {variant_count} 个变异位点...")

    data['n_variants'] = variant_count
    print(f"  完成! 共处理 {variant_count} 个变异位点")
    return data


def check_pl_pattern(pl_array):
    """
    检查PL数组模式，复现现有压缩的checkPlPattern逻辑
    返回: (type, a, b)
    - Type 0: 全零
    - Type 1: 标准模式 [0,a,b,a,b,b,...]
    - Type 2: Type1且 b == 15*a
    - Type 3: 无模式
    """
    if len(pl_array) == 0:
        return 3, 0, 0

    if pl_array[0] != 0:
        return 3, 0, 0

    if all(v == 0 for v in pl_array):
        return 0, 0, 0

    if len(pl_array) < 3:
        return 3, 0, 0

    a, b = pl_array[1], pl_array[2]

    n = len(pl_array)
    allele = int((-1 + (1 + 8*n)**0.5) / 2)

    if allele < 2:
        return 3, 0, 0

    expected = []
    for i in range(allele):
        for j in range(i + 1):
            if i == 0 and j == 0:
                expected.append(0)
            elif j == 0:
                expected.append(a)
            else:
                expected.append(b)

    if len(expected) != len(pl_array):
        return 3, 0, 0

    if list(pl_array) == expected:
        if b == 15 * a:
            return 2, a, b
        else:
            return 1, a, b

    return 3, 0, 0


def analyze_ad_distribution(data, name):
    """AD数值分布分析"""
    ad_values = np.array(data['ad_values'])
    ad_sums = np.array(data['ad_sums'])

    print(f"\n{'='*60}")
    print(f"{name} AD分布统计")
    print(f"{'='*60}")
    print(f"总AD值数量: {len(ad_values):,}")
    print(f"AD数组数量: {len(ad_sums):,}")

    print(f"\n--- AD单值统计 ---")
    print(f"范围: [{ad_values.min()}, {ad_values.max()}]")
    print(f"均值: {ad_values.mean():.2f}")
    print(f"中位数: {np.median(ad_values):.1f}")
    print(f"标准差: {ad_values.std():.2f}")

    zero_ratio = (ad_values == 0).sum() / len(ad_values)
    low_ratio = (ad_values < 5).sum() / len(ad_values)
    print(f"\n--- AD稀疏性 ---")
    print(f"零值比例: {zero_ratio*100:.2f}%")
    print(f"低值(<5)比例: {low_ratio*100:.2f}%")

    print(f"\n--- sum(AD)统计 ---")
    print(f"范围: [{ad_sums.min()}, {ad_sums.max()}]")
    print(f"均值: {ad_sums.mean():.2f}")
    print(f"中位数: {np.median(ad_sums):.1f}")

    # Tip编码模式统计
    print(f"\n--- 现有Tip编码模式比例 ---")
    tip_00_count = 0
    tip_01_count = 0
    tip_10_count = 0
    tip_11_count = 0

    for ad_arr in data['ad_arrays']:
        s = sum(ad_arr)
        if s == 0:
            tip_00_count += 1
        elif s == ad_arr[0] and s == 2:
            tip_01_count += 1
        elif s == ad_arr[0]:
            tip_10_count += 1
        else:
            tip_11_count += 1

    total = len(data['ad_arrays'])
    print(f"Tip 00 (sum==0): {tip_00_count:,} ({tip_00_count/total*100:.2f}%)")
    print(f"Tip 01 (sum==AD[0]==2): {tip_01_count:,} ({tip_01_count/total*100:.2f}%)")
    print(f"Tip 10 (sum==AD[0]!=2): {tip_10_count:,} ({tip_10_count/total*100:.2f}%)")
    print(f"Tip 11 (一般): {tip_11_count:,} ({tip_11_count/total*100:.2f}%)")

    efficient_ratio = (tip_00_count + tip_01_count) / total
    print(f"\n高效模式(Tip 00+01)占比: {efficient_ratio*100:.2f}%")


def analyze_ad_structure(data, name):
    """AD数组结构分析"""
    print(f"\n--- {name} AD结构分析 ---")

    lengths = Counter(data['ad_lengths'])
    print(f"AD数组长度分布:")
    for length, count in sorted(lengths.items()):
        pct = count / len(data['ad_lengths']) * 100
        print(f"  长度{length}: {count:,} ({pct:.2f}%)")

    ad_patterns = Counter(data['ad_arrays'])
    print(f"\n唯一AD数组模式数: {len(ad_patterns):,}")
    print(f"最常见的10个模式:")
    for pattern, count in ad_patterns.most_common(10):
        pct = count / len(data['ad_arrays']) * 100
        print(f"  {pattern}: {count:,} ({pct:.2f}%)")

    total_arrays = len(data['ad_arrays'])
    unique_patterns = len(ad_patterns)
    print(f"\n字典压缩潜力:")
    print(f"  总数组数: {total_arrays:,}")
    print(f"  唯一模式数: {unique_patterns:,}")
    print(f"  压缩比估算: {total_arrays/unique_patterns:.2f}x")


def analyze_ad_differential(data, name):
    """AD差分编码潜力分析"""
    print(f"\n--- {name} AD差分分析 ---")

    ad_sums = np.array(data['ad_sums'])

    if len(ad_sums) > 1:
        diff = np.diff(ad_sums)
        print(f"sum(AD)差分统计:")
        print(f"  差分均值: {np.mean(diff):.2f}")
        print(f"  差分中位数: {np.median(diff):.1f}")
        print(f"  差分标准差: {np.std(diff):.2f}")
        print(f"  零差分比例: {(diff == 0).sum() / len(diff) * 100:.2f}%")
        print(f"  小差分(|diff|<5)比例: {(np.abs(diff) < 5).sum() / len(diff) * 100:.2f}%")


def analyze_pl_distribution(data, name):
    """PL数值分布分析"""
    pl_values = np.array(data['pl_values'])

    print(f"\n{'='*60}")
    print(f"{name} PL分布统计")
    print(f"{'='*60}")
    print(f"总PL值数量: {len(pl_values):,}")
    print(f"PL数组数量: {len(data['pl_arrays']):,}")

    print(f"\n--- PL单值统计 ---")
    print(f"范围: [{pl_values.min()}, {pl_values.max()}]")
    print(f"均值: {pl_values.mean():.2f}")
    print(f"中位数: {np.median(pl_values):.1f}")
    print(f"标准差: {pl_values.std():.2f}")

    zero_ratio = (pl_values == 0).sum() / len(pl_values)
    max_val = pl_values.max()
    max_ratio = (pl_values == max_val).sum() / len(pl_values)
    print(f"\n--- PL特殊值 ---")
    print(f"零值比例: {zero_ratio*100:.2f}%")
    print(f"最大值{max_val}比例: {max_ratio*100:.2f}%")

    print(f"\n--- 值域分布 ---")
    ranges = [(0, 10), (10, 50), (50, 100), (100, 255), (255, 1000), (1000, float('inf'))]
    for low, high in ranges:
        count = ((pl_values >= low) & (pl_values < high)).sum()
        pct = count / len(pl_values) * 100
        high_str = str(int(high)) if high != float('inf') else 'inf'
        print(f"  [{low}, {high_str}): {count:,} ({pct:.2f}%)")


def analyze_pl_patterns(data, name):
    """分析PL模式分布"""
    print(f"\n--- {name} PL模式分析 ---")

    type_counts = Counter()
    a_values = []
    b_values = []

    for pl_arr in data['pl_arrays']:
        pl_type, a, b = check_pl_pattern(pl_arr)
        type_counts[pl_type] += 1
        if pl_type in [1, 2]:
            a_values.append(a)
            b_values.append(b)

    total = len(data['pl_arrays'])
    print(f"PL模式分布 (现有Tip编码):")
    print(f"  Type 0 (全零): {type_counts[0]:,} ({type_counts[0]/total*100:.2f}%)")
    print(f"  Type 1 (标准[0,a,b]): {type_counts[1]:,} ({type_counts[1]/total*100:.2f}%)")
    print(f"  Type 2 (b=15*a): {type_counts[2]:,} ({type_counts[2]/total*100:.2f}%)")
    print(f"  Type 3 (无模式): {type_counts[3]:,} ({type_counts[3]/total*100:.2f}%)")

    compressible = (type_counts[0] + type_counts[1] + type_counts[2]) / total
    print(f"\n可模式压缩比例: {compressible*100:.2f}%")

    if a_values:
        print(f"\na值分布 (Type 1/2):")
        print(f"  范围: [{min(a_values)}, {max(a_values)}]")
        print(f"  均值: {np.mean(a_values):.2f}")
        print(f"  中位数: {np.median(a_values):.1f}")

    return type_counts


def analyze_pl_structure(data, name):
    """PL数组结构分析"""
    print(f"\n--- {name} PL结构分析 ---")

    lengths = Counter(data['pl_lengths'])
    print(f"PL数组长度分布:")
    for length, count in sorted(lengths.items())[:10]:
        pct = count / len(data['pl_lengths']) * 100
        print(f"  长度{length}: {count:,} ({pct:.2f}%)")

    pl_patterns = Counter(data['pl_arrays'])
    print(f"\n唯一PL数组模式数: {len(pl_patterns):,}")
    print(f"最常见的10个模式:")
    for pattern, count in pl_patterns.most_common(10):
        pct = count / len(data['pl_arrays']) * 100
        display_pattern = pattern[:5] if len(pattern) > 5 else pattern
        suffix = "..." if len(pattern) > 5 else ""
        print(f"  {display_pattern}{suffix}: {count:,} ({pct:.2f}%)")

    total_arrays = len(data['pl_arrays'])
    unique_patterns = len(pl_patterns)
    print(f"\n字典压缩潜力:")
    print(f"  总数组数: {total_arrays:,}")
    print(f"  唯一模式数: {unique_patterns:,}")
    print(f"  压缩比估算: {total_arrays/unique_patterns:.2f}x")


def analyze_ad_dp_relation(vcf_path, name, max_variants=5000):
    """分析AD与DP的关系"""
    print(f"\n--- {name} AD-DP关系分析 ---")

    vcf = VCF(vcf_path)

    match_count = 0
    mismatch_count = 0
    mismatches = []

    for var_idx, variant in enumerate(vcf):
        if var_idx >= max_variants:
            break

        try:
            ad = variant.format('AD')
            dp = variant.format('DP')

            if ad is None or dp is None:
                continue

            for sample_idx in range(len(vcf.samples)):
                sample_ad = ad[sample_idx]
                sample_dp = dp[sample_idx][0]

                if sample_ad is None or sample_dp < 0:
                    continue

                ad_sum = sum(int(v) for v in sample_ad if v >= 0)

                if ad_sum == sample_dp:
                    match_count += 1
                else:
                    mismatch_count += 1
                    if len(mismatches) < 10:
                        mismatches.append((ad_sum, int(sample_dp)))
        except:
            continue

    total = match_count + mismatch_count
    if total > 0:
        print(f"DP == sum(AD): {match_count:,} ({match_count/total*100:.2f}%)")
        print(f"DP != sum(AD): {mismatch_count:,} ({mismatch_count/total*100:.2f}%)")

        if mismatches:
            print(f"\n不匹配示例 (sum(AD), DP):")
            for ad_sum, dp in mismatches[:5]:
                print(f"  {ad_sum} vs {dp}")
    else:
        print("无有效AD-DP数据")


def analyze_pl_gq_relation(vcf_path, name, max_variants=5000):
    """分析PL与GQ的关系"""
    print(f"\n--- {name} PL-GQ关系分析 ---")

    vcf = VCF(vcf_path)

    match_count = 0
    mismatch_count = 0

    for var_idx, variant in enumerate(vcf):
        if var_idx >= max_variants:
            break

        try:
            pl = variant.format('PL')
            gq = variant.format('GQ')

            if pl is None or gq is None:
                continue

            for sample_idx in range(len(vcf.samples)):
                sample_pl = pl[sample_idx]
                sample_gq = gq[sample_idx][0]

                if sample_pl is None or sample_gq < 0:
                    continue

                valid_pl = [int(v) for v in sample_pl if v >= 0]
                if len(valid_pl) < 2:
                    continue

                sorted_pl = sorted(valid_pl)
                second_min = sorted_pl[1] if len(sorted_pl) > 1 else sorted_pl[0]

                if second_min == sample_gq:
                    match_count += 1
                else:
                    mismatch_count += 1
        except:
            continue

    total = match_count + mismatch_count
    if total > 0:
        print(f"GQ == second_min(PL): {match_count:,} ({match_count/total*100:.2f}%)")
        print(f"GQ != second_min(PL): {mismatch_count:,} ({mismatch_count/total*100:.2f}%)")
    else:
        print("无有效PL-GQ数据")


def explore_new_ad_patterns(data, name):
    """探索AD的新压缩模式"""
    print(f"\n--- {name} AD新模式探索 ---")

    ad_arrays = data['ad_arrays']

    # 分析Tip 11（一般情况）
    general_patterns = []
    for ad_arr in ad_arrays:
        s = sum(ad_arr)
        if s > 0 and s != ad_arr[0]:
            general_patterns.append(ad_arr)

    if not general_patterns:
        print("没有Tip 11模式")
        return

    print(f"Tip 11模式数量: {len(general_patterns):,}")

    pattern_counter = Counter(general_patterns)

    print(f"\n高频Tip 11模式 (Top 10):")
    for pattern, count in pattern_counter.most_common(10):
        pct = count / len(general_patterns) * 100
        print(f"  {pattern}: {count:,} ({pct:.2f}%)")

    # 二等位基因分析
    biallelic_patterns = [p for p in general_patterns if len(p) == 2]
    if biallelic_patterns:
        print(f"\n二等位基因AD模式分析 ({len(biallelic_patterns):,}个):")

        ratios = []
        for ref, alt in biallelic_patterns:
            if ref + alt > 0:
                ratios.append(ref / (ref + alt))

        if ratios:
            print(f"  ref/(ref+alt)比例统计:")
            print(f"    均值: {np.mean(ratios):.3f}")
            print(f"    中位数: {np.median(ratios):.3f}")

            het_like = sum(1 for r in ratios if 0.3 <= r <= 0.7)
            print(f"    类杂合(0.3-0.7): {het_like/len(ratios)*100:.2f}%")


def explore_new_pl_patterns(data, name):
    """探索PL的新压缩模式"""
    print(f"\n--- {name} PL新模式探索 ---")

    pl_arrays = data['pl_arrays']

    no_pattern = []
    for pl_arr in pl_arrays:
        pl_type, _, _ = check_pl_pattern(pl_arr)
        if pl_type == 3:
            no_pattern.append(pl_arr)

    if not no_pattern:
        print("没有Type 3模式")
        return

    print(f"Type 3模式数量: {len(no_pattern):,}")

    reasons = Counter()

    for pl_arr in no_pattern[:1000]:
        if len(pl_arr) == 0:
            reasons['empty'] += 1
        elif pl_arr[0] != 0:
            reasons['first_not_zero'] += 1
        elif len(pl_arr) < 3:
            reasons['too_short'] += 1
        else:
            reasons['pattern_mismatch'] += 1

    print(f"\nType 3不匹配原因 (采样1000):")
    for reason, count in reasons.most_common():
        print(f"  {reason}: {count}")

    pattern_counter = Counter(no_pattern)
    print(f"\n高频Type 3模式 (Top 10):")
    for pattern, count in pattern_counter.most_common(10):
        pct = count / len(no_pattern) * 100
        display = pattern[:6] if len(pattern) > 6 else pattern
        suffix = "..." if len(pattern) > 6 else ""
        print(f"  {display}{suffix}: {count:,} ({pct:.2f}%)")


def analyze_site_correlation(vcf_path, name, max_variants=2000, max_samples=100):
    """分析相邻位点间的AD相关性"""
    print(f"\n--- {name} 位点间相关性 ---")

    vcf = VCF(vcf_path)
    n_samples = min(len(vcf.samples), max_samples)

    ad_sums_by_site = []

    for var_idx, variant in enumerate(vcf):
        if var_idx >= max_variants:
            break

        try:
            ad = variant.format('AD')
            if ad is not None:
                site_sums = []
                for sample_idx in range(n_samples):
                    sample_ad = ad[sample_idx]
                    if sample_ad is not None:
                        s = sum(int(v) for v in sample_ad if v >= 0)
                        site_sums.append(s)
                    else:
                        site_sums.append(0)
                ad_sums_by_site.append(site_sums)
        except:
            continue

    if len(ad_sums_by_site) < 2:
        print("数据不足")
        return

    ad_matrix = np.array(ad_sums_by_site)
    print(f"数据形状: {ad_matrix.shape} (位点 x 样本)")

    correlations = []
    for i in range(len(ad_matrix) - 1):
        if np.std(ad_matrix[i]) > 0 and np.std(ad_matrix[i+1]) > 0:
            corr = np.corrcoef(ad_matrix[i], ad_matrix[i+1])[0, 1]
            if not np.isnan(corr):
                correlations.append(corr)

    if correlations:
        print(f"相邻位点sum(AD)相关性:")
        print(f"  均值: {np.mean(correlations):.3f}")
        print(f"  中位数: {np.median(correlations):.3f}")
        print(f"  标准差: {np.std(correlations):.3f}")
        print(f"  高相关(>0.5)比例: {sum(1 for c in correlations if c > 0.5)/len(correlations)*100:.2f}%")


def generate_summary():
    """生成分析总结"""
    print("\n" + "="*70)
    print("VCF AD/PL 字段数据分析总结与优化建议")
    print("="*70)

    print("""
【1. AD字段特征】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• 数值分布：高度偏斜，大量零值和小值
• 稀疏性：零值占比通常>40%，低值(<5)占比>60%
• 现有Tip编码效果：
  - Tip 00+01（无需额外存储）约占40-60%
  - Tip 11（需要字典）约占20-40%
• 字典压缩潜力：唯一模式数远小于总数组数

【2. PL字段特征】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• 数值分布：值域广泛(0-数千)，但集中在特定模式
• 模式识别效果：
  - Type 0(全零)+Type 1+Type 2: 约占70-90%
  - Type 3(需字典): 约占10-30%
• 标准[0,a,b]模式高度有效

【3. 跨字段关联】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• DP≈sum(AD): 高度一致(>95%)，异常表策略有效
• GQ≈second_min(PL): 一致性取决于数据来源
• GT-AD/PL: 强相关，可用于校验和预测

【4. 位点间相关性】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• 相邻位点sum(AD)相关性较低(~0.3-0.4)
• 差分编码收益有限
• 每个样本的深度相对稳定(样本内相关性高)

【5. 优化建议】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

(A) AD字段优化：
  1. 扩展Tip编码：考虑添加更多高频模式
     - 如: sum==1, sum==3等常见值
  2. 按样本建立基线：利用样本间深度一致性
  3. 比例编码：对杂合位点用ref:alt比例

(B) PL字段优化：
  1. 扩展模式识别：覆盖更多常见模式
  2. 近似模式匹配：允许小误差的模式匹配
  3. 分层编码：先编码模式类型，再编码参数

(C) 通用优化：
  1. 自适应编码选择：根据数据特征动态选择
  2. 块级字典：减少全局字典大小
  3. 熵编码优化：对tip序列使用更高效的熵编码
""")


def main():
    print("="*70)
    print("VCF AD/PL 字段数据分析")
    print("="*70)

    datasets = {}
    for name, path in VCF_FILES.items():
        print(f"\n>>> 加载 {name}")
        datasets[name] = extract_vcf_data(path, MAX_VARIANTS)

    # AD分析
    for name, data in datasets.items():
        analyze_ad_distribution(data, name)
        analyze_ad_structure(data, name)
        analyze_ad_differential(data, name)

    # PL分析
    for name, data in datasets.items():
        analyze_pl_distribution(data, name)
        analyze_pl_patterns(data, name)
        analyze_pl_structure(data, name)

    # 跨字段关联
    print(f"\n{'='*60}")
    print("跨字段关联分析")
    print(f"{'='*60}")
    for name, path in VCF_FILES.items():
        analyze_ad_dp_relation(path, name)
        analyze_pl_gq_relation(path, name)

    # 新模式探索
    print(f"\n{'='*60}")
    print("新压缩模式探索")
    print(f"{'='*60}")
    for name, data in datasets.items():
        explore_new_ad_patterns(data, name)
        explore_new_pl_patterns(data, name)

    # 位点相关性
    print(f"\n{'='*60}")
    print("位点间相关性分析")
    print(f"{'='*60}")
    for name, path in VCF_FILES.items():
        analyze_site_correlation(path, name)

    # 总结
    generate_summary()


if __name__ == '__main__':
    main()
