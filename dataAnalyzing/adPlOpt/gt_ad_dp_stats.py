#!/usr/bin/env python3
"""
GT/AD/DP consistency stats (biallelic diploid focus).

Motivation:
  - AD is usually correlated with GT:
      GT=0/0 => (n,0) often
      GT=1/1 => (0,n) often
      GT=0/1 => both >0 often
  - DP often equals sum(AD) in many callsets.

If these correlations are strong, AD can be encoded with less payload by
conditioning on GT and/or DP (lossless with sparse "exceptions").

Usage:
  python3 dataAnalyzing/adPlOpt/gt_ad_dp_stats.py toy/final_subset.vcf.gz
  python3 dataAnalyzing/adPlOpt/gt_ad_dp_stats.py /home/binq/data/1000GPi/lc_bams.first50000.vcf.gz --max-variants 20000 --max-samples 200
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

from cyvcf2 import VCF


def _is_special_int(v: int) -> bool:
    return v < 0


@dataclass
class Stats:
    variants: int = 0
    samples_scanned: int = 0
    cells_total: int = 0
    cells_skipped: int = 0
    cells_used: int = 0

    gt_class: Counter = field(default_factory=Counter)

    # Per GT class: (ref_only, alt_only, all_zero, other)
    pat_by_gt: dict[int, Counter] = field(default_factory=lambda: {0: Counter(), 1: Counter(), 2: Counter()})

    dp_present: int = 0
    dp_eq_sumad: int = 0
    dp_ne_sumad: int = 0
    dp_delta_bucket: Counter = field(default_factory=Counter)


def scan(path: str, max_variants: Optional[int], max_samples: Optional[int]) -> tuple[Stats, int]:
    vcf = VCF(path)
    n_samples_total = len(vcf.samples)
    n_samples = n_samples_total if max_samples is None else min(n_samples_total, max_samples)

    st = Stats()
    st.samples_scanned = n_samples

    for rec in vcf:
        st.variants += 1
        if max_variants is not None and st.variants > max_variants:
            break

        st.cells_total += n_samples

        try:
            ad = rec.format("AD")
        except Exception:
            ad = None
        try:
            dp = rec.format("DP")
        except Exception:
            dp = None

        gts = getattr(rec, "genotypes", None)
        if gts is None or ad is None:
            st.cells_skipped += n_samples
            continue

        dp_ok = dp is not None

        for s in range(n_samples):
            ad_row = ad[s]
            dp_row = dp[s] if dp_ok else None
            if ad_row is None or len(ad_row) != 2:
                st.cells_skipped += 1
                continue
            if s >= len(gts) or gts[s] is None or len(gts[s]) < 2:
                st.cells_skipped += 1
                continue

            a0 = int(gts[s][0])
            a1 = int(gts[s][1])
            if a0 < 0 or a1 < 0 or a0 > 1 or a1 > 1:
                st.cells_skipped += 1
                continue
            g = a0 + a1

            r = int(ad_row[0])
            a = int(ad_row[1])
            if _is_special_int(r) or _is_special_int(a):
                st.cells_skipped += 1
                continue

            st.cells_used += 1
            st.gt_class[g] += 1

            if r == 0 and a == 0:
                st.pat_by_gt[g]["(0,0)"] += 1
            elif r > 0 and a == 0:
                st.pat_by_gt[g]["(n,0)"] += 1
            elif r == 0 and a > 0:
                st.pat_by_gt[g]["(0,n)"] += 1
            else:
                st.pat_by_gt[g]["other"] += 1

            if dp_row is not None and len(dp_row) >= 1:
                d = int(dp_row[0])
                if not _is_special_int(d):
                    st.dp_present += 1
                    ssum = r + a
                    if d == ssum:
                        st.dp_eq_sumad += 1
                    else:
                        st.dp_ne_sumad += 1
                        delta = d - ssum
                        if delta == 0:
                            st.dp_delta_bucket["0"] += 1
                        elif -2 <= delta <= 2:
                            st.dp_delta_bucket["[-2,2]"] += 1
                        elif -10 <= delta <= 10:
                            st.dp_delta_bucket["[-10,10]"] += 1
                        else:
                            st.dp_delta_bucket["other"] += 1

    return st, n_samples_total


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", help="Input VCF/VCF.GZ/BCF (must contain GT+AD; DP optional)")
    ap.add_argument("--max-variants", type=int, default=None)
    ap.add_argument("--max-samples", type=int, default=None)
    args = ap.parse_args()

    st, n_samples_total = scan(args.vcf, args.max_variants, args.max_samples)
    used = st.cells_used or 1

    print("=" * 72)
    print(f"Input: {args.vcf}")
    print(f"Scanned: variants={st.variants:,} samples={st.samples_scanned:,} (total samples={n_samples_total:,})")
    print("")
    print("[AD len==2 + GT biallelic diploid]")
    print(f"cells_total:   {st.cells_total:,}")
    print(f"cells_skipped: {st.cells_skipped:,}")
    print(f"cells_used:    {st.cells_used:,}")
    print("")
    print(f"GT class (0/0,0/1,1/1): {st.gt_class.get(0,0):,}/{st.gt_class.get(1,0):,}/{st.gt_class.get(2,0):,}")
    for g in (0, 1, 2):
        total_g = st.gt_class.get(g, 0)
        if total_g <= 0:
            continue
        c = st.pat_by_gt[g]
        print(f"GT={g}: (n,0)={c.get('(n,0)',0):,} (0,n)={c.get('(0,n)',0):,} (0,0)={c.get('(0,0)',0):,} other={c.get('other',0):,}")

    if st.dp_present:
        print("")
        print("[DP vs sum(AD)]")
        print(f"DP present cells: {st.dp_present:,} ({st.dp_present / used:.2%} of used)")
        print(f"DP == sum(AD):    {st.dp_eq_sumad:,} ({st.dp_eq_sumad / st.dp_present:.2%})")
        print(f"DP != sum(AD):    {st.dp_ne_sumad:,} ({st.dp_ne_sumad / st.dp_present:.2%})")
        if st.dp_delta_bucket:
            print(f"DP-sum(AD) delta buckets (only mismatches): {dict(st.dp_delta_bucket)}")


if __name__ == "__main__":
    main()
