#!/usr/bin/env python3
"""
GT/PL consistency stats (biallelic diploid focus).

Motivation:
  - For PL (len==3), the minimum entry should match the called GT most of the time.
  - If mismatches are rare, we can use GT as "free" side information for PL coding
    (e.g., avoid storing minpos, or split residual streams by GT).

Usage:
  python3 dataAnalyzing/adPlOpt/gt_pl_consistency.py toy/final_subset.vcf.gz
  python3 dataAnalyzing/adPlOpt/gt_pl_consistency.py /home/binq/data/output_subset_5000_20000.vcf.gz --max-variants 20000 --max-samples 200
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

from cyvcf2 import VCF


def _is_special_int(v: int) -> bool:
    return v < 0


def _bucket_u32(v: int) -> str:
    if v < 0:
        return "<0"
    if v < 128:
        return "<128"
    if v < 256:
        return "<256"
    if v < 512:
        return "<512"
    if v < 1024:
        return "<1024"
    if v < 2048:
        return "<2048"
    if v < 4096:
        return "<4096"
    if v < 8192:
        return "<8192"
    return ">=8192"


@dataclass
class PlGtStats:
    variants: int = 0
    samples_scanned: int = 0
    cells_total: int = 0
    cells_skipped: int = 0
    cells_used: int = 0

    min_is_zero: int = 0
    gt_is_min: int = 0
    gt_val_is_zero: int = 0
    pl0_is_zero: int = 0

    gt_class: Counter = field(default_factory=Counter)
    mismatch_gt_class: Counter = field(default_factory=Counter)
    mismatch_minpos: Counter = field(default_factory=Counter)

    other_min_bucket_by_gt: dict[int, Counter] = field(default_factory=lambda: {0: Counter(), 1: Counter(), 2: Counter()})
    other_max_bucket_by_gt: dict[int, Counter] = field(default_factory=lambda: {0: Counter(), 1: Counter(), 2: Counter()})
    other_min_max_by_gt: dict[int, int] = field(default_factory=lambda: {0: 0, 1: 0, 2: 0})
    other_max_max_by_gt: dict[int, int] = field(default_factory=lambda: {0: 0, 1: 0, 2: 0})


def scan(path: str, max_variants: Optional[int], max_samples: Optional[int]) -> tuple[PlGtStats, int]:
    vcf = VCF(path)
    n_samples_total = len(vcf.samples)
    n_samples = n_samples_total if max_samples is None else min(n_samples_total, max_samples)

    st = PlGtStats()
    st.samples_scanned = n_samples

    for rec in vcf:
        st.variants += 1
        if max_variants is not None and st.variants > max_variants:
            break

        st.cells_total += n_samples

        try:
            pl = rec.format("PL")
        except Exception:
            pl = None
        gts = getattr(rec, "genotypes", None)
        if gts is None or pl is None:
            st.cells_skipped += n_samples
            continue

        for s in range(n_samples):
            pl_row = pl[s]
            if pl_row is None or len(pl_row) != 3:
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

            p0 = int(pl_row[0])
            p1 = int(pl_row[1])
            p2 = int(pl_row[2])
            if _is_special_int(p0) or _is_special_int(p1) or _is_special_int(p2):
                st.cells_skipped += 1
                continue

            st.cells_used += 1
            st.gt_class[g] += 1
            if p0 == 0:
                st.pl0_is_zero += 1

            minv = p0
            minpos = 0
            if p1 < minv:
                minv = p1
                minpos = 1
            if p2 < minv:
                minv = p2
                minpos = 2

            if minv == 0:
                st.min_is_zero += 1

            gt_val = (p0, p1, p2)[g]
            if gt_val == 0:
                st.gt_val_is_zero += 1

            if gt_val == minv:
                st.gt_is_min += 1
                other = [p0, p1, p2]
                other.pop(g)
                other_min = min(other)
                other_max = max(other)
                st.other_min_bucket_by_gt[g][_bucket_u32(other_min)] += 1
                st.other_max_bucket_by_gt[g][_bucket_u32(other_max)] += 1
                st.other_min_max_by_gt[g] = max(st.other_min_max_by_gt[g], other_min)
                st.other_max_max_by_gt[g] = max(st.other_max_max_by_gt[g], other_max)
            else:
                st.mismatch_gt_class[g] += 1
                st.mismatch_minpos[minpos] += 1

    return st, n_samples_total


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", help="Input VCF/VCF.GZ/BCF (must contain GT+PL)")
    ap.add_argument("--max-variants", type=int, default=None)
    ap.add_argument("--max-samples", type=int, default=None)
    args = ap.parse_args()

    st, n_samples_total = scan(args.vcf, args.max_variants, args.max_samples)

    used = st.cells_used or 1
    print("=" * 72)
    print(f"Input: {args.vcf}")
    print(f"Scanned: variants={st.variants:,} samples={st.samples_scanned:,} (total samples={n_samples_total:,})")
    print("")
    print("[PL len==3 + GT biallelic diploid]")
    print(f"cells_total:   {st.cells_total:,}")
    print(f"cells_skipped: {st.cells_skipped:,}")
    print(f"cells_used:    {st.cells_used:,}")
    print("")
    print(f"min==0:        {st.min_is_zero:,} ({st.min_is_zero / used:.2%})")
    print(f"PL[0]==0:      {st.pl0_is_zero:,} ({st.pl0_is_zero / used:.2%})")
    print(f"PL[GT]==0:     {st.gt_val_is_zero:,} ({st.gt_val_is_zero / used:.2%})")
    print(f"PL[GT]==min:   {st.gt_is_min:,} ({st.gt_is_min / used:.2%})")
    print("")
    print(f"GT class (0/0,0/1,1/1): {st.gt_class.get(0,0):,}/{st.gt_class.get(1,0):,}/{st.gt_class.get(2,0):,}")
    if st.mismatch_gt_class:
        mism = sum(st.mismatch_gt_class.values())
        print(f"mismatch cells: {mism:,} ({mism / used:.2%})")
        print(f"mismatch by GT:  {st.mismatch_gt_class.get(0,0):,}/{st.mismatch_gt_class.get(1,0):,}/{st.mismatch_gt_class.get(2,0):,}")
        print(f"argmin idx (0/1/2) over mismatches: {st.mismatch_minpos.get(0,0):,}/{st.mismatch_minpos.get(1,0):,}/{st.mismatch_minpos.get(2,0):,}")

    print("")
    print("[Other (non-GT) PL values buckets, only when PL[GT]==min]")
    for g in (0, 1, 2):
        total_g = st.gt_class.get(g, 0) - st.mismatch_gt_class.get(g, 0)
        if total_g <= 0:
            continue
        print(f"GT={g}: count={total_g:,} other_min_max={st.other_min_max_by_gt[g]} other_max_max={st.other_max_max_by_gt[g]}")
        print(f"  other_min buckets: {dict(st.other_min_bucket_by_gt[g])}")
        print(f"  other_max buckets: {dict(st.other_max_bucket_by_gt[g])}")


if __name__ == "__main__":
    main()
