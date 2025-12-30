#!/usr/bin/env python3
"""
Quick, streaming stats for AD/PL to guide codec decisions.

Focus:
  - PL (per_sample==3): argmin index distribution and "min==0" rate
  - AD (per_sample==2): (n,0) / (0,n) / other distribution

Usage examples:
  python3 dataAnalyzing/adPlOpt/quick_stats.py toy/final_subset.vcf.gz
  python3 dataAnalyzing/adPlOpt/quick_stats.py /home/binq/data/1000GPi/lc_bams.first50000.vcf.gz --max-variants 2000
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from typing import Optional

from cyvcf2 import VCF


@dataclass
class Pl3Stats:
    total: int = 0
    has_special: int = 0
    min_is_zero: int = 0
    argmin: Counter = None  # type: ignore[assignment]
    pl0_is_zero: int = 0

    def __post_init__(self) -> None:
        self.argmin = Counter()


@dataclass
class Ad2Stats:
    total: int = 0
    has_special: int = 0
    all_zero: int = 0
    ref_only: int = 0  # (n,0), n>0
    alt_only: int = 0  # (0,n), n>0
    other: int = 0  # everything else (including both-nonzero)


def _is_special_int(v: int) -> bool:
    # bcf_int32_missing or bcf_int32_vector_end or negative values
    return v < 0


def scan_vcf(path: str, max_variants: Optional[int], max_samples: Optional[int]) -> tuple[Pl3Stats, Ad2Stats, int, int]:
    vcf = VCF(path)
    n_samples_total = len(vcf.samples)
    n_samples = n_samples_total if max_samples is None else min(n_samples_total, max_samples)

    pl3 = Pl3Stats()
    ad2 = Ad2Stats()

    n_variants = 0
    for rec in vcf:
        n_variants += 1
        if max_variants is not None and n_variants > max_variants:
            break

        # PL
        try:
            pl = rec.format("PL")
        except Exception:
            pl = None
        if pl is not None:
            for s in range(n_samples):
                row = pl[s]
                if row is None or len(row) != 3:
                    continue
                pl3.total += 1
                a0, a1, a2 = int(row[0]), int(row[1]), int(row[2])
                if _is_special_int(a0) or _is_special_int(a1) or _is_special_int(a2):
                    pl3.has_special += 1
                    continue
                if a0 == 0:
                    pl3.pl0_is_zero += 1
                # argmin and min==0
                vals = (a0, a1, a2)
                m = 0
                if a1 < vals[m]:
                    m = 1
                if a2 < vals[m]:
                    m = 2
                pl3.argmin[m] += 1
                if vals[m] == 0:
                    pl3.min_is_zero += 1

        # AD
        try:
            ad = rec.format("AD")
        except Exception:
            ad = None
        if ad is not None:
            for s in range(n_samples):
                row = ad[s]
                if row is None or len(row) != 2:
                    continue
                ad2.total += 1
                v0, v1 = int(row[0]), int(row[1])
                if _is_special_int(v0) or _is_special_int(v1):
                    ad2.has_special += 1
                    continue
                if v0 == 0 and v1 == 0:
                    ad2.all_zero += 1
                elif v0 > 0 and v1 == 0:
                    ad2.ref_only += 1
                elif v0 == 0 and v1 > 0:
                    ad2.alt_only += 1
                else:
                    ad2.other += 1

    return pl3, ad2, n_variants, n_samples


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", nargs="+", help="VCF/VCF.GZ path(s)")
    ap.add_argument("--max-variants", type=int, default=None, help="limit records scanned")
    ap.add_argument("--max-samples", type=int, default=200, help="limit samples scanned (default: 200)")
    args = ap.parse_args()

    for path in args.vcf:
        pl3, ad2, n_variants, n_samples = scan_vcf(path, args.max_variants, args.max_samples)

        print("=" * 72)
        print(f"Input: {path}")
        print(f"Scanned: variants={n_variants:,} samples={n_samples:,}")

        if pl3.total:
            print("\n[PL len==3]")
            print(f"cells: {pl3.total:,} (special: {pl3.has_special:,}, {pl3.has_special / pl3.total:.2%})")
            print(f"PL[0]==0: {pl3.pl0_is_zero:,} ({pl3.pl0_is_zero / pl3.total:.2%})")
            print(f"min==0: {pl3.min_is_zero:,} ({pl3.min_is_zero / pl3.total:.2%})")
            denom = max(1, pl3.total - pl3.has_special)
            argmin0 = pl3.argmin.get(0, 0)
            argmin1 = pl3.argmin.get(1, 0)
            argmin2 = pl3.argmin.get(2, 0)
            print(f"argmin idx (0/1/2): {argmin0:,}/{argmin1:,}/{argmin2:,} (over non-special ~{denom:,})")
        else:
            print("\n[PL len==3] no cells")

        if ad2.total:
            print("\n[AD len==2]")
            print(f"cells: {ad2.total:,} (special: {ad2.has_special:,}, {ad2.has_special / ad2.total:.2%})")
            non_special = max(1, ad2.total - ad2.has_special)
            print(f"(0,0): {ad2.all_zero:,} ({ad2.all_zero / non_special:.2%} of non-special)")
            print(f"(n,0): {ad2.ref_only:,} ({ad2.ref_only / non_special:.2%} of non-special)")
            print(f"(0,n): {ad2.alt_only:,} ({ad2.alt_only / non_special:.2%} of non-special)")
            print(f"other: {ad2.other:,} ({ad2.other / non_special:.2%} of non-special)")
        else:
            print("\n[AD len==2] no cells")


if __name__ == "__main__":
    main()

