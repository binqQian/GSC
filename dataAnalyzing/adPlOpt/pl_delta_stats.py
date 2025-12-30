#!/usr/bin/env python3
"""
Estimate whether delta-coding PL values across nearby variants is promising.

For biallelic diploid sites (PL length == 3) and non-special values:
  - Focus on positions 1 (het) and 2 (hom-alt) in the VCF-defined order.
  - Optionally restrict to records where PL[0]==0 (common normalized case).

We compute per-sample deltas along variant order:
  delta = curr - prev
and report the distribution of |delta| and delta==0 rates.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from typing import Optional

from cyvcf2 import VCF


def _is_special_int(v: int) -> bool:
    return v < 0


@dataclass
class DeltaStats:
    n: int = 0
    zero: int = 0
    abs_lt_2: int = 0
    abs_lt_4: int = 0
    abs_lt_8: int = 0
    abs_lt_16: int = 0
    abs_lt_32: int = 0
    abs_lt_64: int = 0
    abs_lt_128: int = 0
    abs_lt_256: int = 0
    abs_lt_512: int = 0
    abs_lt_1024: int = 0
    abs_lt_2048: int = 0
    abs_lt_4096: int = 0
    abs_ge_4096: int = 0
    sign: Counter = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        self.sign = Counter()

    def add(self, delta: int) -> None:
        self.n += 1
        if delta == 0:
            self.zero += 1
        self.sign["neg" if delta < 0 else "pos" if delta > 0 else "zero"] += 1
        a = abs(delta)
        if a < 2:
            self.abs_lt_2 += 1
        if a < 4:
            self.abs_lt_4 += 1
        if a < 8:
            self.abs_lt_8 += 1
        if a < 16:
            self.abs_lt_16 += 1
        if a < 32:
            self.abs_lt_32 += 1
        if a < 64:
            self.abs_lt_64 += 1
        if a < 128:
            self.abs_lt_128 += 1
        if a < 256:
            self.abs_lt_256 += 1
        if a < 512:
            self.abs_lt_512 += 1
        if a < 1024:
            self.abs_lt_1024 += 1
        if a < 2048:
            self.abs_lt_2048 += 1
        if a < 4096:
            self.abs_lt_4096 += 1
        else:
            self.abs_ge_4096 += 1

    def _pct(self, x: int) -> str:
        return "n/a" if self.n == 0 else f"{x / self.n:.2%}"

    def report(self, title: str) -> None:
        print(f"\n[{title}] deltas={self.n:,}")
        if self.n == 0:
            return
        print(f"delta==0: {self.zero:,} ({self._pct(self.zero)})")
        print(f"|d|<2: {self.abs_lt_2:,} ({self._pct(self.abs_lt_2)})")
        print(f"|d|<8: {self.abs_lt_8:,} ({self._pct(self.abs_lt_8)})")
        print(f"|d|<32: {self.abs_lt_32:,} ({self._pct(self.abs_lt_32)})")
        print(f"|d|<128: {self.abs_lt_128:,} ({self._pct(self.abs_lt_128)})")
        print(f"|d|<512: {self.abs_lt_512:,} ({self._pct(self.abs_lt_512)})")
        print(f"|d|<2048: {self.abs_lt_2048:,} ({self._pct(self.abs_lt_2048)})")
        print(f"|d|>=4096: {self.abs_ge_4096:,} ({self._pct(self.abs_ge_4096)})")
        print(f"sign: {dict(self.sign)}")


def scan(path: str, max_variants: Optional[int], max_samples: int, require_pl0_zero: bool) -> tuple[DeltaStats, DeltaStats, int, int]:
    vcf = VCF(path)
    n_samples_total = len(vcf.samples)
    n_samples = min(n_samples_total, max_samples)

    prev1 = [0] * n_samples
    prev2 = [0] * n_samples
    has_prev = [False] * n_samples

    d1 = DeltaStats()
    d2 = DeltaStats()

    n_variants = 0
    used_variants = 0
    for rec in vcf:
        n_variants += 1
        if max_variants is not None and n_variants > max_variants:
            break

        try:
            pl = rec.format("PL")
        except Exception:
            pl = None
        if pl is None:
            continue

        used_variants += 1
        for s in range(n_samples):
            row = pl[s]
            if row is None or len(row) != 3:
                continue
            v0, v1, v2 = int(row[0]), int(row[1]), int(row[2])
            if _is_special_int(v0) or _is_special_int(v1) or _is_special_int(v2):
                has_prev[s] = False
                continue
            if require_pl0_zero and v0 != 0:
                has_prev[s] = False
                continue

            if has_prev[s]:
                d1.add(v1 - prev1[s])
                d2.add(v2 - prev2[s])
            prev1[s] = v1
            prev2[s] = v2
            has_prev[s] = True

    return d1, d2, n_variants, n_samples


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", nargs="+", help="VCF/VCF.GZ path(s)")
    ap.add_argument("--max-variants", type=int, default=20000)
    ap.add_argument("--max-samples", type=int, default=200)
    ap.add_argument("--require-pl0-zero", action="store_true", help="only consider rows with PL[0]==0; resets delta chain otherwise")
    args = ap.parse_args()

    for path in args.vcf:
        d1, d2, n_variants, n_samples = scan(path, args.max_variants, args.max_samples, args.require_pl0_zero)
        print("=" * 72)
        print(f"Input: {path}")
        print(f"Scanned: variants={n_variants:,} samples={n_samples:,} require_pl0_zero={args.require_pl0_zero}")
        d1.report("PL[1] (het)")
        d2.report("PL[2] (hom-alt)")


if __name__ == "__main__":
    main()

