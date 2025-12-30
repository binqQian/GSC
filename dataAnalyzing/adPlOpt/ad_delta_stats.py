#!/usr/bin/env python3
"""
Estimate whether delta-coding AD across nearby variants is promising.

For biallelic sites (AD length == 2) and non-special values:
  - We focus on the common "ref-only" pattern (AD[1] == 0), where the codec stores AD[0].
  - We compute per-sample deltas of AD[0] along variant order and report |delta| distribution.
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
    abs_ge_256: int = 0
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
        else:
            self.abs_ge_256 += 1

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
        print(f"|d|<256: {self.abs_lt_256:,} ({self._pct(self.abs_lt_256)})")
        print(f"|d|>=256: {self.abs_ge_256:,} ({self._pct(self.abs_ge_256)})")
        print(f"sign: {dict(self.sign)}")


def scan(path: str, max_variants: Optional[int], max_samples: int) -> tuple[DeltaStats, int, int, int]:
    vcf = VCF(path)
    n_samples_total = len(vcf.samples)
    n_samples = min(n_samples_total, max_samples)

    prev = [0] * n_samples
    has_prev = [False] * n_samples
    d = DeltaStats()

    n_variants = 0
    cells = 0
    for rec in vcf:
        n_variants += 1
        if max_variants is not None and n_variants > max_variants:
            break

        try:
            ad = rec.format("AD")
        except Exception:
            ad = None
        if ad is None:
            continue

        for s in range(n_samples):
            row = ad[s]
            if row is None or len(row) != 2:
                continue
            v0, v1 = int(row[0]), int(row[1])
            if _is_special_int(v0) or _is_special_int(v1):
                has_prev[s] = False
                continue
            if v1 != 0:
                has_prev[s] = False
                continue

            cells += 1
            if has_prev[s]:
                d.add(v0 - prev[s])
            prev[s] = v0
            has_prev[s] = True

    return d, n_variants, n_samples, cells


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", nargs="+", help="VCF/VCF.GZ path(s)")
    ap.add_argument("--max-variants", type=int, default=20000)
    ap.add_argument("--max-samples", type=int, default=200)
    args = ap.parse_args()

    for path in args.vcf:
        d, n_variants, n_samples, cells = scan(path, args.max_variants, args.max_samples)
        print("=" * 72)
        print(f"Input: {path}")
        print(f"Scanned: variants={n_variants:,} samples={n_samples:,} ref-only cells={cells:,}")
        d.report("AD[0] when AD[1]==0")


if __name__ == "__main__":
    main()

