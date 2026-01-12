#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import re
import sys


_EXP_NORMALIZE_RE = re.compile(r"0+e\\+")


def _canonicalize_info(info: str) -> str:
    if info == "." or info == "":
        return "."
    parts = info.split(";")
    parts = [p for p in parts if p]
    parts.sort()
    return ";".join(parts) if parts else "."


def _canonicalize_format_and_samples(fmt: str, samples: list[str]) -> tuple[str, list[str]]:
    if fmt == "." or fmt == "":
        return ".", samples
    keys = fmt.split(":")
    if len(keys) <= 1:
        return fmt, samples

    # Canonical key order: lexical, stable.
    order = sorted(range(len(keys)), key=lambda i: keys[i])
    canon_keys = [keys[i] for i in order]

    canon_samples: list[str] = []
    for s in samples:
        vals = s.split(":")
        if len(vals) < len(keys):
            vals = vals + ["."] * (len(keys) - len(vals))
        canon_vals = [vals[i] for i in order]
        canon_samples.append(":".join(canon_vals))

    return ":".join(canon_keys), canon_samples


def main() -> int:
    h = hashlib.sha256()

    for line in sys.stdin:
        if not line or line[0] == "#":
            continue
        line = line.rstrip("\n")
        cols = line.split("\t")
        if len(cols) < 8:
            continue

        cols[7] = _canonicalize_info(cols[7])
        if len(cols) >= 9:
            cols[8], cols[9:] = _canonicalize_format_and_samples(cols[8], cols[9:])

        canon = "\t".join(cols)
        canon = _EXP_NORMALIZE_RE.sub("e+", canon)

        h.update(canon.encode("utf-8"))
        h.update(b"\n")

    sys.stdout.write(h.hexdigest() + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

