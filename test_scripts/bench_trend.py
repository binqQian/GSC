#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import os
import sys
from collections import defaultdict


SPARK_CHARS = "▁▂▃▄▅▆▇█"


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize bench_history trends into markdown + long CSV."
    )
    p.add_argument("--in", dest="input_csv", default="docs/bench/bench_history.csv")
    p.add_argument("--out", dest="out_dir", default="")
    p.add_argument("--last", dest="last_n", type=int, default=10)
    p.add_argument(
        "--filter",
        action="append",
        default=[],
        help="Filter rows by key=value (repeatable).",
    )
    p.add_argument(
        "--group-by",
        default="dataset,mode,step,host,cpu,threads,compressor",
        help="Comma-separated group-by fields.",
    )
    p.add_argument(
        "--metrics",
        default="elapsed_s,rss_kb,throughput_vps,metric_a",
        help="Comma-separated metric fields.",
    )
    return p.parse_args()


def parse_filters(filter_args):
    filters = {}
    for item in filter_args:
        if "=" not in item:
            raise ValueError(f"Invalid filter: {item} (expected key=value)")
        k, v = item.split("=", 1)
        filters[k.strip()] = v.strip()
    return filters


def sparkline(values):
    if not values:
        return ""
    vmin = min(values)
    vmax = max(values)
    if vmin == vmax:
        return SPARK_CHARS[len(SPARK_CHARS) // 2] * len(values)
    span = vmax - vmin
    chars = []
    for v in values:
        idx = int((v - vmin) / span * (len(SPARK_CHARS) - 1))
        chars.append(SPARK_CHARS[idx])
    return "".join(chars)


def fmt_num(x):
    if x is None:
        return "NA"
    if isinstance(x, int):
        return str(x)
    return f"{x:.6g}"


def safe_md(val):
    return str(val).replace("|", "/")


def main():
    args = parse_args()
    try:
        filters = parse_filters(args.filter)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    input_csv = args.input_csv
    if not os.path.isfile(input_csv):
        print(f"ERROR: input not found: {input_csv}", file=sys.stderr)
        return 1

    group_by = [f.strip() for f in args.group_by.split(",") if f.strip()]
    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]

    if not group_by:
        print("ERROR: --group-by is empty", file=sys.stderr)
        return 1
    if not metrics:
        print("ERROR: --metrics is empty", file=sys.stderr)
        return 1

    out_dir = args.out_dir
    if not out_dir:
        stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = f"tmp/bench_trend_{stamp}"
    os.makedirs(out_dir, exist_ok=True)

    with open(input_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print("ERROR: input has no rows", file=sys.stderr)
        return 1

    for key in filters:
        if key not in rows[0]:
            print(f"ERROR: filter key not in header: {key}", file=sys.stderr)
            return 1

    filtered = []
    for i, row in enumerate(rows):
        ok = True
        for k, v in filters.items():
            if row.get(k, "") != v:
                ok = False
                break
        if ok:
            row["_row_index"] = i
            filtered.append(row)

    if not filtered:
        print("ERROR: no rows matched filters", file=sys.stderr)
        return 1

    def row_key(r):
        return tuple(r.get(k, "") for k in group_by)

    groups = defaultdict(list)
    for r in filtered:
        groups[row_key(r)].append(r)

    def row_sort_key(r):
        return (r.get("date", ""), r.get("_row_index", 0))

    for k in groups:
        groups[k].sort(key=row_sort_key)

    long_csv_path = os.path.join(out_dir, "trend_long.csv")
    with open(long_csv_path, "w", newline="") as f:
        fieldnames = ["date"] + group_by + ["metric", "value", "git", "notes"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for key, items in groups.items():
            for r in items:
                for m in metrics:
                    v = r.get(m, "")
                    try:
                        v_num = float(v)
                    except ValueError:
                        continue
                    out = {k: r.get(k, "") for k in group_by}
                    out.update(
                        {
                            "date": r.get("date", ""),
                            "metric": m,
                            "value": f"{v_num:.6g}",
                            "git": r.get("git", ""),
                            "notes": r.get("notes", ""),
                        }
                    )
                    writer.writerow(out)

    summary_path = os.path.join(out_dir, "trend_summary.md")
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(summary_path, "w") as f:
        f.write("# Bench Trend Summary\n\n")
        f.write(f"- input: `{input_csv}`\n")
        f.write(f"- output: `{out_dir}`\n")
        f.write(f"- last_n: `{args.last_n}`\n")
        if filters:
            f.write(f"- filters: `{filters}`\n")
        f.write(f"- generated: `{now}`\n\n")

        for metric in metrics:
            f.write(f"## {metric}\n\n")
            headers = group_by + ["last", "prev", "delta", "delta_pct", "sparkline"]
            f.write("|" + "|".join(headers) + "|\n")
            f.write("|" + "|".join(["---"] * len(headers)) + "|\n")

            for key, items in sorted(groups.items()):
                values = []
                for r in items:
                    v = r.get(metric, "")
                    try:
                        values.append(float(v))
                    except ValueError:
                        continue
                if not values:
                    continue
                tail = values[-args.last_n :]
                last = tail[-1]
                prev = tail[-2] if len(tail) >= 2 else None
                delta = last - prev if prev is not None else None
                delta_pct = (delta / prev * 100) if prev not in (None, 0) else None

                line = [safe_md(v) for v in key]
                line += [
                    fmt_num(last),
                    fmt_num(prev),
                    fmt_num(delta),
                    fmt_num(delta_pct),
                    sparkline(tail),
                ]
                f.write("|" + "|".join(line) + "|\n")
            f.write("\n")

    print(f"Wrote {summary_path}")
    print(f"Wrote {long_csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
