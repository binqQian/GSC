#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import os
import sys
from collections import defaultdict


DEFAULT_THRESHOLDS = {
    "elapsed_s": 10.0,       # % increase
    "rss_kb": 15.0,          # % increase
    "throughput_vps": 10.0,  # % decrease
    "metric_a": 5.0,         # % increase (compressed bytes)
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Detect regressions from bench_history.csv using last two samples."
    )
    p.add_argument("--in", dest="input_csv", default="docs/bench/bench_history.csv")
    p.add_argument("--out", dest="out_dir", default="")
    p.add_argument("--last", dest="last_n", type=int, default=2)
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
    p.add_argument("--elapsed-pct", type=float, default=DEFAULT_THRESHOLDS["elapsed_s"])
    p.add_argument("--rss-pct", type=float, default=DEFAULT_THRESHOLDS["rss_kb"])
    p.add_argument("--throughput-pct", type=float, default=DEFAULT_THRESHOLDS["throughput_vps"])
    p.add_argument("--metric-a-pct", type=float, default=DEFAULT_THRESHOLDS["metric_a"])
    return p.parse_args()


def parse_filters(filter_args):
    filters = {}
    for item in filter_args:
        if "=" not in item:
            raise ValueError(f"Invalid filter: {item} (expected key=value)")
        k, v = item.split("=", 1)
        filters[k.strip()] = v.strip()
    return filters


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
    if args.last_n < 2:
        print("ERROR: --last must be >= 2", file=sys.stderr)
        return 1

    out_dir = args.out_dir
    if not out_dir:
        stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = f"tmp/bench_alert_{stamp}"
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

    thresholds = {
        "elapsed_s": args.elapsed_pct,
        "rss_kb": args.rss_pct,
        "throughput_vps": args.throughput_pct,
        "metric_a": args.metric_a_pct,
    }

    def is_regression(metric, delta_pct):
        if delta_pct is None:
            return False
        thr = thresholds.get(metric)
        if thr is None:
            return False
        if metric in ("elapsed_s", "rss_kb", "metric_a"):
            return delta_pct >= thr
        if metric == "throughput_vps":
            return delta_pct <= -thr
        return False

    csv_path = os.path.join(out_dir, "bench_alert.csv")
    with open(csv_path, "w", newline="") as f:
        fieldnames = group_by + [
            "metric",
            "prev_date",
            "last_date",
            "prev",
            "last",
            "delta",
            "delta_pct",
            "threshold_pct",
            "status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        alert_rows = []
        for key, items in groups.items():
            for metric in metrics:
                values = []
                dates = []
                for r in items:
                    v = r.get(metric, "")
                    try:
                        v_num = float(v)
                    except ValueError:
                        continue
                    values.append(v_num)
                    dates.append(r.get("date", ""))
                if len(values) < 2:
                    continue
                tail_vals = values[-args.last_n :]
                tail_dates = dates[-args.last_n :]
                prev = tail_vals[-2]
                last = tail_vals[-1]
                prev_date = tail_dates[-2]
                last_date = tail_dates[-1]
                delta = last - prev
                delta_pct = None if prev == 0 else (delta / prev * 100)
                status = "REGRESSION" if is_regression(metric, delta_pct) else "OK"
                row = {k: v for k, v in zip(group_by, key)}
                row.update(
                    {
                        "metric": metric,
                        "prev_date": prev_date,
                        "last_date": last_date,
                        "prev": f"{prev:.6g}",
                        "last": f"{last:.6g}",
                        "delta": f"{delta:.6g}",
                        "delta_pct": "NA" if delta_pct is None else f"{delta_pct:.6g}",
                        "threshold_pct": thresholds.get(metric, "NA"),
                        "status": status,
                    }
                )
                writer.writerow(row)
                if status == "REGRESSION":
                    alert_rows.append(row)

    md_path = os.path.join(out_dir, "bench_alert.md")
    with open(md_path, "w") as f:
        f.write("# Bench Regression Report\n\n")
        f.write(f"- input: `{input_csv}`\n")
        f.write(f"- output: `{out_dir}`\n")
        f.write(f"- last_n: `{args.last_n}`\n")
        if filters:
            f.write(f"- filters: `{filters}`\n")
        f.write(
            "- thresholds(%): "
            f"elapsed_s>={thresholds['elapsed_s']}, "
            f"rss_kb>={thresholds['rss_kb']}, "
            f"throughput_vps<=-{thresholds['throughput_vps']}, "
            f"metric_a>={thresholds['metric_a']}\n\n"
        )

        if not alert_rows:
            f.write("No regressions detected.\n")
        else:
            headers = group_by + [
                "metric",
                "prev_date",
                "last_date",
                "prev",
                "last",
                "delta",
                "delta_pct",
                "threshold_pct",
            ]
            f.write("|" + "|".join(headers) + "|\n")
            f.write("|" + "|".join(["---"] * len(headers)) + "|\n")
            for row in alert_rows:
                line = [safe_md(row[h]) for h in headers]
                f.write("|" + "|".join(line) + "|\n")

    print(f"Wrote {md_path}")
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
