#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  test_scripts/bench_merge.sh [options] bench1.csv [bench2.csv ...]

Options:
  -o, --out FILE           History CSV output (default: docs/bench/bench_history.csv)
  --diff FILE              Write diff CSV (requires baseline+compare)
  --baseline FILE          Baseline bench.csv for diff
  --compare FILE           Compare bench.csv for diff
  --date DATE              YYYY-MM-DD (default: today)
  --git SHA                Git short sha (default: auto)
  --host HOST              Hostname (default: auto)
  --cpu MODEL              CPU model (default: auto)
  --threads N              Threads used in bench (default: NA)
  --compressor NAME        Compressor (default: NA)
  --notes TEXT             Notes (default: empty)
  -h, --help               Show this help

Examples:
  THREADS=16 COMPRESSOR=bsc test_scripts/bench_merge.sh tmp/bench_*/bench.csv
  test_scripts/bench_merge.sh --diff tmp/bench_diff.csv --baseline A/bench.csv --compare B/bench.csv
EOF
}

sanitize_csv_field() {
  local s="$1"
  s=${s//$'\n'/ }
  s=${s//,/;}
  printf "%s" "$s"
}

OUT="docs/bench/bench_history.csv"
DIFF=""
BASELINE=""
COMPARE=""
DATE="${DATE:-}"
GIT="${GIT:-}"
HOST="${HOST:-}"
CPU="${CPU:-}"
THREADS="${THREADS:-}"
COMPRESSOR="${COMPRESSOR:-}"
NOTES="${NOTES:-}"

inputs=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--out)
      OUT="$2"
      shift 2
      ;;
    --diff)
      DIFF="$2"
      shift 2
      ;;
    --baseline)
      BASELINE="$2"
      shift 2
      ;;
    --compare)
      COMPARE="$2"
      shift 2
      ;;
    --date)
      DATE="$2"
      shift 2
      ;;
    --git)
      GIT="$2"
      shift 2
      ;;
    --host)
      HOST="$2"
      shift 2
      ;;
    --cpu)
      CPU="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --compressor)
      COMPRESSOR="$2"
      shift 2
      ;;
    --notes)
      NOTES="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      inputs+=("$1")
      shift
      ;;
  esac
done

if [[ ${#inputs[@]} -eq 0 && -z "${DIFF}" ]]; then
  usage >&2
  exit 1
fi

if [[ -z "$DATE" ]]; then
  DATE="$(date +%F)"
fi
if [[ -z "$GIT" ]]; then
  GIT="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
fi
if [[ -z "$HOST" ]]; then
  HOST="$(hostname 2>/dev/null || echo unknown)"
fi
if [[ -z "$CPU" ]]; then
  if command -v lscpu >/dev/null 2>&1; then
    CPU="$(lscpu | awk -F: '/Model name/ {sub(/^[ \t]+/,"",$2); print $2; exit}')"
  fi
  if [[ -z "$CPU" && -r /proc/cpuinfo ]]; then
    CPU="$(awk -F: '/model name/ {sub(/^[ \t]+/,"",$2); print $2; exit}' /proc/cpuinfo)"
  fi
  CPU="${CPU:-unknown}"
fi
THREADS="${THREADS:-NA}"
COMPRESSOR="${COMPRESSOR:-NA}"

DATE="$(sanitize_csv_field "$DATE")"
GIT="$(sanitize_csv_field "$GIT")"
HOST="$(sanitize_csv_field "$HOST")"
CPU="$(sanitize_csv_field "$CPU")"
THREADS="$(sanitize_csv_field "$THREADS")"
COMPRESSOR="$(sanitize_csv_field "$COMPRESSOR")"
NOTES="$(sanitize_csv_field "$NOTES")"

header="date,git,host,cpu,threads,compressor,dataset,mode,step,variants,elapsed_s,rss_kb,throughput_vps,metric_a,metric_b,notes"

if [[ -n "$OUT" && "${#inputs[@]}" -gt 0 ]]; then
  if [[ -f "$OUT" ]]; then
    read -r first_line < "$OUT" || true
    if [[ "$first_line" != "$header" ]]; then
      echo "ERROR: ${OUT} header mismatch; expected:" >&2
      echo "  $header" >&2
      exit 1
    fi
  else
    mkdir -p "$(dirname "$OUT")"
    echo "$header" > "$OUT"
  fi

  for f in "${inputs[@]}"; do
    if [[ ! -f "$f" ]]; then
      echo "ERROR: bench file not found: $f" >&2
      exit 1
    fi
    read -r b_header < "$f" || true
    if [[ "$b_header" != "dataset,mode,step,variants,elapsed_s,rss_kb,throughput_vps,metric_a,metric_b" ]]; then
      echo "ERROR: $f header mismatch; expected bench.csv header." >&2
      exit 1
    fi
    awk -F, -v OFS=, \
      -v date="$DATE" -v git="$GIT" -v host="$HOST" -v cpu="$CPU" \
      -v threads="$THREADS" -v compressor="$COMPRESSOR" -v notes="$NOTES" \
      'NR==1 {next} NF>=9 {print date,git,host,cpu,threads,compressor,$1,$2,$3,$4,$5,$6,$7,$8,$9,notes}' \
      "$f" >> "$OUT"
  done
  echo "Merged ${#inputs[@]} file(s) into ${OUT}"
fi

if [[ -n "$DIFF" ]]; then
  if [[ -z "$BASELINE" || -z "$COMPARE" ]]; then
    if [[ ${#inputs[@]} -lt 2 ]]; then
      echo "ERROR: diff requires --baseline/--compare or at least two bench.csv inputs." >&2
      exit 1
    fi
    BASELINE="${inputs[0]}"
    COMPARE="${inputs[1]}"
  fi
  if [[ ! -f "$BASELINE" || ! -f "$COMPARE" ]]; then
    echo "ERROR: baseline/compare file not found." >&2
    exit 1
  fi
  echo "dataset,mode,step,metric,baseline,compare,delta,delta_pct" > "$DIFF"
  awk -F, -v OFS=, '
    function isnum(x) { return x ~ /^-?[0-9]+([.][0-9]+)?$/ }
    NR==1 { next }
    NR==FNR {
      key=$1 FS $2 FS $3
      base[key,"elapsed_s"]=$5
      base[key,"rss_kb"]=$6
      base[key,"throughput_vps"]=$7
      base[key,"metric_a"]=$8
      next
    }
    NR!=FNR {
      key=$1 FS $2 FS $3
      compare["elapsed_s"]=$5
      compare["rss_kb"]=$6
      compare["throughput_vps"]=$7
      compare["metric_a"]=$8
      for (m in compare) {
        b=base[key,m]
        c=compare[m]
        if (isnum(b) && isnum(c)) {
          d=c-b
          if (b==0) {
            pct="NA"
          } else {
            pct=(d/b)*100
          }
          printf "%s,%s,%s,%s,%s,%s,%.6f,%s\n", $1, $2, $3, m, b, c, d, pct
        }
      }
    }
  ' "$BASELINE" "$COMPARE" >> "$DIFF"
  echo "Diff written to ${DIFF}"
fi
