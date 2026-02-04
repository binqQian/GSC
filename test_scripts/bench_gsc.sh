#!/usr/bin/env bash
set -euo pipefail

GSC_BIN="${GSC_BIN:-./build/gsc}"
OUT_DIR="${OUT_DIR:-tmp/bench_$(date +%Y%m%d_%H%M%S)}"
SMALL_DATA="${SMALL_DATA:-}"
MEDIUM_DATA="${MEDIUM_DATA:-}"
LARGE_DATA="${LARGE_DATA:-}"
THREADS="${THREADS:-}"
COMPRESSOR="${COMPRESSOR:-}"

mkdir -p "${OUT_DIR}"

if [[ ! -x "${GSC_BIN}" ]]; then
  echo "ERROR: GSC binary not found or not executable: ${GSC_BIN}" >&2
  exit 1
fi

if [[ ! -x /usr/bin/time ]]; then
  echo "ERROR: /usr/bin/time not found; required for RSS metrics." >&2
  exit 1
fi

log() {
  echo "[$(date +%H:%M:%S)] $*"
}

variant_count() {
  local f=$1
  if [[ "$f" == *.vcf.gz ]]; then
    zcat "$f" | grep -v '^#' | wc -l
  elif [[ "$f" == *.vcf ]]; then
    grep -v '^#' "$f" | wc -l
  else
    echo 0
  fi
}

raw_size_bytes() {
  local f=$1
  if [[ "$f" == *.vcf.gz ]]; then
    zcat "$f" | wc -c
  elif [[ "$f" == *.vcf ]]; then
    wc -c < "$f"
  else
    wc -c < "$f"
  fi
}

parse_elapsed_seconds() {
  local t="$1"
  local h=0 m=0 s=0
  IFS=':' read -r a b c <<<"$t"
  if [[ -n "${c:-}" ]]; then
    h=$a
    m=$b
    s=$c
  elif [[ -n "${b:-}" ]]; then
    m=$a
    s=$b
  else
    s=$a
  fi
  awk -v h="$h" -v m="$m" -v s="$s" 'BEGIN{printf "%.6f", (h*3600)+(m*60)+s}'
}

run_timed() {
  local label="$1"
  shift
  local time_log="${OUT_DIR}/${label}.time"
  local cmd_log="${OUT_DIR}/${label}.log"
  /usr/bin/time -v "$@" >"$cmd_log" 2>"$time_log"
  local elapsed
  elapsed=$(grep -F "Elapsed (wall clock) time" "$time_log" | awk -F': ' '{print $2}')
  local rss_kb
  rss_kb=$(grep -F "Maximum resident set size" "$time_log" | awk -F': ' '{print $2}')
  local elapsed_s
  elapsed_s=$(parse_elapsed_seconds "$elapsed")
  echo "${elapsed_s} ${rss_kb}"
}

bench_dataset() {
  local label="$1"
  local path="$2"
  if [[ -z "$path" ]]; then
    return
  fi
  if [[ ! -f "$path" ]]; then
    log "SKIP ${label}: missing ${path}"
    return
  fi

  log "Benchmark: ${label} -> ${path}"
  local variants
  variants=$(variant_count "$path")
  local raw_bytes
  raw_bytes=$(raw_size_bytes "$path")

  local base_opts=()
  if [[ -n "${THREADS}" ]]; then
    base_opts+=("-t" "${THREADS}")
  fi
  if [[ -n "${COMPRESSOR}" ]]; then
    base_opts+=("--compressor" "${COMPRESSOR}")
  fi

  local lossless_gsc="${OUT_DIR}/${label}_lossless.gsc"
  local lossless_vcf="${OUT_DIR}/${label}_lossless.vcf.gz"
  local lossy_gsc="${OUT_DIR}/${label}_lossy.gsc"
  local lossy_vcf="${OUT_DIR}/${label}_lossy.vcf.gz"

  local elapsed rss_kb

  # lossless compress
  read -r elapsed rss_kb < <(run_timed "${label}_lossless_compress" "${GSC_BIN}" compress -i "$path" -o "$lossless_gsc" "${base_opts[@]}")
  local comp_bytes
  comp_bytes=$(wc -c < "$lossless_gsc")
  local comp_ratio
  comp_ratio=$(awk -v a="$raw_bytes" -v b="$comp_bytes" 'BEGIN{printf "%.2f", (b/a)*100}')
  local c_vps
  c_vps=$(awk -v v="$variants" -v t="$elapsed" 'BEGIN{printf "%.2f", v/t}')
  printf "%s,lossless,compress,%s,%s,%s,%s,%s,%s\n" \
    "$label" "$variants" "$elapsed" "$rss_kb" "$c_vps" "$comp_bytes" "$comp_ratio" >> "${OUT_DIR}/bench.csv"

  # lossless decompress
  read -r elapsed rss_kb < <(run_timed "${label}_lossless_decompress" "${GSC_BIN}" decompress -i "$lossless_gsc" -o "$lossless_vcf" "${base_opts[@]}")
  local out_variants
  out_variants=$(variant_count "$lossless_vcf")
  local ok="OK"
  if [[ "$variants" -ne "$out_variants" ]]; then
    ok="MISMATCH"
  fi
  local d_vps
  d_vps=$(awk -v v="$out_variants" -v t="$elapsed" 'BEGIN{printf "%.2f", v/t}')
  printf "%s,lossless,decompress,%s,%s,%s,%s,%s,%s\n" \
    "$label" "$out_variants" "$elapsed" "$rss_kb" "$d_vps" "$ok" "$lossless_vcf" >> "${OUT_DIR}/bench.csv"

  # lossy compress
  read -r elapsed rss_kb < <(run_timed "${label}_lossy_compress" "${GSC_BIN}" compress -M -i "$path" -o "$lossy_gsc" "${base_opts[@]}")
  comp_bytes=$(wc -c < "$lossy_gsc")
  comp_ratio=$(awk -v a="$raw_bytes" -v b="$comp_bytes" 'BEGIN{printf "%.2f", (b/a)*100}')
  c_vps=$(awk -v v="$variants" -v t="$elapsed" 'BEGIN{printf "%.2f", v/t}')
  printf "%s,lossy,compress,%s,%s,%s,%s,%s,%s\n" \
    "$label" "$variants" "$elapsed" "$rss_kb" "$c_vps" "$comp_bytes" "$comp_ratio" >> "${OUT_DIR}/bench.csv"

  # lossy decompress
  read -r elapsed rss_kb < <(run_timed "${label}_lossy_decompress" "${GSC_BIN}" decompress -M -i "$lossy_gsc" -o "$lossy_vcf" "${base_opts[@]}")
  out_variants=$(variant_count "$lossy_vcf")
  ok="OK"
  if [[ "$variants" -ne "$out_variants" ]]; then
    ok="MISMATCH"
  fi
  d_vps=$(awk -v v="$out_variants" -v t="$elapsed" 'BEGIN{printf "%.2f", v/t}')
  printf "%s,lossy,decompress,%s,%s,%s,%s,%s,%s\n" \
    "$label" "$out_variants" "$elapsed" "$rss_kb" "$d_vps" "$ok" "$lossy_vcf" >> "${OUT_DIR}/bench.csv"
}

echo "dataset,mode,step,variants,elapsed_s,rss_kb,throughput_vps,metric_a,metric_b" > "${OUT_DIR}/bench.csv"

if [[ -n "${SMALL_DATA}" ]]; then
  bench_dataset "small" "${SMALL_DATA}"
else
  log "SKIP small: SMALL_DATA is empty"
fi
bench_dataset "medium" "${MEDIUM_DATA}"
bench_dataset "large" "${LARGE_DATA}"

log "Benchmark complete. Results: ${OUT_DIR}/bench.csv"
log "Tip: set MEDIUM_DATA/LARGE_DATA env vars to include bigger datasets."
