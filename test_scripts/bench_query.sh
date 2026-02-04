#!/usr/bin/env bash
set -euo pipefail

GSC_BIN="${GSC_BIN:-./build/gsc}"
OUT_DIR="${OUT_DIR:-tmp/bench_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-}"
COMPRESSOR="${COMPRESSOR:-}"

VCF_DATA="${VCF_DATA:-${MEDIUM_DATA:-${SMALL_DATA:-}}}"
GSC_FILE="${GSC_FILE:-}"

RANDOM_SEED="${RANDOM_SEED:-42}"
SAMPLE_SIZES="${SAMPLE_SIZES:-50 100 250 500 1000}"
RANGE_SIZES_BP="${RANGE_SIZES_BP:-1000 5000 10000 50000 100000}"
RANGE_REPEATS="${RANGE_REPEATS:-2}"
COMBINE_MODE="${COMBINE_MODE:-pair}"

mkdir -p "${OUT_DIR}"

if [[ ! -x "${GSC_BIN}" ]]; then
  echo "ERROR: GSC binary not found or not executable: ${GSC_BIN}" >&2
  exit 1
fi

if [[ ! -x /usr/bin/time ]]; then
  echo "ERROR: /usr/bin/time not found; required for RSS metrics." >&2
  exit 1
fi

if [[ -z "${VCF_DATA}" || ! -f "${VCF_DATA}" ]]; then
  echo "ERROR: VCF_DATA not set or missing: ${VCF_DATA}" >&2
  exit 1
fi

if [[ -z "${GSC_FILE}" ]]; then
  if [[ -f "${OUT_DIR}/medium_lossless.gsc" ]]; then
    GSC_FILE="${OUT_DIR}/medium_lossless.gsc"
  fi
fi

if [[ -z "${GSC_FILE}" || ! -f "${GSC_FILE}" ]]; then
  echo "ERROR: GSC_FILE not set or missing: ${GSC_FILE}" >&2
  exit 1
fi

log() {
  echo "[$(date +%H:%M:%S)] $*"
}

vcf_cat() {
  local f=$1
  if [[ "$f" == *.vcf.gz || "$f" == *.bgz || "$f" == *.bgzf ]]; then
    zcat "$f"
  else
    cat "$f"
  fi
}

variant_count() {
  local f=$1
  if [[ "$f" == *.vcf.gz || "$f" == *.bgz || "$f" == *.bgzf ]]; then
    zcat "$f" | grep -v '^#' | wc -l
  elif [[ "$f" == *.vcf ]]; then
    grep -v '^#' "$f" | wc -l
  else
    echo 0
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
  set +e
  /usr/bin/time -v "$@" >"$cmd_log" 2>"$time_log"
  local rc=$?
  set -e
  local elapsed
  elapsed=$(grep -F "Elapsed (wall clock) time" "$time_log" | awk -F': ' '{print $2}')
  local rss_kb
  rss_kb=$(grep -F "Maximum resident set size" "$time_log" | awk -F': ' '{print $2}')
  local elapsed_s
  if [[ -z "${elapsed:-}" ]]; then
    elapsed_s="0"
  else
    elapsed_s=$(parse_elapsed_seconds "$elapsed")
  fi
  echo "${elapsed_s} ${rss_kb:-0} ${rc}"
}

log "Preparing sample list from VCF header: ${VCF_DATA}"
all_samples_file="${OUT_DIR}/samples_all.txt"
python3 - "$VCF_DATA" "$all_samples_file" <<'PY'
import sys, gzip
vcf_path, out_path = sys.argv[1:]
opener = gzip.open if vcf_path.endswith(('.gz', '.bgz', '.bgzf')) else open
with opener(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            parts = line.rstrip('\n').split('\t')
            samples = parts[9:]
            with open(out_path, 'w') as out:
                out.write("\n".join(samples) + "\n")
            break
PY

sample_count=$(wc -l < "$all_samples_file")
if [[ "$sample_count" -eq 0 ]]; then
  echo "ERROR: failed to parse samples from VCF header." >&2
  exit 1
fi

log "Preparing positions list from VCF: ${VCF_DATA}"
positions_file="${OUT_DIR}/positions.tsv"
vcf_cat "$VCF_DATA" | grep -v '^#' | awk -F'\t' '{print $1"\t"$2}' > "$positions_file"

pos_count=$(wc -l < "$positions_file")
if [[ "$pos_count" -eq 0 ]]; then
  echo "ERROR: failed to parse positions from VCF." >&2
  exit 1
fi

log "Generating sample subsets"
mkdir -p "${OUT_DIR}/samples"

sample_sizes=()
for n in ${SAMPLE_SIZES}; do
  if [[ "$n" -gt 0 && "$n" -le "$sample_count" ]]; then
    sample_sizes+=("$n")
  fi
done

for n in "${sample_sizes[@]}"; do
  out_file="${OUT_DIR}/samples/samples_${n}.txt"
  python3 - "$all_samples_file" "$out_file" "$n" "$RANDOM_SEED" <<'PY'
import random, sys
all_path, out_path, n_str, seed_str = sys.argv[1:]
random.seed(int(seed_str) + int(n_str))
with open(all_path) as f:
    samples = [line.strip() for line in f if line.strip()]
random.shuffle(samples)
with open(out_path, 'w') as f:
    f.write("\n".join(samples[:int(n_str)]) + "\n")
PY
  log "  samples_${n}.txt"
done

log "Generating random ranges"
mkdir -p "${OUT_DIR}/ranges"

range_sizes=()
for bp in ${RANGE_SIZES_BP}; do
  if [[ "$bp" -gt 0 ]]; then
    range_sizes+=("$bp")
  fi
done

ranges_tsv="${OUT_DIR}/ranges/ranges.tsv"
python3 - "$positions_file" "$ranges_tsv" "$RANDOM_SEED" "$RANGE_REPEATS" "${range_sizes[@]}" <<'PY'
import random, sys
pos_path, out_path, seed_str, reps_str, *sizes = sys.argv[1:]
random.seed(int(seed_str))
reps = int(reps_str)
positions = []
with open(pos_path) as f:
    for line in f:
        chrom, pos = line.rstrip().split('\t')
        positions.append((chrom, int(pos)))
if not positions:
    raise SystemExit(1)
with open(out_path, 'w') as out:
    out.write("range_id\trange_bp\tchrom\tstart\tend\n")
    for size_str in sizes:
        size = int(size_str)
        for r in range(reps):
            chrom, pos = random.choice(positions)
            start = max(1, pos)
            end = start + size
            range_id = f"bp{size}_r{r+1}"
            out.write(f"{range_id}\t{size}\t{chrom}\t{start}\t{end}\n")
PY

log "Running query benchmarks"
base_opts=()
if [[ -n "${THREADS}" ]]; then
  base_opts+=("-t" "${THREADS}")
fi

query_csv="${OUT_DIR}/query.csv"
echo "test,kind,status,sample_count,range_bp,range,elapsed_s,rss_kb,variants,out_vcf" > "$query_csv"

# Sample-only queries
for n in "${sample_sizes[@]}"; do
  samples_file="${OUT_DIR}/samples/samples_${n}.txt"
  out_vcf="${OUT_DIR}/query_samples_${n}.vcf.gz"
  label="query_samples_${n}"
  read -r elapsed rss_kb rc < <(run_timed "$label" "$GSC_BIN" decompress -i "$GSC_FILE" -o "$out_vcf" --samples "@${samples_file}" "${base_opts[@]}")
  status="OK"
  out_variants=0
  if [[ "$rc" -ne 0 ]]; then
    status="FAIL"
  elif [[ -f "$out_vcf" ]]; then
    out_variants=$(variant_count "$out_vcf")
  fi
  printf "%s,samples,%s,%s,,,%s,%s,%s,%s\n" \
    "$label" "$status" "$n" "$elapsed" "$rss_kb" "$out_variants" "$out_vcf" >> "$query_csv"
done

# Range-only queries
while read -r range_id range_bp chrom start end; do
  if [[ "$range_id" == "range_id" ]]; then
    continue
  fi
  range_str="${chrom}:${start},${end}"
  out_vcf="${OUT_DIR}/query_range_${range_id}.vcf.gz"
  label="query_range_${range_id}"
  read -r elapsed rss_kb rc < <(run_timed "$label" "$GSC_BIN" decompress -i "$GSC_FILE" -o "$out_vcf" --range "$range_str" "${base_opts[@]}")
  status="OK"
  out_variants=0
  if [[ "$rc" -ne 0 ]]; then
    status="FAIL"
  elif [[ -f "$out_vcf" ]]; then
    out_variants=$(variant_count "$out_vcf")
  fi
  printf "%s,range,%s,,%s,%s,%s,%s,%s,%s\n" \
    "$label" "$status" "$range_bp" "$range_str" "$elapsed" "$rss_kb" "$out_variants" "$out_vcf" >> "$query_csv"
done < "$ranges_tsv"

# Combined range + samples
if [[ "${COMBINE_MODE}" == "cross" ]]; then
  for n in "${sample_sizes[@]}"; do
    samples_file="${OUT_DIR}/samples/samples_${n}.txt"
    while read -r range_id range_bp chrom start end; do
      if [[ "$range_id" == "range_id" ]]; then
        continue
      fi
      range_str="${chrom}:${start},${end}"
      out_vcf="${OUT_DIR}/query_range_samples_${range_id}_n${n}.vcf.gz"
      label="query_range_samples_${range_id}_n${n}"
      read -r elapsed rss_kb rc < <(run_timed "$label" "$GSC_BIN" decompress -i "$GSC_FILE" -o "$out_vcf" --range "$range_str" --samples "@${samples_file}" "${base_opts[@]}")
      status="OK"
      out_variants=0
      if [[ "$rc" -ne 0 ]]; then
        status="FAIL"
      elif [[ -f "$out_vcf" ]]; then
        out_variants=$(variant_count "$out_vcf")
      fi
      printf "%s,range+samples,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$label" "$status" "$n" "$range_bp" "$range_str" "$elapsed" "$rss_kb" "$out_variants" "$out_vcf" >> "$query_csv"
    done < "$ranges_tsv"
  done
else
  i=0
  mapfile -t range_lines < <(tail -n +2 "$ranges_tsv")
  range_count=${#range_lines[@]}
  for n in "${sample_sizes[@]}"; do
    if [[ "$range_count" -eq 0 ]]; then
      break
    fi
    line="${range_lines[$((i % range_count))]}"
    IFS=$'\t' read -r range_id range_bp chrom start end <<< "$line"
    range_str="${chrom}:${start},${end}"
    samples_file="${OUT_DIR}/samples/samples_${n}.txt"
    out_vcf="${OUT_DIR}/query_range_samples_${range_id}_n${n}.vcf.gz"
    label="query_range_samples_${range_id}_n${n}"
    read -r elapsed rss_kb rc < <(run_timed "$label" "$GSC_BIN" decompress -i "$GSC_FILE" -o "$out_vcf" --range "$range_str" --samples "@${samples_file}" "${base_opts[@]}")
    status="OK"
    out_variants=0
    if [[ "$rc" -ne 0 ]]; then
      status="FAIL"
    elif [[ -f "$out_vcf" ]]; then
      out_variants=$(variant_count "$out_vcf")
    fi
    printf "%s,range+samples,%s,%s,%s,%s,%s,%s,%s,%s\n" \
      "$label" "$status" "$n" "$range_bp" "$range_str" "$elapsed" "$rss_kb" "$out_variants" "$out_vcf" >> "$query_csv"
    i=$((i + 1))
  done
fi

log "Query benchmark complete. Results: ${query_csv}"
