#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ENV_FILE="${ENV_FILE:-${ROOT_DIR}/docs/bench/bench_env.sh}"
if [[ -f "${ENV_FILE}" ]]; then
  source "${ENV_FILE}"
fi

GSC_BIN="${GSC_BIN:-${ROOT_DIR}/build/gsc}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/tmp/bench_$(date +%Y%m%d_%H%M%S)}"
HISTORY="${HISTORY:-${ROOT_DIR}/docs/bench/bench_history.csv}"
TREND_DIR="${TREND_DIR:-${OUT_DIR}/trend}"
ALERT_DIR="${ALERT_DIR:-${OUT_DIR}/alert}"

TREND_ARGS_STR="${TREND_ARGS:-}"
ALERT_ARGS_STR="${ALERT_ARGS:-}"

read -r -a TREND_ARGS <<< "${TREND_ARGS_STR}"
read -r -a ALERT_ARGS <<< "${ALERT_ARGS_STR}"

GSC_BIN="${GSC_BIN}" OUT_DIR="${OUT_DIR}" "${SCRIPT_DIR}/bench_gsc.sh"

THREADS="${THREADS:-}" COMPRESSOR="${COMPRESSOR:-}" \
  "${SCRIPT_DIR}/bench_merge.sh" --out "${HISTORY}" "${OUT_DIR}/bench.csv"

"${SCRIPT_DIR}/bench_trend.py" --in "${HISTORY}" --out "${TREND_DIR}" "${TREND_ARGS[@]}"
"${SCRIPT_DIR}/bench_alert.py" --in "${HISTORY}" --out "${ALERT_DIR}" "${ALERT_ARGS[@]}"

echo "Pipeline complete:"
echo "  bench:  ${OUT_DIR}/bench.csv"
echo "  hist:   ${HISTORY}"
echo "  trend:  ${TREND_DIR}"
echo "  alert:  ${ALERT_DIR}"
