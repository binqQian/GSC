#!/usr/bin/env bash

export GSC_BIN=./build/gsc
export SMALL_DATA=
export MEDIUM_DATA=/home/binq/data/output_subset_5000_20000.vcf.gz
export LARGE_DATA=

# Optional knobs
export THREADS=32
export COMPRESSOR=bsc
export OUT_DIR=tmp/bench_$(date +%Y%m%d_%H%M%S)

# Query benchmarks
export VCF_DATA=
export GSC_FILE=
export RANDOM_SEED=42
export SAMPLE_SIZES="50 100 250 500 1000 2000"
export RANGE_SIZES_BP="1000 5000 10000 50000 100000"
export RANGE_REPEATS=2
export COMBINE_MODE=pair
