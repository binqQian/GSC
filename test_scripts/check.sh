#!/usr/bin/env bash
BCFTOOLS=bcftools
if [ -x ./bcftools ]; then
  BCFTOOLS=./bcftools
fi
PYTHON_BIN=python3
CANON_HASH_SCRIPT="$(dirname "$0")/vcf_canon_hash.py"

cs1=$($BCFTOOLS view --no-version --threads 4 "$1" | "$PYTHON_BIN" "$CANON_HASH_SCRIPT")
cs2=$($BCFTOOLS view --no-version --threads 4 "$2" | "$PYTHON_BIN" "$CANON_HASH_SCRIPT")

echo "original checksum:$cs1"
echo "restore checksum:$cs2"

if [ $cs1 == $cs2 ] ; then 
    echo OK; 
else 
    echo NOT OK; 
    exit 1;
fi
