#!/bin/bash

set -euo pipefail

# Check if command fails due to incomplete files
grep "IncompleteFilesException" <(snakemake --summary 2>&1 || true) > /dev/null && true
ret=$?

# Incomplete basecalling. Redo.
if [[ $ret -eq 0 ]]; then
    echo "Incomplete files detected. Redoing basecalling."
    snakemake -c1 -np --config rebasecall=True
else
    echo "Basecalling reads."
    snakemake -c1 -np
fi