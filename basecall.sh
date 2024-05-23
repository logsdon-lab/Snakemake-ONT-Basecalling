#!/bin/bash

set -euo pipefail

# Check if command fails due to incomplete files
grep "IncompleteFilesException" <(snakemake -np $@ 2>&1 || true) > /dev/null && true
ret=$?

# Incomplete basecalling. Redo.
if [[ $ret -eq 0 ]]; then
    echo "Incomplete files detected. Redoing basecalling."
    snakemake -p --config rebasecall=True --rerun-incomplete --rerun-triggers mtime $@
else
    echo "Basecalling reads."
    snakemake -p --rerun-triggers mtime $@
fi
