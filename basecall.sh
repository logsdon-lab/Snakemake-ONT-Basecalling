#!/bin/bash

set -euo pipefail

# Check if command fails due to incomplete files
grep "IncompleteFilesException" <(snakemake -np --rerun-triggers mtime -c 1 $@ 2>&1 || true) > /dev/null && true
ret=$?
incomplete_files_list="/tmp/incompletefiles.list"

# Incomplete basecalling. Redo.
if [ $ret -eq 0 ]; then
    echo "Incomplete files detected. Redoing basecalling."

    # Generate list of incomplete files.
    grep -P "^/(.*?)" <(snakemake -np --rerun-triggers mtime -c 1 $@ 2>&1 || true) > $incomplete_files_list
else
    echo "Basecalling reads."
fi

# Then pass as configuration.
snakemake -p \
--config incomplete_files=$incomplete_files_list \
--keep-incomplete \
--use-conda \
--rerun-incomplete \
--rerun-triggers mtime \
-c 256 $@

rm -f "${incomplete_files_list}"
