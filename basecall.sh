#!/bin/bash

set -euo pipefail

WD=$(dirname $0)

# Check if command fails due to incomplete files
grep "IncompleteFilesException" <(snakemake -np --rerun-triggers mtime -c 1 $@ 2>&1 || true) > /dev/null && true
ret=$?
incomplete_files_list="${WD}/incompletefiles.txt"
priority_list="${WD}/priority.txt"

# Check that exists
touch "${priority_list}"
if [ -s "${priority_list}" ]; then
    PRIORITIES="-P "
else
    PRIORITIES=""
fi

# Incomplete basecalling. Redo.
if [ $ret -eq 0 ]; then
    echo "Incomplete files detected. Redoing basecalling."

    # Generate list of incomplete files.
    grep -P "^/(.*?)" <(snakemake -np --rerun-triggers mtime -c 1 $@ 2>&1 || true) > "${incomplete_files_list}"
else
    echo "Basecalling reads."
fi

# Then pass as configuration.
snakemake -p \
--config incomplete_files=$incomplete_files_list \
--keep-incomplete \
"${PRIORITIES}" \
--use-conda \
--rerun-incomplete \
--rerun-triggers mtime \
-c 256 $@

rm -f "${incomplete_files_list}"
