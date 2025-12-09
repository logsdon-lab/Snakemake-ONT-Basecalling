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

# Incomplete basecalling. Redo.
if [ $ret -eq 0 ]; then
    echo "Incomplete files detected. Redoing basecalling."

    # Generate list of incomplete files.
    grep -P "^/(.*?)" <(snakemake -np --rerun-triggers mtime -c 1 $@ 2>&1 || true) > "${incomplete_files_list}"
else
    echo "Basecalling reads."
fi

if [ -s "${priority_list}" ]; then
    # Construct glob. Expects format: "/data/{run_dir}/{subrun_dir}/{flowcell}/pod5/basecalling"
    # NOTE: This should check the configfile for the output format.
    # This script should be rewritten to check at some point.
    FILES_TO_PRIORITIZE=()
    while read -r line; do
        prun_dir="${line}";
        for prun_basecall_dir in $(realpath /data/"${prun_dir}"/*/*/pod5); do
            FILES_TO_PRIORITIZE+=("${prun_basecall_dir}/basecalling/${prun_dir}.bam")
        done
    done < "${priority_list}"
    # --prioritize these BAM files.
    PRIORITIES="-P ${FILES_TO_PRIORITIZE[@]}"
else
    PRIORITIES=""
fi

# Then pass as configuration.
# Don't expand PRIORITIES var in case empty.
snakemake -p \
--config incomplete_files=$incomplete_files_list \
--keep-incomplete ${PRIORITIES} \
--use-conda \
--rerun-incomplete \
--rerun-triggers mtime \
-c 256 $@

rm -f "${incomplete_files_list}"
