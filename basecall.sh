#!/bin/bash

set -euo pipefail

# Check if command fails due to incomplete files
grep "IncompleteFilesException" <(snakemake -np -c 1 $@ 2>&1 || true) > /dev/null && true
ret=$?
incomplete_files_list="/tmp/incompletefiles.list"
args_contain_dry_run=$(echo "Args: $*" | grep -- "-n" || true)

# Incomplete basecalling. Redo.
if [ $ret -eq 0 ] || [ -f "${incomplete_files_list}" ]; then
    echo "Incomplete files detected. Redoing basecalling."
    if [ ! -f "${incomplete_files_list}" ]; then
        # Generate list of incomplete files.
        grep -P "/(.*?)" <(snakemake -np -c 1 $@ 2>&1 || true) > $incomplete_files_list
    fi

    # Then pass as configuration.
    snakemake -p \
    --config incomplete_files=$incomplete_files_list \
    --keep-incomplete \
    --use-conda \
    --rerun-incomplete \
    --rerun-triggers mtime \
    -c 256 $@

    if [ -z "${args_contain_dry_run}" ]; then
        # Remove incomplete files list on completion.
        rm $incomplete_files_list
    fi
else
    echo "Basecalling reads."
    snakemake -p --keep-incomplete --rerun-triggers mtime -c 1 $@
fi
