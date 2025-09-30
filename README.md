# Snakemake-ONT-Basecalling
[![.github/workflows/main.yml](https://github.com/logsdon-lab/Snakemake-ONT-Basecalling/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/Snakemake-ONT-Basecalling/actions/workflows/main.yml)

Basecalling ONT reads workflow using [`dorado`](https://github.com/nanoporetech/dorado).

### Features
* Automatically detects and continues on incomplete files.
* Generates summary plots.

### Requirements
Requires a `dorado` executable.
```yaml
dorado:
  bin: "/home/prom/bin/dorado"
```

And `Snakemake`
```bash
snakemake -np
```

### Data
Expects the following directory structure where:
* `run_name` follows `run_dir_pattern` in `config.yaml`.

```
data/
└── run_name
    └── sample
        └── flowcell
            ├── pod5
            │   └── *.pod5
            └── sequencing_summary.txt
```

### Usage
```bash
./basecall.sh --config config/config.yaml -c 1
```

Each parameter after the script is passed to Snakemake.
