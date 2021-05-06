# COI DB

Download and format database for the COI gene (mitochondrial cytochrome oxidase 
subunit I).

## Overview
This tool downloads sequences + metadata from [GBIF](https://hosted-datasets.gbif.org/)
and formats sequences of interest for use with downstream metabarcoding analyses.

## Installation

1. Install with conda

```bash
conda install -c bioconda coidb
```

2. Clone git repository and build

```bash
git clone <repo>
python setup.py install
```

## Running the software

The `coidb` tool is a wrapper for a small snakemake workflow that downloads,
filters and clusters reference sequences for (mainly) the COI gene.

```
usage: coidb [-h] [-n] [-j CORES] [-f] [-u] [-c [CONFIG_FILE ...]] [--cluster-config CLUSTER_CONFIG] [--workdir WORKDIR] [-p] [targets ...]

positional arguments:
  targets               File(s) to create or steps to run. If omitted, the full pipeline is run.

optional arguments:
  -h, --help            show this help message and exit
  -n, --dryrun          Only print what to do, don't do anything [False]
  -j CORES, --cores CORES
                        Number of cores to run with [4]
  -f, --force           Force workflow run
  -u, --unlock          Unlock working directory
  -c [CONFIG_FILE ...], --config-file [CONFIG_FILE ...]
                        Path to configuration file
  --cluster-config CLUSTER_CONFIG
                        Path to cluster config (for running on SLURM)
  --workdir WORKDIR     Working directory. Defaults to current dir
  -p, --printshellcmds  Print shell commands
```

1. Download reference files

```bash
coidb download
```

This should produce two files `bold_info.tsv` and `bold_seqs.txt` containing
metadata and nucleotide sequences, respectively.

2. Filter reference files

```bash
coidb filter
```

This filters