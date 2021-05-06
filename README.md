# COI DB

Download and format database for the COI gene (mitochondrial cytochrome oxidase 
subunit I).

## Overview
This tool downloads sequences + metadata from [GBIF](https://hosted-datasets.gbif.org/)
and formats sequences of interest for use with downstream metabarcoding analyses.

## Installation

1. Install with conda:

```bash
conda install -c bioconda coidb
```

2. Download a release from the '**Releases**' section, unpack it and install:

```bash
python setup.py install
```

3. Clone git repository and build

```bash
git clone https://github.com/johnne/coidb.git
python setup.py install
```

### Dependencies

This tool depends on software defined in the `environment.yml` file. 

If you install the package with conda this also installs all requirements at the
same time. If you instead install a release or clone the repo, please make sure 
that the required software is available in your `PATH`.

## Quick start
To generate a file `bold_clustered.fasta` with COI-5P sequences, run:

```bash
coidb
```

This will download, filter and cluster sequences from [GBIF Hosted Datasets](https://hosted-datasets.gbif.org/ibol/).

See below for configuration and more options.

## Output

Sequences in the resulting `bold_clustered.fasta` fasta file contain the original 
identifier as their primary id, and a string showing their taxonomic lineage in 
the fasta header:

```bash
>centroid=GBA28357-15 Arthropoda;Insecta;Psocodea;Philotarsidae;Aaroniella;Aaroniella sp.;seqs=1
```

In this example `centroid=` indicates that sequences from this species were 
clustered with `vsearch` and that the representative sequence for the resulting 
cluster is `GBA28357-15`.

## Configuration
There are a few configurable parameters that modifies how sequences are filtered
and clustered. You can modify these parameters using a config file in `yaml` 
format. The default setup looks like this:

```yaml
database:
    # url to download info and sequence files from
    url: "https://hosted-datasets.gbif.org/ibol/ibol.zip"
    # gene of interest (will be used to filter sequences)
    gene:
        - "COI-5P"
    # phyla of interest (omit this in order to include all phyla)
    phyla: []
    # Percent identity to cluster seqs in the database by
    pid: 1.0
```

### Gene types
By default, only sequences named 'COI-5P' are included in the 
final output. To modify this behaviour you can supply a config file in `yaml`
format via `-c <path-to-configfile.yaml>`. For example, to also include 
'COI-3P' sequences you can create a config file, _e.g._ named `config.yaml` with 
these contents:

```yaml
database:
  gene:
    - 'COI-5P'
    - 'COI-3P' 
```

Then run `coidb` as:

```bash
coidb -c config.yaml
```

Typical gene names and their occurrence in the database are shown in 
[this table](#Common-gene-types).

### Phyla

The default is to include sequences from all taxa. However, you can filter the 
resulting sequences to only those from one or more phyla. For instance, to only
include sequences from the phyla 'Arthropoda' and 'Chordata' you supply a 
config file with these contents:

```yaml
database:
  phyla:
    - 'Arthropoda'
    - 'Chordata' 
```

Typical phyla and their occurrence in the database are shown in 
[this table](#Common-phyla).

### Clustering

After sequences have been filtered to the genes and phyla of interest they are
clustered on a per-species (or BOLD `BIN` id where applicable) basis using 
`vsearch`. By default this clustering is performed at 100% identity. To change
this behaviour, to _e.g._ 95% identity make sure your config file contains:

```yaml
database:
  pid: 0.95
```

## Command line options

The `coidb` tool is a wrapper for a small snakemake workflow that handles
all the downloading, filtering and clustering.

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

Explanation:

`-n, --dryrun`: Only print what will be done, don't actually do anything.

`-j, --cores`: The number of cores to run the workflow with. Because the download 
and filtering steps have to be run in sequential order this only affects the
clustering step using `vsearch`.

`-f, --force`: Force the execution of the workflow even though files already
exist.

`-u, --unlock`: Release a working directory lock (which could result from a
previously interrupted run)

`-c, --configfile`: Supply a configuration file to alter the behaviour of the
tool.

`--workdir`: Specify the directory in which to read/write output files. 
Defaults to the current directory.

`-p, --printshellcmds`: Shows the actual commands as they are being executed. 

### Step-by-step

You can also run the `coidb` tool in steps, _e.g._ if you are only interested
in some of the files or if you want to inspect the results before proceeding 
to the next step. This is done using the positional argument `targets`. 

Valid targets are `download`, `filter` and `cluster`. 

#### Step 1: Download
For example, to only
download files from GBIF you can run:

```bash
coidb download
```

This should produce two files `bold_info.tsv` and `bold_seqs.txt` containing
metadata and nucleotide sequences, respectively.

#### Step 2: Filter

To also filter the `bold_info.tsv` and `bold_seqs.txt` files (according to the 
default 'COI-5P' gene or any other genes/phyla you've defined in the optional 
config file) you can run:

```bash
coidb filter
```

This filters sequences in `bold_seqs.txt` and entries in `bold_info.tsv` to 
potential genes and phyla of interest, respectively. Entries are then merged
so that only sequences with relevant information are kept. Output files from
this step are `bold_filtered.fasta` and `bold_info_filtered.tsv`.


#### Step 3: Clustering

The final step clusters sequences in `bold_filtered.fasta` on a per-species 
basis. This means that for each species, the sequences are gathered, 
clustered with `vsearch` and only the representative sequences are kept. In this 
step sequences can either have a species name or a BOLD `BIN` ID 
(_e.g._ `BOLD:AAY5017`) and are treated as being equivalent.

To run the clustering step, do:

```bash
coidb cluster
```

The end result is a file `bold_clustered.fasta`.

## Common gene types

| #seqs   | gene   |
| ------- | ------ |
| 6074566 | COI-5P |
| 153409 | COI-3P |
| 146758 | ITS |
| 114124 | matK |
| 110798 | ITS2 |
| 86915 | rbcL |
| 66793 | rbcLa |
| 14192 | 16S |
| 13496 | CYTB |
| 10675 | trnH-psbA |
| 9787 | COII |
| 9140 | 28S |
| 9066 | COXIII |
| 6166 | ND2 |
| 5872 | ND1 |
| 5868 | ND5-0 |
| 5868 | ND3 |
| 5867 | ND4 |
| 5863 | ND4L |
| 5843 | ND6 |
| 5772 | ITS1 |
| 4866 | 28S-D2 |
| 3940 | 12S |
| 3751 | 18S |
| 3547 | atp6 |
| 3459 | 5-8S |
| 3135 | trnL-F |
| 3027 | D-loop |
| 2870 | EF1-alpha |
| 1991 | Wnt1 |
| 1822 | Rho |
| 1722 | COI-PSEUDO |
| 1716 | H3 |
| 1326 | CAD |
| 1241 | rpoC1 |
| 1236 | atpF-atpH |
| 968 | tufA |
| 944 | COI-LIKE |
| 865 | rpoB |
| 865 | UPA |
| 749 | psbK-psbI |
| 597 | 28S-D2-D3 |
| 470 | CAD4 |
| 449 | PSBA |
| 431 | PGD |
| 393 | DBY-EX7-8 |
| 383 | GAPDH |
| 353 | RpS5 |
| 336 | ycf1 |
| 309 | AATS |
| 296 | 28S-D1-D2 |
| 240 | 28S-D9-D10 |
| 223 | MDH |
| 194 | TPI |
| 192 | trnD-trnY-trnE |
| 188 | LWRHO |
| 186 | RAG1 |
| 172 | H4 |
| 171 | COII-COI |
| 168 | ND6-ND3 |
| 167 | RAG2 |
| 167 | 16S-ND2 |
| 166 | IDH |
| 154 | RpS2 |
| 144 | 18S-V4 |
| 141 | 28S-D3-D5 |
| 137 | RNF213 |
| 132 | MC1R |
| 132 | MB2-EX2-3 |
| 125 | fbpA |
| 124 | ND4L-MSH |
| 124 | ArgKin |
| 120 | CADH |
| 117 | CHD-Z |
| 107 | ENO |
| 103 | 28S-D3 |
| 101 | CHOLC |
| 99 | VDAC |
| 98 | ADR |
| 95 | RPB2 |
| 94 | atpB-rbcL |
| 94 | atp6-atp8 |
| 92 | DYN |
| 91 | H3-NUMT |
| 88 | COI-NUMT |
| 86 | PSA |
| 86 | CYTB-NUMT |
| 81 | AOX-fmt |
| 72 | trnK |
| 69 | matR |
| 65 | CsIV |
| 64 | nucLSU |
| 64 | EF2 |
| 61 | TYR |
| 61 | ARK |
| 56 | ATP1A |
| 55 | petD-intron |
| 55 | matK-trnK |
| 53 | PLAGL2 |
| 47 | psbA-3P |
| 38 | PER |
| 31 | matK-like |
| 31 | FL-COI |
| 30 | CAD1 |
| 30 | 18S-3P |
| 25 | rbcL-like |
| 24 | DDC |
| 21 | HfIV |
| 20 | R35 |
| 17 | COII-COIII |
| 16 | RBM15 |
| 16 | NGFB |
| 16 | CK1 |
| 15 | WSP |
| 14 | psaB |
| 14 | TULP |
| 10 | rpL32-trnL |
| 10 | PY-IGS |
| 9 | EF1-alpha-5P |
| 7 | NBC-COI-5P |
| 4 | COI-5PNMT1 |
| 2 | TMO-4C4 |
| 2 | PKD1 |
| 1 | S7 |
| 1 | RPL37 |
| 1 | RPB1 |
| 1 | RBCL-5P |
| 1 | COI-5PNMT2 |
| 1 | Beta-tubulin

## Common phyla

| #entries | Phylum |
| -------- | ------ | 
| 5886491 | Arthropoda |
| 505704 | Chordata |
| 270743 | Magnoliophyta |
| 180809 | Mollusca |
| 76301 | Ascomycota |
| 58536 | Annelida |
| 48727 | Basidiomycota |
| 29817 | Rhodophyta |
| 28723 | Echinodermata |
| 28105 | Platyhelminthes |
| 21786 | Nematoda |
| 19453 | Cnidaria |
| 16321 | Bryophyta |
| 9116 | Rotifera |
| 8368 | Pteridophyta |
| 7122 | Chlorophyta |
| 5863 | Pinophyta |
| 4877 | Porifera |
| 4863 | Heterokontophyta |
| 3770 | Nemertea |
| 3516 | Glomeromycota |
| 2934 | Zygomycota |
| 2095 | Acanthocephala |
| 1787 | Bryozoa |
| 1671 | Tardigrada |
| 1512 | Pyrrophycophyta |
| 1339 | Chaetognatha |
| 1248 | Onychophora |
| 952 | Lycopodiophyta |
| 711 | Gastrotricha |
| 640 | Sipuncula |
| 573 | Ciliophora |
| 393 | Kinorhyncha |
| 370 | Nematomorpha |
| 276 | Chytridiomycota |
| 273 | Cycliophora |
| 223 | Myxomycota |
| 202 | Brachiopoda |
| 153 | Ctenophora |
| 149 | Hemichordata |
| 104 | Priapulida |
| 102 | Phoronida |
| 61 | Chlorarachniophyta |
| 48 | Rhombozoa |
| 21 | Entoprocta |
| 16 | Xenacoelomorpha |
| 16 | Gnathostomulida |
| 12 | Placozoa