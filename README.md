# COI DB

## Overview
This tool downloads sequences + metadata from [GBIF](https://hosted-datasets.gbif.org/)
and formats sequences of interest for use with downstream metabarcoding analyses.

## Installation options

**Option 1**: Install with conda:

```bash
conda install -c bioconda coidb
```

**Option 2**: Download a release from the '**Releases**' section, unpack it and then create and activate the conda environment. Finally, install the software:

```bash
conda env create -f environment.yml
conda activate coidb
python -m pip install .
```

## Quick start
To see the steps that will be run, without actually running them, do:

```bash
coidb -n
```

Remove the `-n` flag to actually run the steps. 

## Output

The primary outputs of the tool are:

1. bold_clustered_cleaned.fasta: A fasta file with sequences clustered at whatever threshold set in the config file (default is 1.0 which means 100% identity). The header of each sequence in this file has the format

```
>GMGMN070-14 Animalia;Arthropoda;Insecta;Lepidoptera;Pieridae;Gonepteryx;Gonepteryx rhamni;BOLD:AAA9222
```

In this example `GMGMN070-14` is the representative id for the sequence and can be viewed in the BOLD database at https://www.boldsystems.org/index.php/Public_RecordView?processid=GMGMN070-14.

2. bold_clustered.sintax.fasta: This fasta file is compatible with the SINTAX classification tool implemented in [vsearch](https://github.com/torognes/vsearch) and has headers with the format:

```
>GMGMN070-14;tax=d:Animalia,k:Arthropoda,p:Insecta,c:Lepidoptera,o:Pieridae,f:Gonepteryx,g:Gonepteryx rhamni,s:BOLD:AAA9222
```

> [!NOTE]
> In the SINTAX formatted headers, the taxonomic ranks are shifted to allow classification down to BOLD_bin. Since SINTAX only allows for ranks prefixed with 'd' (for domain) 'k' (kingdom), 'p' (phylum), 'c' (class), 'o' (order), 'f' (family), 'g' (genus), or 's' (species) we shift the taxonomy so that kingdom becomes domain, etc., and prefix the BOLD bin id with 's'.

3. bold_clustered.assigntaxonomy.fasta and bold_clustered.addSpecies.fasta: These fasta files are compatible with the assignTaxonomy and addSpecies functions implemented in [DADA2](https://github.com/benjjneb/dada2/). For the assignTaxonomy file the headers have the format:

```
>Animalia;Arthropoda;Insecta;Lepidoptera;Pieridae;Gonepteryx;Gonepteryx rhamni;
```

and for the addSpecies file the headers have the format:

```
>GMGMN070-14 Gonepteryx rhamni
```

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
usage: coidb [-h] [-n] [-j CORES] [-f] [-u] [-c [CONFIG_FILE ...]] [--cluster-config CLUSTER_CONFIG] [--workdir WORKDIR] [-p] [-t]
             [targets ...]

positional arguments:
  targets               File(s) to create or steps to run. If omitted, the full pipeline is run.

options:
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
  -t, --touch           Touch output files (mark them up to date without really changing them) instead of running their commands.
```

## How it works

Firstly sequence and taxonomic information for records in the BOLD database is 
downloaded from the [GBIF Hosted Datasets](https://hosted-datasets.gbif.org/ibol/).
GBIF processes taxonomic information from BOLD in order to resolve ambiguous 
assignments for BOLD BINs. When there are conflicting assignments at a taxonomic 
rank an 80% consensus rule is applied to keep _e.g._ a species level assignment
if four out of five names in the BIN are equal [Kõljalg et al 2020](https://www.mdpi.com/2076-2607/8/12/1910/htm).
This data is then filtered to only keep records annotated as 'COI-5P' and assigned
to a BIN ID and duplicate entries are removed. 

#### Taxonomy
The taxonomic information obtained from GBIF is then parsed in order to extract
species names to BOLD BINs. This is done by:
1. find all BOLD BINs with a taxonomic assignment at genus level, these likely have
species names assigned from GBIF (see methods for species assignment [here](https://www.mdpi.com/2076-2607/8/12/1910/htm))
2. obtain all parent taxonomic ids for BOLD BINs from step 1 and use these to 
look up the species name for the BOLD BINs. 
3. For BOLD BINs where species name look-up failed in step 2, try to obtain 
species name using the [GBIF API](https://www.gbif.org/developer/summary).

The taxonomic data is then searched for rows where missing values for ranks are 
filled with the last known higher level rank, suffixed with `_X`. For instance,

| BOLD BIN     | kingdom   | phylum          | class | order       | family | genus | species |
|--------------|-----------|-----------------|-------|-------------|--------|-------|---------|
| BOLD:ACX1129 | Animalia  | Platyhelminthes | NaN   | Polycladida | NaN    | NaN   | NaN     |
| BOLD:ACX6548 | Chromista | Ochrophyta      | NaN   | NaN         | NaN    | NaN   | NaN     |

becomes:

| BOLD BIN     | kingdom   | phylum          | class             | order         | family         | genus           | species          |
|--------------|-----------|-----------------|-------------------|---------------|----------------|-----------------|------------------|
| BOLD:ACX1129 | Animalia  | Platyhelminthes | Platyhelminthes_X | Polycladida   | Polycladida_X  | Polycladida_XX  | Polycladida_XXX  |
| BOLD:ACX6548 | Chromista | Ochrophyta      | Ochrophyta_X      | Ochrophyta_XX | Ochrophyta_XXX | Ochrophyta_XXXX | Ochrophyta_XXXXX |

As you can see, an `X` is appended for each downstream rank with a missing assignment.

BOLD BINs are then screened for cases where there are more than 1 unique parent 
lineage for the same taxonomic assignment. For example, the following taxonomic 
information may be found for BOLD BINs with assignment 'Aphaenogaster' at the
genus level.

| kingdom  | phylym     | class       | order         | family        | genus         |
|----------|------------|-------------|---------------|---------------|---------------|
| Animalia | Animalia_X | Animalia_XX | Animalia_XXX  | Animalia_XXXX | Aphaenogaster |
| Animalia | Arthropoda | Insecta     | Hymenoptera   | Formicidae    | Aphaenogaster |               

A check is first made to see if unique parent lineages can be obtained by 
removing BINs that only have missing assignments for parent ranks up to and including 
phylum. If that doesn't result in a unique parent lineage, the conflicting rank
assignments are prefixed with the lowest assigned parent rank. 

For example, BOLD BINs with genus level assignment 'Paralagenidium' have both 
`k_Chromista;p_Oomycota;c_Peronosporea;o_Peronosporales;f_Pythiaceae` and 
`k_Chromista;p_Ochrophyta;c_Ochrophyta_X;o_Ochrophyta_XX;f_Ochrophyta_XXX` as parent
lineages. Since these conflicts cannot be resolved by removing BINs (all BINs have
assignments at phylum level), the taxa labels at genus and species level are prefixed
with either `Pythiaceae_` or `Ochrophyta_XXX_`.

#### Sequence processing
Sequences are then processed to remove gap characters and leading and trailing 
`N`s. After this, any sequences with remaining non-standard characters are removed.
Sequences are then clustered at 100% identity using [vsearch](https://github.com/torognes/vsearch) 
(Rognes _et al._ 2016). This clustering is done separately for sequences assigned 
to each BIN ID.   

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

#### Step 4: Clean headers

The `clean` step removes extra information from sequence headers generated as part of clustering. To run this step, do:

```bash
coidb clean
```

#### Step 5: Generate SINTAX/DADA2 formatted fasta

To also get the SINTAX and/or DADA2 formatted fasta file, do:

```bash
coidb format_sintax
```

or 

```bash
coidb format_dada2
```