# COI reference sequences from BOLD DB

## General information

Author: SBDI molecular data team  
Contact e-mail: john.sundh@scilifelab.se
DOI: 10.17044/scilifelab.20514192
License: CC BY 4.0
Categories: Bioinformatics and computational biology not elsewhere classified, Computational ecology and phylogenetics, Ecology not elsewhere classified
Item type: Dataset
Keywords: COI sequence analysis, Ampliseq  
Funding: Swedish Research Council (VR), grant number 2019-00242. National Bioinformatic Infrastructure Sweden

This README file was last updated: 2022-08-23

Please cite as: Swedish Biodiversity Data Infrastructure (SBDI; 2022). COI reference sequences from BOLD DB

## Dataset description

This repository contains COI (mitochondrial cytochrome oxidase subunit I) sequences 
collected from the [BOLD database](https://boldsystems.org/). The fasta file
bold_clustered_cleaned.fasta.gz has record ids that can be queried in the [Public
Data Portal](https://boldsystems.org/index.php/Public_BINSearch?searchtype=records)
and each fasta header contains the taxonomic ranks + the BIN ID assigned to the
record. The taxonomic information for each record is also given in the tab-separated
file bold_info_filtered.tsv.gz.

The dataset was last created on February 18, 2022.

### Methods
The code used to generate this dataset consists of a snakemake workflow wrapped
into a python package that can be installed with [conda](https://docs.conda.io/en/latest/miniconda.html)
(`conda install -c bioconda coidb`).

Firstly sequence and taxonomic information for records in the BOLD database is 
downloaded from the [GBIF Hosted Datasets](https://hosted-datasets.gbif.org/ibol/).
This data is then filtered to only keep records annotated as 'COI-5P' and assigned
to a BIN ID. The taxonomic information is parsed in order to assign species names
and resolve higher level ranks for each BIN ID. Sequences are processed to remove
gap characters and leading and trailing `N`s. After this, any sequences with
remaining non-standard characters are removed.
Sequences are then clustered at 100% identity using [vsearch](https://github.com/torognes/vsearch) 
(Rognes _et al._ 2016). This clustering is done separately for sequences assigned 
to each BIN ID.   

## References

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584