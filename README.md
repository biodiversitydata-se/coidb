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

## References

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584