#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import shutil
import os
import urllib.request
import datetime
import json

def api_match_species(bin_id):
    """
    This function uses the GBIF API to match a BIN id to a species name.
    It's only used as a fallback in case the add_species function does
    not work, for instance due to improper taxonomic info added for a BIN.

    :param bin_id: The BOLD bin id to match
    :return: Species name or NaN
    """
    try:
        with urllib.request.urlopen(f"https://api.gbif.org/v1/species/match?name={bin_id}") as response:
            json_text = response.read()
        response_dict = json.loads(json_text)
        return response_dict["species"]
    except:
        return np.nan


def add_species(species_bins, bin_tax_df, parent_df):
    """
    This function goes through putative species BINs, i.e. BINs that are
    classified down to genus level and which may have species annotations also.
    These have to be identified by looking up the name of the parent for each BIN.

    :param species_bins: list of BOLD BIN ids that may have species assignments
    :param bin_tax_df: the taxonomic dataframe
    :param parent_df: dataframe with parents to the species bins
    :return: the taxonmomic dataframe with species names added
    """
    for bold_id in tqdm(species_bins, unit=" BINs", desc="adding species"):
        parent = bin_tax_df.loc[bold_id, "parentNameUsageID"]
        if len(parent_df.loc[parent_df.taxonID == parent]) == 0:
            continue
        try:
            # Try to get the species name from the parent backbone
            bin_tax_df.loc[bold_id, "species"] = parent_df.loc[parent_df.taxonID==parent, "canonicalName"].values[0]
        except IndexError:
            # If that fails for some reason, do one attempt with the API
            bin_tax_df.loc[bold_id, "species"] = api_match_species(bold_id)
    return bin_tax_df


def fill_unassigned(df, ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]):
    """
    This function fills in 'the blanks' for each row in a dataframe
    For example:
    BOLD:ACQ8069 	Animalia 	Mollusca 	Gastropoda 	NaN 	Hermaeidae 	Mourgona 	NaN
    BOLD:AAN1572 	Animalia 	Mollusca 	Gastropoda 	NaN 	Hermaeidae 	NaN 	    NaN

    becomes:
    BOLD:ACQ8069 	Animalia 	Mollusca 	Gastropoda 	Gastropoda_X 	Hermaeidae 	Mourgona 	Mourgona_X
    BOLD:AAN1572 	Animalia 	Mollusca 	Gastropoda 	Gastropoda_X 	Hermaeidae 	Hermaeidae_X 	Hermaeidae_XX

    :param df: Dataframe to fill
    :param ranks: Ranks to iterate through
    :return: A dataframe with filled ranks
    """
    d = {}
    for row in tqdm(df.iterrows(), unit=" rows",
                         total=df.shape[0],
                         desc="filling unassigned ranks"):
        unknowns = 0
        for rank in ranks:
            # If an NaN entry is found
            if row[1][rank] != row[1][rank]:
                # Get the previous rank and its classification
                prev_rank = ranks.index(rank) - 1
                prev = ranks[prev_rank]
                # If this is the first time we're adding a suffix include the underscore
                if unknowns == 0:
                    row[1][rank] = row[1][prev] + "_X"
                else:
                    row[1][rank] = row[1][prev] + "X"
                unknowns += 1
            else:
                # Reset the unknowns counter in case there are intermediate unassigned ranks
                unknowns = 0
        d[row[0]] = row[1][ranks].to_dict()
    return pd.DataFrame(d).T


def logg(f):
    """
    Decorator for dataframe processing
    :param f:
    :return:
    """
    def wrapper(dataf, *args, **kwargs):
        """
        This wrapper outputs stats on running dataframe processing functions

        :param dataf: input dataframe
        :param args: arguments
        :param kwargs: keyword arguments
        :return: processed dataframes
        """
        # Get rows before processing
        rows_before = dataf.shape[0]
        # Get start time
        tic = datetime.datetime.now()
        # Perform the processing
        result = f(dataf, *args, **kwargs)
        # Get end time
        toc = datetime.datetime.now()
        # Get rows after processing
        rows_after = result.shape[0]
        sys.stderr.write(f"{toc-tic} seconds to {f.__name__}, {rows_before-rows_after} rows removed, {rows_after} rows remaining\n")
        return result
    return wrapper


@logg
def start(dataf):
    """
    Dummy function to generate a copy of the starting dataframe

    :param dataf: Starting dataframe
    :return: copy of dataframe
    """
    return dataf.copy()


@logg
def extract_bold_bins(dataf):
    """
    Returns the dataframe filtered to only rows with a BOLD BIN id, i.e.
    excluding NaN values

    :param dataf: Input dataframe
    :return: Filtered dataframe
    """
    return dataf.loc[dataf.bold_id==dataf.bold_id]


@logg
def fillna(dataf):
    """
    Fills NaN values in dataframe and returns it

    :param dataf: Input dataframe
    :return: Processed dataframe
    """
    return dataf.fillna("")


@logg
def filter_dataframe(dataf, filter_vals=[], filter_col=None):
    """
    Filter input dataframe based on column->value settings

    :param dataf: Input dataframe
    :param filter_vals: Value to filter with
    :param filter_col: Column in dataframe to filter on
    :return: Filtered dataframe
    """
    if len(filter_vals) > 0:
        sys.stderr.write(f"Filtering dataframe to {len(filter_vals)} items in {filter_col}\n")
        return dataf.loc[dataf[filter_col].isin(filter_vals)]
    return dataf

def write_seqs(seq_df, outfile, tmpfile, ranks):
    """
    Writes sequences to file, with taxonomic info in header

    :param seq_df: Sequence dataframe, including BOLD BIN ids and taxonomic info
    :param outfile: Output file to write to
    :param tmpfile: Temporary file
    :param ranks: list of ranks to write in header
    :return: Dataframe of taxonomic information that remains
    """
    old_index = []
    new_index = []
    # Sort sequences by BOLD IDs
    sys.stderr.write("Sorting sequences by BOLD IDs\n")
    seq_df = seq_df.sort_values("bold_id")
    tmpfile = os.path.expandvars(tmpfile)
    outfile = os.path.abspath(outfile)
    with open(tmpfile, 'w') as fhout:
        for record_id in tqdm(seq_df.index,
                              desc=f"Writing sequences to temporary directory",
                              unit=" seqs"):
            row = seq_df.loc[record_id]
            seq = seq_df.loc[record_id, "seq"]
            desc = ";".join([seq_df.loc[record_id, x] for x in ranks])
            fhout.write(f">{record_id} {desc}\n{seq}\n")
    sys.stderr.write(f"Moving {tmpfile} to {outfile}\n")
    shutil.move(tmpfile, outfile)
    return seq_df.drop("seq", axis=1)


def filter(sm):
    """
    Function to filter the BOLD information from GBIF

    :param sm:
    :return:
    """
    genes = sm.params.genes
    taxa = sm.params.filter_taxa
    filter_rank = sm.params.filter_rank
    ranks = sm.params.ranks
    nrows = sm.params.nrows
    ####################################
    ### Read and process occurrences ###
    ####################################
    sys.stderr.write(f"Reading {sm.input[0]}\n")
    occurrences = pd.read_csv(sm.input[0], sep="\t", usecols=[0, 4],
                     names=["record_id", "bold_id"],
                     dtype={'bold_id': str}, nrows=nrows)
    sys.stderr.write(f"{occurrences.shape[0]} records read\n")
    occurrences = (occurrences
        .pipe(start)
        .pipe(extract_bold_bins)
        .pipe(fillna))
    ##################################
    ### Read and process sequences ###
    ##################################
    sys.stderr.write(f"Reading {sm.input[1]}\n")
    seqs = pd.read_csv(sm.input[1], header=None, sep="\t", index_col=0,
                       names=["record_id", "gene", "seq"], usecols=[0,1,2],
                       nrows=nrows)
    sys.stderr.write(f"{seqs.shape[0]} records read\n")
    seqs = (seqs
            .pipe(start)
            .pipe(filter_dataframe, genes, "gene"))
    ##########################################
    ### Read and process backbone taxonomy ###
    ##########################################
    sys.stderr.write(f"Reading {sm.input[2]}\n")
    backbone = pd.read_csv(sm.input[2], header=0, sep="\t", nrows=nrows,
                           usecols=[0, 2, 5, 7, 11, 17, 18, 19, 20, 21, 22])
    sys.stderr.write(f"{backbone.shape[0]} lines read\n")
    backbone = (backbone
            .pipe(start)
            .pipe(filter_dataframe, taxa, filter_rank))
    #######################################
    ### Merge occurrences and sequences ###
    #######################################
    sys.stderr.write(
        f"Merging occurrences ({occurrences.shape[0]}) and sequence ({seqs.shape[0]}) records\n")
    seq_df = pd.merge(occurrences, seqs, left_on="record_id", right_index=True,
                      how="inner")
    sys.stderr.write(f"{seq_df.shape[0]} records remaining\n")
    ################################
    ### Remove duplicate records ###
    ################################
    sys.stderr.write(f"Removing duplicate records\n")
    seq_df_nr = seq_df.groupby("record_id").first().reset_index()
    sys.stderr.write(f"{seq_df.shape[0] - seq_df_nr.shape[0]} rows removed, {seq_df_nr.shape[0]} rows remaining\n")
    ##################################
    ### Assign taxonomy to records ###
    ##################################
    # Create new dataframe with scientific name as index
    bin_tax_df = backbone.set_index("scientificName")
    # Extract only BOLD IDs
    bin_tax_df = bin_tax_df.loc[bin_tax_df.index.str.startswith("BOLD:")]
    # Assign default species column
    bin_tax_df = bin_tax_df.assign(
        species=pd.Series([np.nan] * bin_tax_df.shape[0],
                          index=bin_tax_df.index))
    # Extract BOLD ids assigned down to genus level, putative species bins
    species_bins = list(
        bin_tax_df.loc[bin_tax_df.genus == bin_tax_df.genus].index)
    sys.stderr.write(f"{len(species_bins)} BINs assigned to genus level\n")
    # Extract ids of parents
    parent_ids = bin_tax_df.loc[species_bins, "parentNameUsageID"].values
    # Extract parent dataframe from parent ids
    parent_df = backbone.loc[(backbone.taxonID.isin(parent_ids))&(backbone.taxonRank=="species")]
    sys.stderr.write("Adding species names\n")
    # Attempt to add species to species_bins using parent dataframe
    bin_tax_df = add_species(species_bins, bin_tax_df, parent_df)
    bins_with_species_names = bin_tax_df.loc[bin_tax_df.species==bin_tax_df.species].shape[0]
    unique_species_names = len(bin_tax_df.loc[bin_tax_df.species==bin_tax_df.species, "species"].unique())
    sys.stderr.write(f"Added {unique_species_names} unique species names to {bins_with_species_names} BINS\n")
    # Fill unassigned ranks
    sys.stderr.write("Filling unassigned ranks\n")
    bin_tax_df = fill_unassigned(bin_tax_df, ranks)
    ################################################
    ### Merge BIN taxonomy with record dataframe ###
    ################################################
    sys.stderr.write(f"Merging BIN taxonomy with records\n")
    df = pd.merge(seq_df_nr, bin_tax_df.loc[:, ranks], left_on="bold_id", right_index=True)
    df.set_index("record_id", inplace=True)
    #####################
    ### Write to file ###
    #####################
    # Write seqs to file
    info_df = write_seqs(df, sm.output.fasta, sm.params.tmpf, sm.params.ranks)
    # Write info to file
    info_df.to_csv(sm.output.info, header=True, index=True, sep="\t")


def clean_fasta(sm):
    """
    Reformats fasta headers to strip vsearch specific strings

    :param sm: snakemake object
    :return:
    """
    from Bio import SeqIO
    import re
    with open(sm.input.fasta, 'r') as fhin, open(sm.output.fasta, 'w') as fhout:
        for record in SeqIO.parse(fhin, "fasta"):
            desc = (record.description).lstrip("centroid=")
            desc = re.split(";seqs=\d+", desc)[0]
            fhout.write(f">{desc}\n{record.seq}\n")


def format_fasta(sm):
    """
    Format a fasta file into two output files, one for use with the
    assignTaxonomy function and one with addSpecies in DADA2.

    The file for assignTaxonomy should have:
    >Level1;Level2;Level3;Level4;Level5;Level6;
    ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTC
    >Level1;Level2;Level3;Level4;Level5;
    CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCA

    while assignSpecies should have:
    >ID Genus species
    ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTC
    >ID Genus species
    CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCA

    for example:
    >Arthropoda;Insecta;Lepidoptera;Depressariidae;Acraephnes;
    and
    >AANIC216-11 Acraephnes nivea
    :param sm: snakemake object
    :return:
    """
    ranks = sm.params.ranks
    if "species" in ranks:
        ranks.remove("species")
    from Bio import SeqIO
    info = pd.read_csv(sm.input.info, sep="\t", index_col=0, header=0)
    with open(sm.output.assignTaxaFasta, 'w') as fh1, open(sm.output.addSpeciesFasta, 'w') as fh2:
        for record in SeqIO.parse(sm.input.fasta, "fasta"):
            names = []
            id = record.id
            rec_info = info.loc[id.lstrip("centroid=")]
            # Iterate ranks and add names as long as they are not NaN
            for r in ranks:
                n = rec_info[r]
                if n == n:
                    names.append(n)
                else:
                    break
            id_tax = ";".join(names)+";"
            fh1.write(f">{id_tax}\n{record.seq}\n")
            species = rec_info["species"]
            bold_id = rec_info["bold_id"]
            id_spec = f"{id} {species}"
            fh2.write(f">{id_spec}\n{record.seq}\n")


def main(sm):
    toolbox = {'filter_data': filter,
               'format': format_fasta,
               'clean': clean_fasta}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake) # noqa: F821
