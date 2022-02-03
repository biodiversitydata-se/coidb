#!/usr/bin/env python

import pandas as pd
from tqdm import tqdm
import sys
import shutil
import os
import urllib.request
import json


def logg(f):
    def wrapper(dataf, *args, **kwargs):
        tic = datetime.datetime.now()
        rows_before = dataf.shape[0]
        result = f(dataf, *args, **kwargs)
        rows_after = result.shape[0]
        toc = datetime.datetime.now()
        sys.stderr.write(f"{toc-tic} seconds to {f.__name__}, {rows_before-rows_after} rows removed, {rows_after} rows remaining\n")
        return result
    return wrapper


@logg
def start(dataf):
    return dataf.copy()


@logg
def extract_bold_bins(dataf):
    return dataf.loc[dataf.bold_id==dataf.bold_id]


@logg
def fillna(dataf):
    return dataf.fillna("")


@logg
def filter_dataf(dataf, filter_vals=[], filter_col=None):
    if len(filter_vals) > 0:
        sys.stderr.write(f"Filtering dataframe to {len(filter_vals)} items in {filter_col}\n")
        return dataf.loc[dataf[filter_col].isin(filter_vals)]
    return dataf

def write_seqs(seq_df, outfile, tmpfile):
    ranks = ["phylum", "class", "order", "family", "genus", "species", "bold_id"]
    old_index = []
    new_index = []
    # Filter out sequences with non-DNA characters
    seq_df, dropped_ids = filter_non_standard(seq_df)
    sys.stderr.write(f"{dropped_ids} sequences dropped\n")
    # Sort sequences by BOLD IDs
    sys.stderr.write("Sorting sequences by BOLD IDs\n")
    seq_df = seq_df.sort_values("bold_id")
    tmpfile = os.path.expandvars(tmpfile)
    outfile = os.path.abspath(outfile)
    with open(tmpfile, 'w') as fhout:
        for record_id in tqdm(seq_df.index,
                              desc=f"Writing sequences to temporary directory",
                              unit=" seqs"):
            seq = seq_df.loc[record_id, "seq"]
            desc = ";".join([seq_df.loc[record_id, x] for x in ranks])
            fhout.write(f">{record_id} {desc}\n{seq}\n")
    sys.stderr.write(f"Moving {tmpfile} to {outfile}\n")
    shutil.move(tmpfile, outfile)
    return seq_df.drop("seq", axis=1)


def filter_non_standard(df):
    """
    Removes sequences with non-standard nucleotides

    :param df: Dataframe with fasta sequences
    :return: Dataframe with sequences with non-standard characters removed
    """
    drop_ids = []
    for record_id in tqdm(df.index, unit=" records",
                          desc="removing non-standard nucleotide seqs"):
        seq = df.loc[record_id, "seq"]
        seq = seq.replace("-", "").strip("N")
        letters = set([x for x in seq])
        for l in letters:
            if l.upper() not in ["A", "C", "G", "T"]:
                drop_ids.append(record_id)
                break
        df.loc[record_id, "seq"] = seq
    return df.drop(drop_ids), len(drop_ids)


def filter(sm):
    genes = sm.params.genes
    taxa = sm.params.filter_taxa
    filter_rank = sm.params.filter_rank
    # Read info
    sys.stderr.write(f"Reading info file {sm.input[0]}\n")
    df = pd.read_csv(sm.input[0], sep="\t", usecols=[0, 4],
                     names=["record_id", "bold_id"],
                     dtype={'bold_id': str})
    # Process info dataframe
    df = (df
        .pipe(start)
        .pipe(extract_bold_bins)
        .pipe(fillna))
    # Read fasta
    sys.stderr.write(f"Reading fasta file {sm.input[1]}\n")
    seqs = pd.read_csv(sm.input[1], header=None, sep="\t", index_col=0,
                       names=["record_id", "gene", "seq"], usecols=[0,1,2])
    # Process seq dataframe
    seqs = (seqs
            .pipe(start)
            .pipe(filter_dataf, genes, "gene"))
    # Process backbone
    sys.stderr.write(f"Reading GBIF backbone {sm.input[2]}\n")
    backbone = pd.read_csv(sm.input[2], header=0, sep="\t", nrows=None,
                           usecols=[0, 2, 5, 7, 11, 17, 18, 19, 20, 21, 22])
    backbone = (backbone
            .pipe(start)
            .pipe(filter_dataf, taxa, filter_rank))
    # Merge in order to get the intersection
    sys.stderr.write("Merging info and sequences\n")
    seq_df = pd.merge(df, seqs, left_on="record_id", right_index=True,
                      how="inner")
    sys.stderr.write(f"{seq_df.shape[0]} sequence records remaining\n")
    # Reindex seq df
    seq_df = seq_df.reset_index().drop("index", axis=1)
    # Count record ids
    rec_id_counts = seq_df.groupby("record_id").count()["gene"]
    multi = rec_id_counts.loc[rec_id_counts > 1]
    if multi.shape[0] > 0:
        for multi_record in tqdm(multi.index,
                             desc=f"Fixing {multi.shape[0]} non-unique ids",
                            unit=" ids"):
            recs = seq_df.loc[seq_df.record_id == multi_record]
            recs = handle_nonuniq(recs.copy())
            seq_df.loc[recs.index] = recs
    seq_df.set_index("record_id", inplace=True)
    # Write seqs to file
    info_df = write_seqs(seq_df, sm.output.fasta, sm.params.tmpf)
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
    ranks = ["phylum","class","order","family","genus"]
    from Bio import SeqIO
    info = pd.read_csv(sm.input.info, sep="\t", index_col=0, header=0)
    with open(sm.output.assignTaxaFasta, 'w') as fh1, open(sm.output.addSpeciesFasta, 'w') as fh2:
        for record in SeqIO.parse(sm.input.fasta, "fasta"):
            id = record.id
            rec_info = info.loc[id.lstrip("centroid=")]
            names = ["Eukaryota"]
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
            if species == species:
                # If species name includes more info than bold_id, add
                # species name as suffix
                if species != bold_id:
                    species = f"{bold_id} ({species})"
                id_spec = f"{id} {species}"
                fh2.write(f">{id_spec}\n{record.seq}\n")


def main(sm):
    toolbox = {'filter': filter,
               'format': format_fasta,
               'clean': clean_fasta}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake) # noqa: F821
