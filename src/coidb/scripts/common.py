#!/usr/bin/env python

import pandas as pd
from tqdm import tqdm
import sys
import shutil
import os


def handle_nonuniq(recs):
    newids = []
    for i, r in enumerate(recs.record_id, start=1):
        newid = f"{r}_{i}"
        newids.append(newid)
    recs.loc[:, "record_id"] = newids
    return recs


def write_seqs(seq_df, outfile, tmpfile):
    ranks = ["phylum", "class", "order", "family", "genus", "species"]
    old_index = []
    new_index = []
    # Sort sequences by species
    sys.stderr.write("Sorting sequences by species\n")
    seq_df = seq_df.sort_values("species")
    tmpfile = os.path.expandvars(tmpfile)
    outfile = os.path.abspath(outfile)
    with open(tmpfile, 'w') as fhout:
        for record_id in tqdm(seq_df.index,
                              desc=f"Writing sequences to temporary directory",
                              unit=" seqs"):
            seq = seq_df.loc[record_id, "seq"]
            seq = seq.replace("-", "").strip("N")
            if "N" in seq:
                continue
            desc = ";".join([seq_df.loc[record_id, x] for x in ranks])
            fhout.write(f">{record_id} {desc}\n{seq}\n")
    sys.stderr.write(f"Moving {tmpfile} to {outfile}\n")
    shutil.move(tmpfile, outfile)
    return seq_df.drop("seq", axis=1)


def filter(sm):
    genes = sm.params.genes
    phyla = sm.params.phyla
    # Read info
    sys.stderr.write(f"Reading info file {sm.input[0]}\n")
    df = pd.read_csv(sm.input[0], sep="\t", usecols=[0, 4, 7, 8, 9, 10, 11, 12],
                     names=["record_id", "bold_id", "phylum", "class", "order",
                            "family", "genus", "species"],
                     dtype={'bold_id': str})
    df.fillna("", inplace=True)
    sys.stderr.write(f"{df.shape[0]} records read\n")
    if len(phyla) > 0:
        # Filter dataframe to phyla
        sys.stderr.write("Filtering info file to phyla of interest\n")
        df = df.loc[df.phylum.isin(phyla)]
        sys.stderr.write(f"{df.shape[0]} records remaining\n")
    # Read fasta
    sys.stderr.write(f"Reading fasta file {sm.input[1]}\n")
    seqs = pd.read_csv(sm.input[1], header=None, sep="\t", index_col=0,
                       names=["record_id", "gene", "seq"])
    sys.stderr.write(f"{seqs.shape[0]} sequences read\n")
    if len(genes) > 0:
        # Filter fasta to gene(s)
        sys.stderr.write("Filtering sequences to gene(s) of interest\n")
        seqs = seqs.loc[seqs.gene.isin(genes)]
        sys.stderr.write(f"{seqs.shape[0]} sequences remaining\n")
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


def main(sm):
    toolbox = {'filter': filter}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake) # noqa: F821
