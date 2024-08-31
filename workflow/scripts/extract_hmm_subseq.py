#!/usr/bin/env python

"""
parse-hmmsearch.py

This script is meant to find a region on a protein sequence matching a region of interest in a hmm profile. Taking a tab-separated
file with translation coordinates used to translate DNA sequences, it will extract the DNA sequence for the region of interest.

Use-case:
Say you have ASVs from an amplicon study of a protein coding gene, e.g. COX1, and you a database of reference sequences that may or may 
not overlap with the ASV sequences. You want to extract the region on each reference sequence that matches the ASV sequence. If we know what
region of the HMM that the ASV sequence matches, we can use hmmsearch to find the same region on the reference sequence, then extract the
corresponding DNA region.

Parse hmmsearch output and extract the CDS sequence for a given region of the hmm profile. 

AATAAATAACATAAGATTTTGATTATTACCCCCTTCTTTATCTTTACTATTAATTAGAAGAATAGTTGAAACTGGAACAGGTACCGGATGAACTGTTTACCCA
.xxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyyxxxyyy

1. Get the seq_from and seq_to coordinates from the hmmsearch output
2. Calculate the start of the CDS sequence as cds_start = (seq_from - 1)*3+1
3. Calculate the end of the CDS sequence as cds_end = seq_to*3
4. Calculate the length of the CDS sequence as cds_len = cds_end - cds_start + 1
5. Read each untranslated DNA sequence from the fasta file
6. Read the translated coordinates
7. For each sequence, 
  a) if strand is -, reverse complement the DNA sequence
  b) get translation start and stop coordinates
  c) extract translated subseq as dna[translation_start-1:cds_len+1]
"""


from argparse import ArgumentParser
import pandas as pd
import tqdm
from Bio.SeqIO import parse
import sys

def read_seqs(f):
    records = {}
    with open(f, 'r') as fhin:
        for record in tqdm.tqdm(parse(fhin, "fasta"), unit=" records", desc=f"reading {f}", leave=False):
                seqid = record.description
                records[seqid] = record
    return records

def read_hmmout(f):
    """
    Reads the --domtblout output from hmmsearch
    """
    df = pd.read_csv(f, comment="#", sep=" +", engine="python", header=None, usecols=[0,2,5,15,16,17,18,19,20,22], names=["target","tlen","qlen","hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","description"])
    df["description"] = df.description.replace({"-": ""})
    df.index = df["target"]+" "+df["description"]
    df.rename(index=lambda x: x.rstrip(), inplace=True)
    df.drop(["target","description"], axis=1, inplace=True)
    return df

def get_hmm_region(df, hmm_from=78, hmm_to=202, shift_start=2, pad_end=5):
    """
    Filters each query sequence to only keep queries that have a match in the hmm region of interest
    aa sequence example
    Example:
    ASV aa seq:        MNNMSFWLLPPSLSLLLISSMVETGTGTGWTVYPPLSSIIAHTGSSVDFSIFSLHIAGIS
    ref aa seq: nnnnnnnMNNMSFWMLPPSFSLLLASSAVEAGVGTGWTVYPPLSSNIAHAGASVDLAIFSLHLAGVSnnnnnnnnn
    HMM:         -----{                     hmm region of interest                 }-----
    
    """
    df_filtered = df.loc[(df.hmm_from<=hmm_from)&(df.hmm_to>=hmm_to)]
    # distance from hmm start to hmm region of interest
    df_filtered = df_filtered.assign(aa_start_offset = hmm_from - df_filtered["hmm_from"])
    # distance from hmm end to hmm region of interest
    df_filtered = df_filtered.assign(aa_end_offset = df_filtered["hmm_to"] - hmm_to)
    # add start offset to env_from
    df_filtered = df_filtered.assign(aa_seq_from = df_filtered["env_from"] + df_filtered["aa_start_offset"])
    # subtract end offset from env_to
    df_filtered = df_filtered.assign(aa_seq_to = df_filtered["env_to"] - df_filtered["aa_end_offset"])
    # calculate start of CDS sequence
    df_filtered = df_filtered.assign(cds_start = (df_filtered["aa_seq_from"] - 1)*3+1+shift_start*3)
    # calculate end of CDS sequence
    df_filtered = df_filtered.assign(cds_end = df_filtered["aa_seq_to"]*3 + pad_end*3)
    # calculate length of CDS sequence
    df_filtered = df_filtered.assign(cds_len = df_filtered["cds_end"] - df_filtered["cds_start"] + 1)
    return df_filtered

def read_coords(f):
    coord = pd.read_csv(f, sep="\t", header=None, index_col=0, names=["seqid","translation_start","translation_stop","strand"])
    return coord

def get_subseq_coord(df, coord):
    df = pd.merge(df, coord, left_index=True, right_index=True)
    df["subseq_start"] = df["cds_start"]-df["translation_start"]
    df["subseq_end"] = df["subseq_start"]+df["cds_len"]
    return df

def main(args):
    # set params (for debugging)
    hmmout = args.hmmout
    dnafile = args.dna
    hmm_from = args.hmm_from
    hmm_to = args.hmm_to
    coordsfile = args.coords
    pad_end = args.pad_end
    shift_start = args.shift_start
    min_len = args.min_len
    df = read_hmmout(hmmout)
    df_filtered = get_hmm_region(df, hmm_from, hmm_to, shift_start,pad_end)
    coord = read_coords(coordsfile)
    df = get_subseq_coord(df_filtered, coord)
    dnarecs = read_seqs(dnafile)
    for seqid, d in df.iterrows():
        if d["strand"] == "-":
            seq = dnarecs[seqid].seq.reverse_complement()
        else:
            seq = dnarecs[seqid].seq
        subseq = seq[d["subseq_start"]:d["subseq_end"]]
        if len(subseq) >= min_len:
            sys.stdout.write(f">{seqid}\n{subseq}\n")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("hmmout", type=str, help="hmmsearch output file")
    parser.add_argument("coords", type=str, help="coordinates file")
    parser.add_argument("dna", type=str, help="DNA sequence file")
    parser.add_argument("-p", "--pad_end", type=int, default=5, help="when calculating end of cds, add this many amino acids to the end to account for insertions")
    parser.add_argument("-s", "--shift_start", type=int, default=0, help="when calculating start of cds, shift this many amino acids to the left")
    parser.add_argument("--hmm_from", type=int, default=78, help="hmm region start")
    parser.add_argument("--hmm_to", type=int, default=202, help="hmm region end")
    parser.add_argument("-m", "--min_len", type=int, default=300, help="minimum length of subseq")
    args = parser.parse_args()
    main(args)

