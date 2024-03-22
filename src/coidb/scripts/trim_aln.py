#!/usr/bin/env python

"""
Trims an alignment to match an ASV sequence. This ASV sequence is assumed to have the name "ASV" in the alignment file
"""


from argparse import ArgumentParser
from Bio.SeqIO import parse
import sys
import pandas as pd

def read_seqs(f):
    records = {}
    with open(f, 'r') as fhin:
        for record in parse(fhin, "fasta"):
                seqid = record.description
                records[seqid] = record
    return records

def get_start(seq):
    start = 0
    for i, x in enumerate(seq):
        if x != "-":
            return i
    return start

def get_end(seq):
    end = -1
    for i, x in enumerate(seq[::-1]):
        if x != "-":
            return -i
    return end

def get_start_end(records):
    pos_count = {}
    for record in records.values():
        seq = str(record.seq)
        for i, x in enumerate(seq):
            if x != "-":
                if i not in pos_count:
                    pos_count[i] = 0
                pos_count[i] += 1
    df = pd.DataFrame(pos_count, index=["n"]).T
    df = df.sort_index()
    for row in df.iterrows():
        if row[1]["n"] == len(records):
            start = row[0]
            break
    df = df.sort_index(ascending=False)
    for row in df.iterrows():
        if row[1]["n"] == len(records):
            end = row[0]
            break
    return start, end


def main(args):
    records = read_seqs(args.fasta)
    if "ASV" in records.keys():
        asv = str(records["ASV"].seq)
        start = get_start(asv)
        end = get_end(asv)
    else:
        start, end = get_start_end(records)
    sys.stderr.write(f"Found start {start} and end {end}\n")
    del records["ASV"]
    sys.stderr.write(f"Trimming alignments\n")
    for seqid, record in records.items():
        subseq = str(record.seq[start:end])
        subseq = subseq.replace("-", "")
        if len(subseq) >= args.minlen:
            print(f">{seqid}\n{subseq}")
        else:
            sys.stderr.write(f"Skipping {seqid} because it is too short ({len(subseq)} < {args.minlen})\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", help="fasta file")
    parser.add_argument("-m", "--minlen", type=int, default=300, help="minimum length of subseq")
    args = parser.parse_args()
    main(args)
    