#!/usr/bin/env python

"""
Trims an alignment to match an ASV sequence. This ASV sequence is assumed to have the name "ASV" in the alignment file
"""


from argparse import ArgumentParser
from Bio.SeqIO import parse
import sys

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


def main(args):
    records = read_seqs(args.fasta)
    asv = str(records["ASV"].seq)
    start = get_start(asv)
    end = get_end(asv)
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
    