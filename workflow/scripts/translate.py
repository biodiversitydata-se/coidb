#!/usr/bin/env python

from Bio.SeqIO import parse
from argparse import ArgumentParser
import tqdm
import gzip
import sys


def translate(record, table=5):
    """
    Translate a record, test all six frames and return the longest translation without stop codons
    """
    translation = ""
    start = 0
    stop = 0
    strand = ""
    for j, seq in enumerate([record.seq, record.seq.reverse_complement()]):
        for i in [0,1,2]:
            # check if sequence is divisible by 3
            if len(seq[i:]) % 3 != 0:
                continue
            s = seq[i:].translate(table=table)
            if "*" not in s and len(s) > len(translation):
                translation = s
                start = i+1
                stop = len(s)*3+i
                if j == 0:
                    strand = "+"
                elif j == 1:
                    strand = "-"
    return translation, start, stop, strand

def main(args):
    """
    Read a fasta file and translate
    """
    coords = {}
    if (args.fasta).endswith(".gz"):
        openfn = gzip.open
    else:
        openfn = open
    with openfn(args.fasta, 'rt') as fhin:
        for record in tqdm.tqdm(parse(fhin, "fasta"), unit=" records", desc=f"reading {args.fasta}", leave=False):
            seqid = record.description
            aa, start, stop, strand = translate(record)
            if len(aa) > 0:
                sys.stdout.write(f">{seqid}\n{aa}\n")
                coords[seqid] = [start, stop, strand]
    if args.coord:
        with open(args.coord, 'w') as fhout:
            for seqid, coord in coords.items():
                fhout.write(f"{seqid}\t{coord[0]}\t{coord[1]}\t{coord[2]}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", type=str, help="Input fasta file")
    parser.add_argument("--coord", type=str, help="Output coordinates")
    args = parser.parse_args()
    main(args)
