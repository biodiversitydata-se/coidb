#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import pandas as pd
from tqdm import tqdm


def read_info(f, replacements=None):
    info = pd.read_csv(f, sep="\t", index_col=0)
    d = {}
    if replacements:
        for item in replacements:
            key, value = item.split("=")
            d[key] = value
    return info.rename(columns=d)


def main(args):
    info = read_info(args.info, args.replace_rank)
    assert len(args.ranks) == len(
        set(info.columns).intersection(args.ranks)
    ), "not all ranks found in info file"
    with open(args.outfile, "w") as fhout:
        for record in tqdm(
            parse(args.fasta, "fasta"), unit=" records", desc="formatting records"
        ):
            taxlabels = [info.loc[record.id][x] for x in args.ranks]
            taxprefix = [x[0] for x in args.ranks]
            desc = ",".join(
                [f"{taxprefix[x]}:{taxlabels[x]}" for x in range(0, len(taxlabels))]
            )
            header = f">{record.id};tax={desc}"
            fhout.write(f"{header}\n{record.seq}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", type=str, help="Input fasta")
    parser.add_argument("info", type=str, help="Input info file")
    parser.add_argument("outfile", type=str, help="Output fasta")
    parser.add_argument(
        "--ranks",
        nargs="*",
        help="Ranks to assign",
        default=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
        choices=[
            "domain",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    )
    parser.add_argument(
        "--replace_rank", nargs="*", help="Replacements for header in info file"
    )
    args = parser.parse_args()
    main(args)
