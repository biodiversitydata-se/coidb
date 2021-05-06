#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse, write as write_fasta
from tempfile import NamedTemporaryFile
import subprocess
import sys
import os


def cluster_records(records, pid, threads):
    """
    Takes a list of sequence records, writes to a temporary file and
    clusters them with vsearch.

    :param records: Sequence records
    :param pid: Percent id to cluster by
    :return:
    """
    clustered_records = []
    if len(records) == 1:
        return records
    # Write to file
    f = NamedTemporaryFile(mode="w", delete=False)
    with open(f.name, 'w') as fh:
        write_fasta(records, fh, "fasta")
    f_null = open(os.devnull, 'w')
    cons_out = NamedTemporaryFile(mode="w", delete=False)
    # Run vsearch on tempfile
    subprocess.call(
        ['vsearch', '--cluster_fast', f.name, '--id',
         str(pid), '--consout', cons_out.name, '--notrunclabels', "--threads",
         str(threads)], stdout=f_null, stderr=f_null)
    # Read file with consensus sequences
    for record in parse(cons_out.name, 'fasta'):
        clustered_records.append(record)
    cons_out.close()
    f.close()
    return clustered_records


def get_species_clusters(f, pid, threads):
    """
    Iterates a fasta file sorted by species and clusters sequence for each
    species

    :param f: input fasta file
    :param pid: percent identity to cluster by
    :return: cluster dictionary (species as keys, records as values)
    """
    sp_groups = {}
    clusters = {}
    seqs = 0
    sp = ""
    sys.stderr.write(f"Reading {f} and clustering with vsearch...\n")
    for i, record in enumerate(parse(f, "fasta")):
        species = (record.description).split(";")[-1]
        if i == 0:
            sp = species
            sp_groups[species] = [record]
            continue
        # If the next record is from the same species, add it to the list
        if species == sp:
            sp_groups[species].append(record)
        # If not the same, attempt to cluster the stored sequences
        else:
            clusters[sp] = cluster_records(sp_groups[sp], pid, threads)
            seqs += len(clusters[sp])
            sp = species
            sp_groups[species] = [record]
    # Cluster the final sequences
    clusters[sp] = cluster_records(sp_groups[sp], pid, threads)
    seqs += len(clusters[sp])
    sys.stderr.write(f"Records read: {i+1}\n"
                     f"Species clustered: {len(clusters)}\n"
                     f"Sequence clusters: {seqs}\n")
    return clusters


def main(args):
    seq_clusters = get_species_clusters(args.fasta, args.pid, args.threads)
    with open(args.outfile, 'w') as fhout:
        for species, records in seq_clusters.items():
            write_fasta(records, fhout, "fasta")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", type=str, help="Fasta file input")
    parser.add_argument("outfile", type=str, help="Fasta file output")
    parser.add_argument("--pid", type=float, default=1.0,
                        help="Percent identity to cluster sequences by (1.0)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads to use (passed to vsearch) (1)")
    args = parser.parse_args()
    main(args)
