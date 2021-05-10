from importlib_resources import path as resource_path
from argparse import ArgumentParser
import yaml
import pandas as pd
from coidb.__main__ import run as coidb_run
from Bio import SeqIO


def populate_args():
    parser = ArgumentParser()
    args = parser.parse_args()
    args.targets = ['bold_filtered.fasta', 'bold_info_filtered.tsv']
    args.show_failed_logs = True
    args.workdir = "tests/data"
    args.dryrun = False
    args.cores = 1
    args.config_file = []
    args.cluster_config = None
    args.printshellcmds = False
    args.unlock = False
    args.force = True
    return args


def load_config(config_file):
    with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
    return config


def test_url():
    """
    Tests the default database url
    """
    from urllib import request
    with resource_path('coidb', "config.yaml") as configfile_path:
        config = load_config(configfile_path)
    url = config["database"]["url"]
    status_code = request.urlopen(url).getcode()
    assert status_code == 200


def test_filter():
    args = populate_args()
    assert coidb_run(args)


def test_filter_genes():
    args = populate_args()
    args.config_file = ["tests/config.filter_genes.yaml"]
    config = load_config(args.config_file[0])
    assert coidb_run(args)
    df = pd.read_csv(f"{args.workdir}/bold_info_filtered.tsv", sep="\t",
                     header=0)
    assert sorted(list(df.gene.unique())) == sorted(
        list(config["database"]["gene"]))


def test_filter_phyla():
    args = populate_args()
    args.config_file = ["tests/config.filter_phyla.yaml"]
    config = load_config(args.config_file[0])
    assert coidb_run(args)
    df = pd.read_csv(f"{args.workdir}/bold_info_filtered.tsv", sep="\t",
                     header=0)
    assert sorted(list(df.phylum.unique())) == sorted(
        list(config["database"]["phyla"]))


def test_cluster100():
    args = populate_args()
    args.config_file = ["tests/config.cluster_100.yaml"]
    args.workdir = "tests/clust_data"
    args.targets = ["bold_clustered.fasta"]
    assert coidb_run(args)
    fastafile = f"{args.workdir}/bold_clustered.fasta"
    records = list(SeqIO.parse(fastafile, "fasta"))
    assert len(records) == 3


def test_cluster95():
    args = populate_args()
    args.config_file = ["tests/config.cluster_95.yaml"]
    args.workdir = "tests/clust_data"
    args.targets = ["bold_clustered.fasta"]
    assert coidb_run(args)
    fastafile = f"{args.workdir}/bold_clustered.fasta"
    records = list(SeqIO.parse(fastafile, "fasta"))
    assert len(records) == 2
