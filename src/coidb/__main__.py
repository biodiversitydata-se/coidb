"""
Download and format COI database
"""
from importlib_resources import path as resource_path
from argparse import ArgumentParser
import logging
import sys

from snakemake import snakemake
from snakemake.utils import available_cpu_count


class SnakemakeError(Exception):
    pass


def run(args):
    with resource_path("coidb", "Snakefile") as snakefile_path:
        forcerun = []
        if args.force:
            forcerun = args.targets
        success = snakemake(
            snakefile_path,
            targets=args.targets,
            dryrun=args.dryrun,
            cores=args.cores,
            configfiles=args.config_file,
            cluster_config=args.cluster_config,
            workdir=args.workdir,
            printshellcmds=args.printshellcmds,
            unlock=args.unlock,
            forcerun=forcerun,
            rerun_triggers="mtime",
        )
        return success
    if not success:
        raise SnakemakeError()


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "targets",
        nargs="*",
        default=[],
        help="File(s) to create or steps to run. If omitted, "
        "the full pipeline is run.",
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        action="store_true",
        help="Only print what to do, don't do anything [False]",
    )
    parser.add_argument(
        "-j", "--cores", type=int, default=4, help="Number of cores to run with [4]"
    )
    parser.add_argument("-f", "--force", action="store_true", help="Force workflow run")
    parser.add_argument(
        "-u", "--unlock", action="store_true", help="Unlock working directory"
    )
    parser.add_argument(
        "-c", "--config-file", nargs="*", default=[], help="Path to configuration file"
    )
    parser.add_argument(
        "--cluster-config",
        type=str,
        help="Path to cluster config (for running on SLURM)",
    )
    parser.add_argument(
        "--workdir", type=str, help="Working directory. Defaults to current dir"
    )
    parser.add_argument(
        "-p", "--printshellcmds", action="store_true", help="Print shell commands"
    )
    args = parser.parse_args()
    run(args)
