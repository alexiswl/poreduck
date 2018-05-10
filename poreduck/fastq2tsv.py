#!/usr/bin/env python3

"""Simple script that grabs functions from other scripts to user the fastq headers to make a tsv dataframe"""

from albacore_server_scaled import get_fastq_dataframe
import argparse
import pandas as pd
import os
import sys


def get_args():
    parser = argparse.ArgumentParser(description="Get info from fastq header")
    parser.add_argument("--fastq_file", required=True, type=str)
    parser.add_argument("--output_file", required=True, type=str)
    parser.add_argument("--gzipped", action='store_true', default=False, dest='gzipped')
    return parser.parse_args()


def set_args(args):
    if not os.path.isdir(os.path.dirname(args.output_file)):
        os.mkdir(os.path.dirname(args.output_file))
    if not os.path.isfile(args.fastq_file):
        sys.exit("Fastq file %s does not exist" % args.fastq_file)


def main():
    args = get_args()
    set_args(args)
    df = get_fastq_dataframe(args.fastq_file, is_gzipped=args.gzipped)
    df.to_csv(args.output_file, index=False, header=True, sep="\t")


if __name__ == "__main__":
    main()