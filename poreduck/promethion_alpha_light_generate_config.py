#!/usr/bin/env python3

import argparse
import os
import yaml
import re
import pandas as pd
import json

#"""
#Generate a yaml file used to run each of the alpha_light python commands
#"""


def get_args():
    parser = argparse.ArgumentParser(description="Generate a yaml file")
    parser.add_argument('--sequencing_summary_dir',
                        help="Path to sequencing_summary folder", required=True)
    parser.add_argument('--fastq_dir',
                        help="Path to fastq files", required=True)
    parser.add_argument("--fast5_dir",
                        help="Path to folder of fast5 files", required=True)
    parser.add_argument("--flowcellID",
                        help="Flowcell ID to add to suffix of tar file", required=True)
    parser.add_argument("--rnumber",
                        help="Add random number to folder", required=True)
    parser.add_argument("--output_md5sum_fast5", required=True,
                        help="Where do you want the md5 for the fast5 files")
    parser.add_argument("--output_md5sum_fastq", required=True,
                        help="Where do you want the md5 for the fastq files")
    parser.add_argument("--output_yaml_file",
                        help="Yaml file to create", required=True)

    return parser.parse_args()


def get_all_files(sequencing_summary_dir, fastq_dir, fast5_dir):
    """
    :rtype: pd.DataFrame
    :param sequencing_summary_dir: string
    :param fastq_dir: string
    :param fast5_dir: string
    """
    sequencing_summary_files = [os.path.join(sequencing_summary_dir, sequencing_summary_file)
                                for sequencing_summary_file in os.listdir(sequencing_summary_dir)
                                if re.match('sequencing_summary_\d+.txt', sequencing_summary_file)]

    fastq_files = [os.path.join(fastq_dir, fastq_file)
                   for fastq_file in os.listdir(sequencing_summary_dir)
                   if re.match('fastq_\d+.fastq', fastq_file)]

    fast5_dirs = [os.path.join(fast5_dir, fast5_folder)
                  for fast5_folder in os.listdir(fast5_dir)
                  if os.path.isdir(os.path.join(fast5_dir, fast5_folder))
                  and re.match("^\d+$", fast5_folder)]

    sequencing_summary_df = pd.DataFrame(sequencing_summary_files, columns=["sequencing_summary_file"])
    fastq_df = pd.DataFrame(fastq_files, columns=["fastq_file"])
    fast5_df = pd.DataFrame(fast5_dirs, columns=["fast5_dir"])

    # Append number onto each dataframe
    sequencing_summary_df['number'] = sequencing_summary_df['sequencing_summary_file'].apply(
        lambda x: int(re.match("sequencing_summary_(\d+).txt", os.path.basename(x)).group(1)))
    fastq_df['number'] = fastq_df['fastq_file'].apply(
        lambda x: int(re.match("fastq_(\d+).fastq", os.path.basename(x)).group(1)))
    fast5_df['number'] = fast5_df['fast5_dir'].apply(
        lambda x: int(re.match('(\d+)', os.path.basename(x)).group(1)))

    # Sort dataframes by number
    sequencing_summary_df.sort_values(by=['number'], inplace=True)
    fastq_df.sort_values(by=['number'], inplace=True)
    fast5_df.sort_values(by=['number'], inplace=True)

    # Set number as index for each dataframe
    sequencing_summary_df.set_index("number", inplace=True)
    fastq_df.set_index("number", inplace=True)
    fast5_df.set_index("number", inplace=True)

    # Now merge each using pd.concat
    return pd.concat([sequencing_summary_df, fastq_df, fast5_df], axis='columns', join='inner', sort=True)


def output_yaml(yaml_file, dataset):
    with open(yaml_file, 'w') as file:
        yaml.dump(json.loads(dataset.to_json(orient='records')), file, default_flow_style=True)


def main():
    args = get_args()
    dataset = get_all_files(sequencing_summary_dir=args.sequencing_summary_dir,
                            fastq_dir=args.fastq_dir, fast5_dir=args.fast5_dir)

    # Get flowcell id
    dataset["FlowcellID"] = args.flowcellID

    # Get Rnumber
    dataset['rnumber'] = args.rnumber

    # Add md5sum outputs
    dataset['md5_fast5'] = args.output_md5sum_fast5
    dataset['md5_fastq'] = args.output_md5sum_fastq

    # Output the yaml file
    output_yaml(args.output_yaml_file, dataset)


if __name__ == "__main__":
    main()
