#!/usr/bin/env python3

import argparse
import yaml
import json
import pandas as pd
import subprocess
import logging
import os

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def get_args():
    parser = argparse.ArgumentParser(description="Get the a yaml file created in the config script")
    parser.add_argument("--config", required=True,
                        help="Path to config file")
    return parser.parse_args()


def read_config(config):
    with open(config) as f:
        config_data = yaml.load(f)
    return config_data


def run_process(config_data):
    here = os.path.dirname(os.path.realpath(__file__))
    tar_proc = subprocess.run([os.path.join(here, "promethion_alpha_light.py"),  # Then come the options.
                               "--sequencing_summary_path=%s" % config_data.sequencing_summary_file,
                               "--fastq_path=%s" % config_data.fastq_file,
                               "--fast5_path=%s" % config_data.fast5_dir,
                               "--flowcellID=%s" % config_data.FlowcellID,
                               "--rnumber=%s" % config_data.rnumber,
                               "--md5_fast5=%s" % config_data.md5_fast5,
                               "--md5_fastq=%s" % config_data.md5_fastq,
                               "--inplace"],
                              capture_output=True)

    if tar_proc.returncode == 0:
        logging.warning("Process completed successfully")
        logging.warning("Stdout = %s" % tar_proc.stdout.decode())
        logging.warning("Stderr = %s" % tar_proc.stderr.decode())
    else:
        logging.warning("Process returned non-zero exit code.")
        logging.warning("Stdout = %s" % tar_proc.stdout.decode())
        logging.warning("Stderr = %s" % tar_proc.stderr.decode())


def main():
    args = get_args()
    dataframe = pd.DataFrame(read_config(args.config))

    for row in dataframe.itertuples():
        run_process(row)


if __name__ == "__main__":
    main()    
