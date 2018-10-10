#!/usr/bin/env python3

import argparse
import os

from promethion_alpha_light_plotter_helper import plot_data
from promethion_alpha_light_plotter_reader import get_summary_files
from promethion_alpha_light_plotter_reader import read_datasets

"""
Sequencing summary file columns
filename	read_id	run_id	channel	start_time	duration	num_events	template_start	num_events_template
template_duration	sequence_length_template	mean_qscore_template	strand_score_template
"""


def get_args():
    """Get arguments from commandline"""
    """
    Two simple arguments.
    1. Path to Sequencing Summary Directory
    """
    parser = argparse.ArgumentParser(description="Plot a run as it is going")
    parser.add_argument("--summary_dir", type=str, required=True,
                        help="Contains the txt files")
    parser.add_argument("--plots_dir", type=str, required=True,
                        help="Where do the plots go")
    parser.add_argument("--name", type=str, required=True,
                        help="Titles for plots")
    args = parser.parse_args()
    return args


def main():
    # Get args
    args = get_args()

    # Get summary files
    summary_files = get_summary_files(args.summary_dir)

    # Read in dataset
    dataset = read_datasets(summary_files)

    # Check plots_dir exists
    if not os.path.isdir(args.plots_dir):
        os.mkdir(args.plots_dir)

    # Plot yields and histograms
    plot_data(dataset, args.name, args.plots_dir)


if __name__ == "__main__":
    main()
