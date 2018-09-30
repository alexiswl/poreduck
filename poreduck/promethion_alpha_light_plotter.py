#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import humanfriendly
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig

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


def get_summary_files(summary_dir):
    summary_files = [os.path.join(summary_dir, summary_file)
                     for summary_file in os.listdir(summary_dir)
                     if summary_file.endswith(".txt")
                     and summary_file.startswith("sequencing")]
    return summary_files


def read_datasets(sequencing_summary_files):
    # Read in each of the csv files.
    dataset = [pd.read_csv(sequencing_summary_file, sep="\t", header=0)
                for sequencing_summary_file in sequencing_summary_files]

    # Merge the list of datasets
    dataset = merge_dataset(dataset)

    # Reset the dtypes for the time columns
    dataset = set_time_dtypes(dataset)

    # Sort the dataset by the template_start time
    dataset = sort_dataset(dataset)

    # Get yield column
    dataset['yield'] = get_yield(dataset)

    # Get pass column
    dataset['pass'] = get_pass(dataset)

    # Get duration ratio
    dataset['pore_speed'] = get_duration_ratio(dataset)

    # Get events ratio
    dataset['events_ratio'] = get_events_ratio(dataset)

    # Return the object
    return dataset


def merge_dataset(dataset):
    # Merge a list of datasets into a single pandas dataframe
    return pd.concat(dataset, sort=True, ignore_index=True)

def set_time_dtypes(dataset):

    # Return the time datasets as appropriate
    dataset.rename(columns={"start_time": "start_time_float",
                            "template_start": "template_time_float"}, inplace=True)

    # Add start_time_float template_time_float
    dataset['start_time_timedelta'] = pd.to_timedelta(dataset['start_time_float'], unit='s')
    dataset['template_start_time_timedelta'] = pd.to_timedelta(dataset['template_time_float'], unit='s')

    return dataset


def sort_dataset(dataset):
    # Sort the dataset by start_time (in seconds)
    return dataset.sort_values(['start_time_timedelta'])


def get_yield(dataset):
    # Get the yield datset
    return dataset.sequence_length_template.cumsum()


# plt.gcf().autofmt_xdate()
def get_pass(dataset):
    # Determine if sequence passed quality
    return dataset['mean_qscore_template'].apply(lambda x: True if x > 9 else False)


def get_duration_ratio(dataset):
    # Return the length in bases over time in seconds.
    return dataset.apply(lambda x: x.sequence_length_template/x.template_duration, axis='columns')


def get_events_ratio(dataset):
    # Return the events-per-base ratio
    return dataset.apply(lambda x: x.num_events / x.sequence_length_template, axis='columns')


# Plot yield
def plot_yield(dataset, name, plots_dir):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""
    # Plot total yield for the sample
    # Yield plot
    # Set up plotting structure
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Plot setting start_time_float as axis index
    dataset.set_index("start_time_float")["yield"].plot(ax=ax)
    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels
    ax.set_title("Read Distribution Graph for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative Yield")
    # Format nicely
    fig.tight_layout()
    savefig(os.path.join(plots_dir, "yield.png"))


# Plot histogram
def plot_hist(dataset, name, plots_dir):
    # Set globals
    num_bins = 50
    max_quantile = 0.999

    # Open up plotting frame
    plt.close('all')
    fig, ax = plt.subplots(1)

    # Get linspacing of histogram
    max_length = dataset['sequence_length_template'].quantile(max_quantile)
    trimmed = dataset.query("sequence_length_template < %s" % max_length)['sequence_length_template']

    # Develop the linspace
    bins = np.linspace(start=0, stop=trimmed.max(), num=num_bins)
    bin_width = bins[1] - bins[0]

    # Plot weighted histogram
    trimmed.plot(kind="hist", ax=ax, normed=1, bins=bins, alpha=0.6, weights=trimmed)

    # Set the axis formatters
    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        if y == 0:
            return 0
        s = humanfriendly.format_size(bin_width * trimmed.sum() * y, binary=False)
        return reformat_human_friendly(s)

    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable_seq))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))

    # Set titles
    ax.set_title("Read Distribution Graph for %s" % name)
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_ylabel("Bases per bin")
    ax.set_xlabel("Read length")

    # Ensure labels are not missed.
    fig.tight_layout()
    savefig(os.path.join(plots_dir, "hist.png"))


def reformat_human_friendly(s):
    """
    humanfriendly module returns with a few quirks
    1 = 1 byte ==> 1 b
    2 = 2 bytes ==> 2 bytes
    1000 = 1 KB ==> 1 Kb
    """
    s = s.replace(" byte", "")
    s = s.replace(" bytes", "")
    s = s.replace("B", "")
    s = s.replace("s", "")
    s += "b"
    return s


def y_yield_to_human_readable(y, position):
    # Convert distribution to base pairs
    if y == 0:
        return 0
    s = humanfriendly.format_size(y, binary=False)
    return reformat_human_friendly(s)


def x_yield_to_human_readable(x, position):
    # Convert time in seconds to hours or minutes
    hours = int(x // 3600)
    minutes = int((x % 3600) // 60)
    seconds = int(x % 60)
    if x == 0:
        return 0
    s = f"{hours:02d}:{minutes:02d}"
    return s


def x_hist_to_human_readable(x, position):
    # Convert distribution to base pairs
    if x == 0:
        return 0
    s = humanfriendly.format_size(x, binary=False)
    return reformat_human_friendly(s)


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
    plot_yield(dataset, args.name, args.plots_dir)
    plot_hist(dataset, args.name, args.plots_dir)

if __name__ == "__main__":
    main()
