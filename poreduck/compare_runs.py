#!/usr/bin/env python

import os
import platform
import pandas as pd
import matplotlib
if platform.system() == 'Linux':
    matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
import numpy as np
import sys

from poreduck.plot_yields import Read_set
from poreduck.plot_yields import x_hist_to_human_readable
from poreduck.plot_yields import y_yield_to_human_readable
from poreduck.plot_yields import x_yield_to_human_readable

CSV_DIR = ""
PLOTS_DIR = ""
RUNS = []
CLIP = False

"""
Take two runs, 
create a histogram and a yield comparison between the two runs.
pore duration to come.
Definitely scope for comparing quality as well.
"""


class Run:
    def __init__(self, fastq_dir, name):
        self.all_data = None
        self.yield_data = None
        self.read_sets = []
        self.fastq_dir = fastq_dir
        self.name = name

    def get_read_sets(self):
        # FASTQ_DIR used by readset class.
        fastq_files = [fastq_file for fastq_file in os.listdir(self.fastq_dir)
                       if fastq_file.endswith(".fastq")]
        for fastq_file in fastq_files:
            self.read_sets.append(Read_set(fastq_dir=self.fastq_dir, fastq_file=fastq_file))

    def get_fastq_data(self):
        for read_set in self.read_sets:
            read_set.read_fastq()

    def aggregate_dataframes(self):
        first_dataframe = True
        for read_set in self.read_sets:
            if first_dataframe:
                first_dataframe = False
                columns = list(read_set.df.columns)
                self.all_data = pd.DataFrame(columns=columns)
            self.all_data = self.all_data.append(read_set.df, ignore_index=True)
        self.all_data = self.all_data.sort_values(['time'], ascending=[True])

    def assign_yield_data(self):
        """We use this for both yield plots"""

        # Read in seq length and time from ALL_READS dataframe
        self.yield_data = self.all_data[['time', "seq_length"]]

        # Aggregate seqlength for each minute of sequencing. I love this resample command!
        self.yield_data.set_index(pd.DatetimeIndex(self.yield_data['time']), inplace=True)
        self.yield_data = self.yield_data.resample("1T").sum()
        self.yield_data.reset_index(inplace=True)#, drop=True)
        # Generate a cumulative sum of sequence data
        self.yield_data['cumsum_bp'] = self.yield_data['seq_length'].cumsum()
        # Convert time to timedelta format and then to float format, in hours.
        self.yield_data['duration_tdelta'] = self.yield_data['time'].apply(lambda t: t - self.yield_data['time'].min())
        self.yield_data['duration_float'] = self.yield_data['duration_tdelta'].apply(lambda t: t.total_seconds())


def plot_read_length_hist():
    seq_df_1 = RUNS[0].all_data["seq_length"]
    seq_df_2 = RUNS[1].all_data["seq_length"]
    # Define how many plots we want (1)
    fig, ax = plt.subplots(1)
    if CLIP:
        # Filter out the top 1000th percentile.
        seq_df_1 = seq_df_1[seq_df_1 < seq_df_1.quantile(0.9999)]
        seq_df_2 = seq_df_2[seq_df_2 < seq_df_2.quantile(0.9999)]
    # Set the axis formatters
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
    # Set labels of axis.
    ax.set_xlabel("Read length")
    ax.set_ylabel("")
    ax.get_xaxis().set_ticklabels([])

    # Plot the histogram
    ax.hist(seq_df_1, 50, weights=seq_df_1,
            normed=1, facecolor='blue', alpha=1, label=RUNS[0].name)
    ax.hist(seq_df_2, 50, weights=seq_df_2,
            normed=1, facecolor='red', alpha=0.5, label=RUNS[1].name)
    # Set the titles and add a legend.
    ax.set_title(f"Read Distribution Graph for {RUNS[0].name} and {RUNS[1].name}")
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    plt.legend()

    savefig(os.path.join(PLOTS_DIR, f"{RUNS[0].name}_{RUNS[1].name}_read_length_hist.png"))


def plot_yield_general():
    # Set subplots.
    fig, ax = plt.subplots(1)
    # Create ticks using numpy linspace. Ideally will create 6 points between 0 and 48 hours.
    num_points = 6
    min_x = min([RUNS[0].yield_data['duration_float'].min(), RUNS[1].yield_data['duration_float'].min()])
    max_x = max([RUNS[0].yield_data['duration_float'].max(), RUNS[1].yield_data['duration_float'].max()])
    x_ticks = np.linspace(min_x,
                          max_x,
                          num_points)
    ax.set_xticks(x_ticks)

    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels and limits.
    ax.set_xlabel("Duration (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_xlim(min_x, max_x)
    ax.set_title(f"Yield for {RUNS[0].name} and {RUNS[1].name} (B/Hour)")
    ax.plot(RUNS[0].yield_data['duration_float'], RUNS[0].yield_data['cumsum_bp'],
            linestyle="solid", markevery=[], label=RUNS[0].name)
    ax.plot(RUNS[1].yield_data['duration_float'], RUNS[1].yield_data['cumsum_bp'],
            linestyle="solid", markevery=[], label=RUNS[1].name)
    plt.legend()
    savefig(os.path.join(PLOTS_DIR, f"{RUNS[0].name}_{RUNS[1].name}_general_yield_plot.png"))


def set_args(args):
    global PLOTS_DIR, CLIP
    # Check to ensure that both fastq folders are there
    if not os.path.isdir(args.fastq_1):
        sys.exit(f"Error, could not find directory {args.fastq_1}")
    if not os.path.isdir(args.fastq_2):
        sys.exit(f"Error, could not find directory {args.fastq_2}")
    # Make plots directory if it doesn't exist
    if not os.path.isdir(args.plots_dir):
        os.mkdir(args.plots_dir)
    PLOTS_DIR = args.plots_dir
    if args.clip:
        CLIP = True


def get_runs(args):
    global RUNS
    RUNS.append(Run(args.fastq_1, args.name_1))
    RUNS.append(Run(args.fastq_2, args.name_2))
    # Now load up dataframes.
    for run in RUNS:
        run.get_read_sets()
        run.get_fastq_data()
        run.aggregate_dataframes()
        run.assign_yield_data()


def main(args):
    set_args(args)
    get_runs(args)
    plot_read_length_hist()
    plot_yield_general()

