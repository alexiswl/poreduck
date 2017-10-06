#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse
import sys
import os
import statistics
import matplotlib
import platform
if platform.system() == 'Linux':
    matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import humanize
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
from itertools import chain
import seaborn as sns
import time


"""
This script will create a yield plot of the data that has been created by the
albacore_server_scaled.py script and placed in the folder transfer_fast5_to_server.py
As it takes some time to handle and format the data such that the yield plots can be used...
we might as well make a bunch of plots to go along with it!
"""

# Before we begin, are we using python 3.6 or greater?
try:
    assert sys.version_info >= (3, 6)
except AssertionError:
    sys.exit("Error: Python version out of date. Require 3.6 or higher.")

# Set global variables
CSV_DIR = ""
FASTQ_DIR = ""
PLOTS_DIR = ""
CWD = os.path.abspath(os.getcwd())
CSV_FILES = []
FASTQ_FILES = []
READ_SETS = {}
ALL_READS = None
DATEPARSE_CSV = lambda dates: [pd.datetime.strptime(d, '%a %b %d %H:%M:%S %Y') for d in dates]
DATEPARSE_FASTQ = lambda dates: [pd.datetime.strptime(d, "%Y-%m-%dT%H:%M:%SZ") for d in dates]
SAMPLE_NAME = ""
YIELD_DATA = None
QUALITY_DESCRIPTIONS = ["bad", "alright", "better", "best"]
QUALITY_COLOURS = ['#e51400', '#fa6800', '#a4c400', '#60A917']
CLIP = False

# Import arguments.
"""
    csv directory
    fastq directory
    output directory
"""


class Read_Set:
    """
    The class readset is assigned to each set of 4000 reads
    It contains attributes such as fastq file, csv file, random_number and mux boolean.
    """
    def __init__(self, fastq_dir, fastq_file):
        # Use fastq file to get meta data info
        self.fastq_path = os.path.join(fastq_dir, fastq_file)
        self.fastq_file = fastq_file
        self.number = None
        self.random = None
        self.id = None
        self.csv_df = None

        if not CSV_DIR == "":
            self.number = self.fastq_file.split("_")[0]
            self.random = self.fastq_file.split("_")[1]
            self.csv_file = [csv_file for csv_file in os.listdir(CSV_DIR)
                             if csv_file.split("_")[0] == self.number
                             and csv_file.split("_")[1] == self.random
                             and csv_file.endswith(".csv")][0]
            self.csv_path = os.path.join(CSV_DIR, self.csv_file)

            if "mux" in fastq_file:
                self.id = self.number + self.random + "_mux"
            else:
                self.id = self.number + self.random + "_seq"
        self.df = None
        self.added_fastq_data = False
        self.added_csv_data = False
        self.aggregated_to_global_dataframe = False

    def read_csv(self):
        self.df = pd.read_csv(self.csv_path, header=0, parse_dates=['ctime'], date_parser=DATEPARSE_CSV)

    def read_fastq(self):
        if not os.path.isfile(self.fastq_path):
            print("Could not find path")
            return
        # Use the SeqIO package to iterate through the records:
        # Fastq header looks like this.
        # 'b1814d98-a01e-4fc6-a32c-60e1a062a957 runid=0608745933a900777aad6c8d9636f25227fe54a1 read=111140
        #  ch=351 start_time=2017-08-24T07:04:36Z
        # Create the columns we will write to fastq_id, and seq_length
        self.df = pd.DataFrame(data=None, columns=["fastq_id", "read", "channel", "time", "seq_length", "av_qual"])
        # Run through fastq file and add attributes to dataframe.
        for record in SeqIO.parse(self.fastq_path, "fastq"):
            fastq_id = record.id.split()[0]
            row_as_dict = dict(x.split("=") for x in record.description.split()[1:])
            # Get length and quality of read
            fastq_length = len(record.seq)
            fastq_quality = statistics.mean(record.letter_annotations["phred_quality"])
            # Order row in order of fastq, read, channel, time, length, quality
            row_as_list = [fastq_id, row_as_dict["read"], row_as_dict["ch"],
                           row_as_dict["start_time"],
                           fastq_length, fastq_quality]
            self.df.loc[-1] = row_as_list
            self.df.index += 1
        # Convert columns to correct format
        self.df['seq_length'] = pd.to_numeric(self.df["seq_length"])
        self.df['av_qual'] = pd.to_numeric(self.df['av_qual'])
        # self.df['time'] = pd.to_datetime(self.df['time'], format="%Y-%m-%dT%H:%M:%SZ")
        self.df['time'] = self.df['time'].apply(lambda x: pd.to_datetime(x, format="%Y-%m-%dT%H:%M:%SZ"))
        # Tick the box that we have added the fastq to a dataframe
        self.added_fastq_data = True

    def append_csv_data(self):
        # Set default columns
        self.df["mux"] = self.df.apply(lambda _: 0, axis=1)
        self.df["duration"] = self.df.apply(lambda _: 0, axis=1)
        muxs = {}
        durations = {}
        # Open csv file: Columns are filename, channel, read_no, ctime, mux number
        self.csv_df = pd.read_csv(self.csv_path, header=0)
        for csv_row in self.csv_df.itertuples():
            # Find the index of the respective csv file
            channel_csv = str(csv_row.channel)
            read_csv = str(csv_row.read_no)
            try:
                df_index = self.df.query("channel==@channel_csv & read==@read_csv").index.tolist()[0]
            except IndexError: 
                continue
            # Write value to dictionary with index of our fastq dataframe as the key
            muxs[df_index] = csv_row.mux
            durations[df_index] = csv_row.duration
        # Now write the value of the mux and duration to fastq dataframe
        for index, mux in muxs.items():
            self.df.set_value(index, "mux", mux)
        for index, duration in durations.items():
            self.df.set_value(index, "duration", durations)


def get_arguments():
    parser = argparse.ArgumentParser(
        description="This plot takes in the csv files along with the fastq files to produce a yield plot of the data")
    parser.add_argument("--csv_dir", type=str, required=True,
                        help="/path/to/csv_dir" +
                             "Should have a bunch of csv files in it.")
    parser.add_argument("--fastq_dir", type=str, required=True,
                        help="/path/to/fastq/files")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="/path/to/plots_dir" +
                             "By default will be created in current working directory")
    parser.add_argument("--sample_name", type=str, required=False,
                        help="Name to add onto each of the plots")
    args = parser.parse_args()
    return args


def set_arguments(args):
    global CSV_DIR, FASTQ_DIR, PLOTS_DIR
    global CSV_FILES, SAMPLE_NAME, CLIP
    if not args.no_csv:
        CSV_DIR = args.csv_dir
    FASTQ_DIR = args.fastq_dir
    if args.output_dir:
        PLOTS_DIR = args.output_dir
        if not os.path.isdir(os.path.abspath(os.path.join(PLOTS_DIR, os.pardir))):
            sys.exit(f"Error: Cannot create {PLOTS_DIR}. Parent directory does not exist.")
    else:
        PLOTS_DIR = os.path.join(CWD, "plots")
    if not os.path.isdir(PLOTS_DIR):
        os.mkdir(PLOTS_DIR)
    # Check and set CSV_DIR
    if not CSV_DIR == "":
        if not os.path.isdir(CSV_DIR):
            sys.exit(f"Error: {CSV_DIR} could not be found")
        CSV_DIR = os.path.abspath(CSV_DIR)
    if not os.path.isdir(FASTQ_DIR):
        sys.exit(f"Error: {FASTQ_DIR} could not be found")
    FASTQ_DIR = os.path.abspath(FASTQ_DIR)
    if args.sample_name:
        SAMPLE_NAME = args.sample_name
    if args.clip:
        CLIP = True


def import_fastq():
    global SAMPLE_NAME
    for fastq_file in FASTQ_FILES:
        # If we don't have a sample name we can guess it!
        if SAMPLE_NAME == "":
            # Unless it's the mux scan!
            if "mux_scan" in fastq_file:
                continue
            SAMPLE_NAME = '_'.join(fastq_file.split(".")[0].split("_")[2:]).replace(".fastq", "")
        if "mux_scan" in fastq_file:
            fastq_id = fastq_file + "_mux"
        else:
            fastq_id = fastq_file + "_seq"
        if fastq_id not in READ_SETS.keys():
            READ_SETS[fastq_id] = Read_Set(FASTQ_DIR, fastq_file)
            READ_SETS[fastq_id].read_fastq()


def add_csv_data_to_dataframes():
    for bin_number, read_set in READ_SETS.items():
        if not read_set.added_csv_data:
            read_set.append_csv_data()


def aggregate_dataframes():
    global ALL_READS
    for bin_number, read_set in READ_SETS.items():
        if read_set.aggregated_to_global_dataframe:
            continue
        if ALL_READS is None:
            columns = list(read_set.df.columns)
            ALL_READS = pd.DataFrame(columns=columns)
        ALL_READS = ALL_READS.append(read_set.df, ignore_index=True)
        read_set.aggregated_to_global_dataframe = True
    ALL_READS = ALL_READS.sort_values(['time'], ascending=[True])


def assign_yield_data():
    """We use this for both yield plots"""
    global YIELD_DATA
    # Read in seq length and time from ALL_READS dataframe
    YIELD_DATA = ALL_READS[['time', "seq_length"]]
    # Aggregate seq length for each minute of sequencing. I love this resample command!
    YIELD_DATA.set_index('time', inplace=True)
    YIELD_DATA = YIELD_DATA.resample("1T").sum()
    YIELD_DATA.reset_index(inplace=True)
    # Generate a cumulative sum of sequence data
    YIELD_DATA['cumsum_bp'] = YIELD_DATA['seq_length'].cumsum()
    # Convert time to timedelta format and then to float format, in seconds.
    YIELD_DATA['duration_tdelta'] = YIELD_DATA['time'].apply(lambda t: t - YIELD_DATA['time'].min())
    # Create a duration float.
    YIELD_DATA['duration_float'] = YIELD_DATA['duration_tdelta'].apply(lambda t: t.total_seconds())


def plot_yield_general():
    # Set subplots.
    fig, ax = plt.subplots(1)
    # Create ticks using numpy linspace. Ideally will create 6 points between 0 and 48 hours.
    num_points = 6
    x_ticks = np.linspace(YIELD_DATA['duration_float'].min(), YIELD_DATA['duration_float'].max(), num_points)
    ax.set_xticks(x_ticks)
    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels and limits.
    ax.set_xlabel("Duration (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_xlim(YIELD_DATA['duration_float'].min(), YIELD_DATA['duration_float'].max())
    ax.set_title(f"Yield for {SAMPLE_NAME} over run duration")
    ax.plot(YIELD_DATA['duration_float'], YIELD_DATA['cumsum_bp'],
            linestyle="solid", markevery=[])
    savefig(os.path.join(PLOTS_DIR, f"{SAMPLE_NAME}_yield_plot.png"))


def plot_yield_by_quality():
    # Read in seqlength and time from ALL_READS dataframe
    new_yield_data = ALL_READS[['time', "seq_length", "av_qual"]]
    # Bin qualities
    qual_bins = [0, 9, 15, 22, new_yield_data["av_qual"].max()]
    new_yield_data["descriptive_quality"] = pd.cut(new_yield_data["av_qual"], qual_bins,
                                                   labels=QUALITY_DESCRIPTIONS)
    new_yield_data.set_index(pd.DatetimeIndex(new_yield_data['time']), inplace=True)
    new_yield_data.drop('av_qual', axis=1, inplace=True)
    yield_data_grouped = new_yield_data.groupby("descriptive_quality").apply(lambda d: d.resample("1T").sum())[
        'seq_length']
    yield_data_by_quality = {description: yield_data_grouped[description].to_frame().reset_index()
                             for description in
                             QUALITY_DESCRIPTIONS}

    for description, yield_df in yield_data_by_quality.items():
        yield_df.reset_index(inplace=True)
        yield_df.set_index("time", inplace=True)
        yield_df = yield_df.reindex(index=YIELD_DATA.time, fill_value=0)
        yield_df.reset_index(inplace=True)
        # Generate a cumulative sum of sequence data
        yield_df['cumsum_bp'] = yield_df['seq_length'].cumsum()
        # Convert time to timedelta format and then to float format, in hours.
        yield_df['duration_tdelta'] = yield_df['time'].apply(lambda t: t - yield_df['time'].min())
        yield_df['duration_float'] = yield_df['duration_tdelta'].apply(lambda t: t.total_seconds() / 3600)
        yield_data_by_quality[description] = yield_df

    # Set subplots.
    fig, ax = plt.subplots(1)
    # Create ticks using numpy linspace. Ideally will create 6 points between 0 and 48 hours.
    num_points = 6
    x_ticks = np.linspace(YIELD_DATA['duration_float'].min(), YIELD_DATA['duration_float'].max(), num_points)
    ax.set_xticks(x_ticks)
    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels and limits.
    ax.set_xlabel("Duration (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_xlim(YIELD_DATA['duration_float'].min(), YIELD_DATA['duration_float'].max())
    ax.set_title(f"Yield for {SAMPLE_NAME}")
    ax.stackplot(YIELD_DATA['duration_float'],
                 [yield_data_by_quality[description]['cumsum_bp']
                  for description in reversed(QUALITY_DESCRIPTIONS)],
                 colors=reversed(QUALITY_COLOURS))
    # Convert bases to appropriate level
    formatter_y = FuncFormatter(y_yield_to_human_readable)

    # creating the legend manually
    ax.yaxis.set_major_formatter(formatter_y)
    ax.legend([mpatches.Patch(color=colour)
               for colour in reversed(QUALITY_COLOURS)],
              reversed(QUALITY_DESCRIPTIONS), loc=2)
    savefig(os.path.join(PLOTS_DIR, f"{SAMPLE_NAME}_yield_plot_by_quality.png"))


def plot_read_length_hist():
    num_bins = 50
    seq_df = ALL_READS["seq_length"]
    if CLIP:
        # Filter out the top 1000th percentile.
        seq_df = seq_df[seq_df < seq_df.quantile(0.9999)]

    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        num_bins = 50
        if y == 0:
            return 0
        s = humanize.naturalsize(seq_df.sum() * seq_df.count() * y / num_bins, gnu=True)
        return s + "b"

    # Define how many plots we want (1)
    fig, ax = plt.subplots(1)
    # Set the axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable_seq))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
    # Plot the histogram
    ax.hist(seq_df, num_bins, weights=seq_df,
            normed=1, facecolor='blue', alpha=0.76)
    # Set the titles and axis labels
    ax.set_title(f"Read Distribution Graph for {SAMPLE_NAME}")
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_xlabel("Read length")
    ax.set_ylabel("Bases per bin")
    savefig(os.path.join(PLOTS_DIR, f"{SAMPLE_NAME}_hist_read_length_by_basepair.png"))


def plot_poremap():
    def minknow_column_order(i):
        return chain(range(i + 33, i + 41), reversed(range(i + 1, i + 9)))
    # Channels are not in order, 121 is in the topleft, 89 in the top right.
    # The bottom left is 33 while the bottom right is channel 1.
    # The following four lines of code create an array that is a top-down, left-right 2D array of MinKNOW.

    # Split into chunks of 64 (rows of 4)
    chunks = [1, 2, 3, 4, 5, 6, 7, 0]
    # In which each row has the follow multiplication factor
    row_factors = [3, 2, 1, 0]
    # Create the values that make up the numbers on the far-right column of the grid.
    rh_values = [64 * chunk + 8 * row_factor for chunk in chunks for row_factor in row_factors]
    # Use the minknow_column_order function which reference the far-right column for a given row
    # to fill in the rest of the values for each row.
    channels_by_order_array = np.array([[j for j in minknow_column_order(i)] for i in rh_values])
    # Create an array of the same dimensions but filled with zeroes.
    channels_by_yield_array = np.zeros(channels_by_order_array.shape)
    # Sum the values for each channel.
    channels_by_yield_df = pd.DataFrame(ALL_READS.groupby("channel")["seq_length"].sum())
    # Reset the index and have channel as a column instead of the index.
    channels_by_yield_df.reset_index(level=0, inplace=True)
    # Iterate through each row of the yield by channel dataframe.
    for yield_row in channels_by_yield_df.itertuples():
        channel_index = [(ix, iy)
                         for ix, row in enumerate(channels_by_order_array)
                         for iy, i in enumerate(row)
                         if int(i) == int(yield_row.channel)][0]
        # Assign channel yield to position in MinKNOW
        channels_by_yield_array[channel_index] = yield_row.seq_length

    # The documentation for seaborn is pretty poor.
    # I will comment what I've done as best as possible.
    fig, ax = plt.subplots()

    fig.set_size_inches(15, 7)
    # Use the formatter we used for the yield plots.
    formatter_y = FuncFormatter(y_yield_to_human_readable)

    sns.heatmap(channels_by_yield_array,
                # Remove labels from side, they're not useful in this context.
                xticklabels=False,
                yticklabels=False,
                ax=ax,
                # Prevent extreme values from over-scaling the sidebar.
                robust=True,
                # Use the greens scale but in reverse, similar to MinKNOW.
                cmap="Greens_r",
                # Format keyword args for the side bar.
                cbar_kws={"format": formatter_y,
                          "label": "Bases per channel"})
    # Create line down the middle as shown in MinKNOW.
    ax.axvline([8], color='white', lw=15)
    # Nice big title!
    ax.set_title("Map of Yield by Channel", fontsize=25);
    savefig(os.path.join(PLOTS_DIR, f"{SAMPLE_NAME}_yield_map_by_pore.png"))


def plot_pore_yield_hist():
    num_bins = 50
    new_yield_data = ALL_READS.groupby(["channel", "mux"])['seq_length'].sum()
    fig, ax = plt.subplots(1)
    (n, bins, patches) = ax.hist(new_yield_data, num_bins, weights=None,
                                 # [1],#channels_by_yield_df['seq_length'],
                                 normed=1, facecolor='blue', alpha=0.76)
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))

    def y_muxhist_to_human_readable(y, position):
        # Get numbers of reads per bin in the histogram
        s = humanize.naturalsize((bins[1]-bins[0])*y*new_yield_data.count(), gnu=True)
        return s.replace("B", "")
    ax.yaxis.set_major_formatter(FuncFormatter(y_muxhist_to_human_readable))

    # Set the titles and axis labels
    ax.set_title(f"Yield by pore {SAMPLE_NAME}")
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_xlabel("Yield in single pore")
    ax.set_ylabel("Pores per bin")
    savefig(os.path.join(PLOTS_DIR, f"{SAMPLE_NAME}_hist_yield_by_pore.png"))


def y_hist_to_human_readable(y, position):
    # Convert distribution to base pairs
    num_bins = 50
    if y == 0:
        return 0
    s = humanize.naturalsize(ALL_READS["seq_length"].sum()*ALL_READS["seq_length"].count()*y/num_bins, gnu=True)

    return s + "b"


def x_hist_to_human_readable(x, position):
    # Convert distribution to base pairs
    if x == 0:
        return 0
    s = humanize.naturalsize(x, gnu=True)
    return s + "b"


def y_yield_to_human_readable(y, position):
    # Convert distribution to base pairs
    if y == 0:
        return 0
    s = humanize.naturalsize(y, gnu=True)

    return s + "b"


def x_yield_to_human_readable(x, position):
    # Convert time in seconds to hours or minutes
    hours = int(x/3600)
    minutes = int((x - 3600*hours)/60)
    if x == 0:
        return 0
    s = f"{hours:02d}:{minutes:02d}"
    return s


def is_basecalling():
    # Does the file basecalling exist (note must be in run directory).
    if os.path.isfile(os.path.join(os.getcwd(), "BASECALLING")):
        return True
    else:
        return False


def have_a_break():
    time.sleep(15)


def run_plot_commands():
    plot_yield_general()
    plot_yield_by_quality()
    plot_read_length_hist()
    plot_poremap()
    # CSV specific plots
    if not CSV_DIR == "":
        plot_pore_yield_hist()


def get_fastq_files():
    global FASTQ_FILES
    FASTQ_FILES = [fastq_file
                   for fastq_file in os.listdir(FASTQ_DIR)
                   if fastq_file.endswith(".fastq")]


def get_csv_files():
    global CSV_FILES
    CSV_FILES = [os.path.join(CSV_DIR, csv_file)
                 for csv_file in os.listdir(CSV_DIR)
                 if csv_file.endswith(".csv")]


def main(args):
    set_arguments(args)
    basecalling = True
    while basecalling:
        current_fastq_files = FASTQ_FILES.copy()
        get_fastq_files()
        if CSV_DIR is not "":
            get_csv_files()
        if len(list(set(FASTQ_FILES).difference(current_fastq_files))) == 0:
            have_a_break()
        import_fastq()
        if CSV_DIR is not "":
            add_csv_data_to_dataframes()
        aggregate_dataframes()
        assign_yield_data()
        run_plot_commands()
        if not is_basecalling():
            basecalling = False





