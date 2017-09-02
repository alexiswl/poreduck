#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
import argparse
import sys
import os
import statistics
from matplotlib import pyplot as plt
from matplotlib import dates as mdates
import numpy as np
import datetime
import matplotlib.patches as mpatches
import humanize
import matplotlib.mlab as mlab
from matplotlib.ticker import FuncFormatter
import matplotlib


"""
This script will create a yield plot of the data that has been created by the
albacore_server_scaled.py script and placed in the folder transfer_fast5_to_server.py
As it takes some time to handle and format the data such that the yield plots can be used...
we might as well make a bunch of plots to go along with it!
"""

# Before we begin, are we using python 3.6 or greater?
try:
    assert sys.version_info >= (3, 6)
    my_python_3p6_string = f"If you see a syntax error here you need to update your python version"
except (AssertionError, SyntaxError):
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
DATEPARSE = lambda dates: [pd.datetime.strptime(d, '%a %b %d %H:%M:%S %Y') for d in dates]
SAMPLE_NAME = ""
YIELD_DATA = None
QUALITY_DESCRIPTIONS = ["bad", "alright", "better", "best"]

# Import arguments.
"""
    csv directory
    fastq directory
    output directory
"""

class Read_set:
    """
    The class readset is assigned to each set of 4000 reads
    It contains attributes such as fastq file, csv file, random_number and mux boolean.
    """
    def __init__(self, csv_file):
        # Use fastq file to get meta data info
        self.csv_path = csv_file
        self.csv_file = os.path.basename(os.path.normpath(csv_file))
        self.number = self.csv_file.split("_")[0]
        self.random = self.csv_file.split("_")[1]
        if "mux" in self.csv_file:
            self.id = self.number + self.random + "_mux"
        else:
            self.id = self.number + self.random + "_seq"
        self.fastq_file = self.csv_file.replace(".csv", ".0.fastq")
        self.fastq_path = os.path.join(FASTQ_DIR, self.fastq_file)
        self.df = None
        self.added_fastq_data = False

    def read_csv(self):
        self.df = pd.read_csv(self.csv_path, header=0, parse_dates=['ctime'], date_parser=DATEPARSE)

    def append_fastq_data(self):
        if not os.path.isfile(self.fastq_path):
            print("Could not find path")
            return
        # Use the SeqIO package to iterate through the records:
        # Fastq header looks like this.
        # 'b1814d98-a01e-4fc6-a32c-60e1a062a957 runid=0608745933a900777aad6c8d9636f25227fe54a1 read=111140
        #  ch=351 start_time=2017-08-24T07:04:36Z
        # Create the columns we will write to fastq_id, and seq_length
        self.df["fastq_id"] = self.df.apply(lambda _: "DATA_NOT_FOUND", axis=1)
        self.df["seq_length"] = self.df.apply(lambda _: 0, axis=1)
        self.df["av_qual"] = self.df.apply(lambda _: 0, axis=1)
        for record in SeqIO.parse(self.fastq_path, "fastq"):
            fastq_id, run_id, read_no, channel, start_time = [description.split('=')[-1]
                                                              for description in record.description.split()]
            fastq_length = len(record.seq)
            fastq_av_quality = statistics.mean(record.letter_annotations["phred_quality"])
            df_index = self.df.query("channel==@channel & read_no==@read_no").index
            # Insert the fastq id
            self.df.set_value(df_index, "fastq_id", fastq_id)
            # Insert the fastq_file length
            self.df.set_value(df_index, "seq_length", fastq_length)
            # Insert the average quality
            self.df.set_value(df_index, "av_qual", fastq_av_quality)
        self.added_fastq_data = True

    def append_set_id_to_csv(self):
        self.df["read_set"] = self.df.apply(lambda _: self.id, axis=1)


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
    return(args)


def set_arguments(args):
    global CSV_DIR, FASTQ_DIR, PLOTS_DIR
    global CSV_FILES, FASTQ_FILES, SAMPLE_NAME
    CSV_DIR = args.csv_dir
    FASTQ_DIR = args.fastq_dir
    if args.output_dir:
        PLOTS_DIR = args.output_dir
        print(os.path.abspath(os.path.join(PLOTS_DIR, os.pardir)))
        if not os.path.isdir(os.path.abspath(os.path.join(PLOTS_DIR, os.pardir))):
            sys.exit(f"Error: Cannot create {PLOTS_DIR}. Parent directory does not exist.")
    else:
        PLOTS_DIR = os.path.join(CWD, "plots")
    if not os.path.isdir(PLOTS_DIR):
        os.mkdir(PLOTS_DIR)
    # Check and set CSV_DIR
    if not os.path.isdir(CSV_DIR):
        sys.exit(f"Error: {CSV_DIR} could not be found")
    CSV_DIR = os.path.abspath(CSV_DIR)
    CSV_FILES = [os.path.join(CSV_DIR, csv_file)
                 for csv_file in os.listdir(CSV_DIR)
                 if csv_file.endswith(".csv")]
    if not os.path.isdir(FASTQ_DIR):
        sys.exit(f"Error: {FASTQ_DIR} could not be found")
    FASTQ_DIR = os.path.abspath(FASTQ_DIR)
    FASTQ_FILES = [os.path.join(FASTQ_DIR, fastq_file)
                   for fastq_file in os.listdir(FASTQ_DIR)
                   if fastq_file.endswith(".fastq")]
    if args.sample_name:
        SAMPLE_NAME = args.sample_name


def import_csvs():
    global SAMPLE_NAME
    for csv_path in CSV_FILES:
        # If we don't have a sample name we can guess it!
        csv_file = os.path.basename(os.path.normpath(csv_path))
        if SAMPLE_NAME == "":
            if "mux_scan" in csv_path:
                continue
            SAMPLE_NAME = '_'.join(csv_file.split("_")[2:]).replace(".csv", "")
        csv_number = csv_file.split("_")[0]

        if "mux_scan" in csv_path:
            csv_id = csv_number + "_mux"
        else:
            csv_id = csv_number + "_seq"
        READ_SETS[csv_number] = Read_set(csv_path)
        READ_SETS[csv_number].read_csv()
        READ_SETS[csv_number].append_set_id_to_csv()



def add_fastq_data_to_dataframes():
    for bin_number, read_set in READ_SETS.items():
        read_set.append_fastq_data()


def aggregate_dataframes():
    global ALL_READS
    first_dataframe = True
    for bin_number, read_set in READ_SETS.items():
        if first_dataframe:
            first_dataframe = False
            columns = list(read_set.df.columns)
            ALL_READS = pd.DataFrame(columns=columns)
        ALL_READS = ALL_READS.append(read_set.df, ignore_index=True)
    ALL_READS = ALL_READS.sort_values(['ctime'], ascending=[True])


def assign_yield_data():
    """We use this for both yield plots"""
    global YIELD_DATA
    # Read in seqlength and time from ALL_READS dataframe
    YIELD_DATA = ALL_READS[['ctime', "seq_length"]]

    # Aggregate seqlength for each minute of sequencing. I love this resample command!
    YIELD_DATA.set_index('ctime', inplace=True)
    YIELD_DATA = YIELD_DATA.resample("1T").sum()
    YIELD_DATA.reset_index(inplace=True)
    # Generate a cumulative sum of sequence data
    YIELD_DATA['cumsum_bp'] = YIELD_DATA['seq_length'].cumsum()
    # Convert time to timedelta format and then to float format, in hours.
    YIELD_DATA['duration_tdelta'] = YIELD_DATA['ctime'].apply(lambda t: t - YIELD_DATA['ctime'].min())
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
    ax.set_title(f"Yield for {SAMPLE_NAME} (MB/Hour)")
    ax.plot(YIELD_DATA['duration_float'], YIELD_DATA['cumsum_mb'],
            linestyle="solid", markevery=[])
    plt.show()


def plot_yield_by_quality():
    # Read in seqlength and time from ALL_READS dataframe
    new_yield_data = ALL_READS[['ctime', "seq_length", "av_qual"]]
    # Bin qualities
    qual_bins = [0, 9, 15, 22, new_yield_data["av_qual"].max()]
    print(new_yield_data["av_qual"].max())
    new_yield_data["descriptive_quality"] = pd.cut(new_yield_data["av_qual"], qual_bins,
                                                   labels=QUALITY_DESCRIPTIONS)
    new_yield_data.set_index('ctime', inplace=True)
    new_yield_data.drop('av_qual', axis=1, inplace=True)
    yield_data_grouped = new_yield_data.groupby("descriptive_quality").apply(lambda d: d.resample("1T").sum())[
        'seq_length']
    yield_data_by_quality = {description: yield_data_grouped[description].to_frame().reset_index()
                             for description in
                             QUALITY_DESCRIPTIONS}

    for description, yield_df in yield_data_by_quality.items():
        yield_df.reset_index(inplace=True)
        yield_df.reindex(index=YIELD_DATA.index, fill_value=0)
        # Generate a cumulative sum of sequence data
        yield_df['cumsum_bp'] = yield_df['seq_length'].cumsum()
        # Convert time to timedelta format and then to float format, in hours.
        yield_df['duration_tdelta'] = yield_df['ctime'].apply(lambda t: t - yield_df['ctime'].min())
        yield_df['duration_float'] = yield_df['duration_tdelta'].apply(lambda t: t.total_seconds() / 3600)

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
    colours = ['red', 'orange', 'blue', 'green']

    ax.stackplot(YIELD_DATA['duration_float'],
                 [yield_data_by_quality[description]['cumsum_bp'] for description in QUALITY_DESCRIPTIONS],
                 colors=colours)
    # Convert bases to appropriate level
    formatter_y = FuncFormatter(y_yield_to_human_readable)

    # creating the legend manually
    ax.yaxis.set_major_formatter(formatter_y)
    ax.legend([mpatches.Patch(color=colour) for colour in reversed(colours)],
              reversed(QUALITY_DESCRIPTIONS), loc=2)

    plt.show()


def plot_read_length_hist():
    # Define how many plots we want (1)
    fig, ax = plt.subplots(1)
    # Set the axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
    # Plot the histogram
    ax.hist(ALL_READS["seq_length"], 50, weights=ALL_READS["seq_length"],
            normed=1, facecolor='blue', alpha=0.76)
    # Set the titles and axis labels
    ax.set_title(f"Read Distribution Graph for {SAMPLE_NAME}")
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_xlabel("Read length")
    ax.set_ylabel("Bases per bin")

    plt.show()


def y_hist_to_human_readable(y, position):
    # Convert distribution to mega basepairs
    if y == 0:
        return 0
    s = humanize.naturalsize(ALL_READS["seq_length"].sum()*ALL_READS["seq_length"].count()*y, gnu=True)

    return s + "b"


def x_hist_to_human_readable(x, position):
    # Convert distribution to mega basepairs
    if x == 0:
        return 0
    s = humanize.naturalsize(x, gnu=True)
    return s + "b"


def y_yield_to_human_readable(y, position):
    # Convert distribution to mega basepairs
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


def main():
    args = get_arguments()
    set_arguments(args)
    import_csvs()
    add_fastq_data_to_dataframes()
    aggregate_dataframes()
    assign_yield_data()
    plot_yield_general()
    plot_yield_by_quality()
    plot_read_length_hist()


main()


# Get csv's
# Required to generate the yield plot.
# Each csv is written to a pandas dataframe.


# Use SeqIO to map reads to length
# Add this to end of dataframe.

# Concatenate csvs

# Sort by date column.

# Group by minute, aggregate read-lengths

# Aggregate yield by minute

# Plot yield using matplot lib

# Plot quality using matplotlib