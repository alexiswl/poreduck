#!/usr/bin/env python3

"""
Given a sample sheet and a run directory, tar up a PromethION run.
"""

import argparse
import h5py
from datetime import datetime, timedelta
import pandas as pd
import os
import shutil
import subprocess
import time
# Use the python tar module to pipe data into gzip.
import tarfile
import sys
# Import matplotlib and friends
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import humanfriendly
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
import seaborn as sns
import numpy as np
from itertools import chain
import re

"""
Class types

Run:
A run exists for sample name.
Each run has a subset of the sample sheet dataframes.
We imply the read paths from each flowcell being used.

Subfolder:
Folder specified as an int value.
By default each folder holds 4000 fast5 files.
Sets up a dataframe for each fastq file


Fastq File:
Grabs meta data about each fastq file.
Along with the expected finish time of the run.

Issues:

Testing:

Solved:

To do:

"""


class Fast5file:
    def __init__(self, filename, input_folder, is_mux=False):
        self.filename = filename
        self.file_path = os.path.join(input_folder, self.filename)
        # To get the rest of the attributes from the filename,
        # Split the filename into that which is before 'sequencing_run'
        # And that which is after sequencing run
        if is_mux:
            pre_seq_pivot, post_seq_pivot = self.filename.rsplit("mux_scan", 1)
        else:
            pre_seq_pivot, post_seq_pivot = self.filename.rsplit("sequencing_run", 1)
        pre_seq_pivot = pre_seq_pivot.split("_")
        post_seq_pivot = post_seq_pivot.split("_")
        # Slurm node, date and f_id are all premux
        # Not needed as of yet
        # Channel, read, rnumber and sample_id are all post pivot
        self.channel = post_seq_pivot[-2]
        self.read = post_seq_pivot[-4]
        self.rnumber = post_seq_pivot[-6]
        self.sample_id = post_seq_pivot[0:-6]
        self.corrupted = False
        # Now get inside the fast5 file
        with h5py.File(self.file_path) as f:
            # Get attributes from /Raw/Reads/
            try:
                read_attributes = dict(f['Raw/Reads/Read_%s' % self.read].attrs.items())
            except KeyError:
                self.corrupted = True
                print("%s is corrupted" % self.filename)
                return
            # Get values from inside the fast5 value
            self.mux_id = read_attributes["start_mux"]
            self.read_id = read_attributes["read_id"]
            # Get the Experiment duration from the context tags.
            try:
                context_tags = dict(f['UniqueGlobalKey/context_tags'].attrs.items())
                tracking_id = dict(f['UniqueGlobalKey/tracking_id'].attrs.items())
                channel_id = dict(f['UniqueGlobalKey/channel_id'].attrs.items())
            except KeyError:
                self.corrupted = True
                print("%s is corrupted" % self.filename)
                return
            # Time (even mux have start times)
            self.exp_start_time = datetime.strptime(tracking_id['exp_start_time'].decode(), 
                                                    "%Y-%m-%dT%H:%M:%SZ") 
            # Minute format and in byte strings (not found in mux)
            if not is_mux:
                self.exp_duration_set = timedelta(minutes=int(context_tags['experiment_duration_set'].decode()))
            else:
                self.exp_duration_set = timedelta(minutes=10)
            # Read start time = start_time / sampling rate + exp_start_time
            self.duration_time = int(read_attributes["duration"])/int(channel_id["sampling_rate"])
            start = int(read_attributes["start_time"])
            sampling_rate = int(channel_id["sampling_rate"])
            self.pore_start_time = timedelta(seconds=start/sampling_rate) + self.exp_start_time
            self.pore_end_time = self.pore_start_time + timedelta(seconds=self.duration_time)

    def to_series(self):
        series_columns = ["Name", "Channel", "Read", "RNumber", # From file name
                          "MuxID", "StartTime", "EndTime"]  # From fast5 file - requried for plots

        return pd.Series(data=[self.filename, self.channel, self.read, self.rnumber,
                               self.mux_id, self.pore_start_time, self.pore_end_time],
                         index=series_columns)


class Subfolder:
    def __init__(self, reads_path, number, metadata_dir, run, is_mux=False, threshold=4000):
        # Int name of the fast5 file
        self.number = number
        self.standard_int = self.number.zfill(4)
        self.is_mux = is_mux
        self.pardir = reads_path
        self.path = os.path.join(self.pardir, number) 
        self.metadata_dir = metadata_dir
        self.metadata_path = ""
        # Initialise the stage parameters. 
        self.is_full = False
        self.is_tarred = False
        # Initialise fast5_file list and dataframe
        self.fast5_files = []
        self.pd = None
        # Initialise start and end times
        self.rnumber = ""
        self.start_time = None
        self.end_time = None
        self.threshold = threshold
        self.new_folder_name = ""
        self.new_folder_path = ""
        self.tar_file = ""
        self.tar_path = ""
        self.num_fast5_files = 0
        self.run = run
        self.md5sum = None

    def get_new_folder_name(self): 
        # Create the new folder name
        if self.is_mux:
            mux_seq = "mux_scan"
        else:
            mux_seq = "sequencing_run"
        self.new_folder_name = '_'.join([self.standard_int, mux_seq, self.rnumber])
        self.new_folder_path = os.path.join(self.pardir, self.new_folder_name)
        self.tar_file = self.new_folder_name + ".fast5.tar.gz"
        self.tar_path = os.path.join(self.pardir, self.tar_file)
        self.metadata_path = os.path.join(self.metadata_dir, self.new_folder_name+".tsv")

    def get_fast5_files(self):
        # Fast5 class
        self.fast5_files = [Fast5file(fast5_file, self.path, is_mux=self.is_mux)
                            for fast5_file in os.listdir(path=self.path)
                            if fast5_file.endswith(".fast5")]
        self.num_fast5_files = len(self.fast5_files)
        if self.num_fast5_files == 0:
            # We get in here when empty folders exist post run.
            self.is_full = False
    
    def get_dataframe(self):
        # Generate dataframe for subfolder
        self.pd = None
        for fast5_file in self.fast5_files:
            if fast5_file.corrupted:
                continue
            if self.pd is None:
                first_series = fast5_file.to_series()
                
                self.pd = pd.DataFrame(data=[first_series])

            else:
                self.pd = self.pd.append(fast5_file.to_series(),
                                         ignore_index=True)

    def write_dataframe(self):
        self.pd.to_csv(self.metadata_path, header=True, index=False, sep="\t")

    def folder_exists(self):
        if os.path.isdir(self.path):
            return True
        else:
            return False

    def is_empty(self):
        # Are there any files in this directory at all?
        all_files_count = len([any_file
                               for any_file in os.listdir(path=self.path)])
        if all_files_count == 0:
            return True
        else:
            return False

    def check_if_full(self):    
        if self.is_full:
            # We shouldn't be calling this twice.
            return
        # Check the current status of the directory.
        raw_fast5_count = len([fast5
                               for fast5 in os.listdir(path=self.path)
                               if fast5.endswith(".fast5")])

        # Determine if this folder needs deleting
        if self.is_empty():
            return

        # Determine if this folder is still being written to
        if raw_fast5_count < self.threshold and not self.run.complete:
            if self.is_mux:
                self.is_full = True
            else:
                self.is_full = False
                return
        self.is_full = True
        # Get fast5 files
        self.get_fast5_files()
        # Get dataframe
        self.get_dataframe()
        # Redetermine if bin is full
        if not self.is_full and self.run.complete:
            # Get fast5 files
            self.get_fast5_files()
            # Get data frame
            self.get_dataframe()

        # Set values of rnumber, start_time and end_time
        try:
            self.rnumber = self.pd.RNumber.unique().item()
        except ValueError:
            print(self.pd.RNumber.unique())
        self.start_time = self.pd.StartTime.min()
        self.end_time = self.pd.EndTime.max()
        # Get new folder name
        self.get_new_folder_name()
        self.write_dataframe()

    def tar_folder(self):
        # Always ensure the previous stage has been completed
        if not self.is_full and not self.run.complete:
            return
        # Always ensure the this stage has been
        # First move the folder to the new folder path location
        os.rename(self.path, self.new_folder_path)

        """Tar folder using pigz"""
        with tarfile.open(name=self.tar_path+".tmp", mode="w:gz") as tar_h:
            tar_h.add(name=self.new_folder_path,  # Set arc name as relative path
                      arcname=os.path.basename(os.path.normpath(self.new_folder_path)),
                      recursive=True)
        # Just let the filesystem catch up
        time.sleep(1)
        # Now delete the path that exists
        shutil.rmtree(self.new_folder_path)

        # On the rare occasion, the system may not be able to find this file even if it exists.
        # Generate a while loop and exit after fifteen seconds if no file is found.
        index = 0
        while True:
            if index > 3:
                sys.exit("Error, could not found %s.tmp" % self.tar_path)
            try:
                os.rename(self.tar_path+".tmp", self.tar_path)
                break
            except IOError:
                print("Warning: %s.tmp not found." % self.tar_path)
                print("Attempting to find again in fifteen seconds")
                time.sleep(15)
                index += 1
        self.is_tarred = True

    def get_tar_md5(self):
        md5_command = ["md5sum", self.tar_path]
        md5_proc = subprocess.run(md5_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if not md5_proc.returncode == 0:
            print("Warning, md5_proc returned error code %s" % md5_proc.returncode)
        stdout = md5_proc.stdout.decode()
        self.md5sum = stdout.split(" ", 1)[0].strip()


class Run:
    def __init__(self, path, name, start_date, start_time, is_mux=False):
        self.path = path
        self.name = name
        self.fast5_path = os.path.join(self.path, "fast5")
        self.is_mux = is_mux
        self.complete = False
        self.start_time = start_time
        self.start_date = start_date
        self.completion_time = None
        self.subfolders = []
        self.metadata_dir = os.path.join(self.path, "metadata")
        self.plots_dir = os.path.join(self.path, "plots")
        self.checksum = os.path.join(self.path, "checksum.md5")
        self.df = None
        if not os.path.isdir(self.metadata_dir):
            os.mkdir(self.metadata_dir)
        if not os.path.isdir(self.plots_dir):
            os.mkdir(self.plots_dir)

    def get_subfolders(self):
        # Get all folders currently in the directory
        folders = sorted([folder
                          for folder in os.listdir(path=self.fast5_path)
                          if folder.isdigit()],
                          key=lambda x: int(x))
        for folder in folders:
            # Don't read in current folders
            if folder in [subfolder.number
                          for subfolder in self.subfolders]:
                continue
            # Append new folders
            self.subfolders.append(Subfolder(self.fast5_path, folder, self.metadata_dir, self, is_mux=self.is_mux))

    def tar_subfolders(self):
        for folder in self.subfolders:
            if folder.is_tarred:
                continue # Only the first folder will get to here.
            if not folder.is_full:
                folder.check_if_full()
            if folder.is_full and not folder.is_tarred:
                folder.tar_folder()
                folder.get_tar_md5()

    def slim_tarred_subfolders(self):
        # For each subfolder, unlink the list of fast5 files.
        # Keep only the pandas dataframe.
        for subfolder in self.subfolders:
            if subfolder == self.subfolders[0]:
                continue
            elif not subfolder.is_tarred:
                # In case we need to reference the time again.
                continue
            else:
                try:
                    del subfolder.fast5_files
                except AttributeError:
                    # Fast5 files already deleted
                    pass
                    
    def get_run_finish_time(self):
        # Get standard fast5 file (not that simple)
        if len(self.subfolders) == 0:
            print("No subfolders, using folder names to get run finish time")
            return self.get_default_finishtime()
        for subfolder in self.subfolders:
            fast5_files_iter = iter(subfolder.fast5_files)
            fast5_file = None
            while fast5_file is None or fast5_file.corrupted:
                try:
                    fast5_file = next(fast5_files_iter)
                except StopIteration:
                    break
            if fast5_file is not None and not fast5_file.corrupted:
                self.start_time = fast5_file.exp_start_time
                minutes = fast5_file.exp_duration_set
                return self.start_time + minutes
            # If we are here, it means no folder is full. We expect run is complete
            return self.get_default_finishtime()

    def get_default_finishtime(self):
        # Use the start_date and start_time to get the proposed finish times
        print("Getting default run finish times")
        start_string = '_'.join([self.start_date, self.start_time])
        start_date_ob = datetime.strptime(start_string, "%Y%m%d_%H%M")
        if self.is_mux:
            end_date_ob = start_date_ob + timedelta(minutes=8)
        else:
            end_date_ob = start_date_ob + timedelta(hours=48)
        return end_date_ob

    def is_run_complete(self):
        # Get expected finish time from fast5 file
        if self.completion_time is None:
            self.completion_time = self.get_run_finish_time()
        # Determine if run is complete.
        current_time = datetime.utcnow()
        # If difference is less than zero, run is finished
        diff = self.completion_time - current_time 
        if diff.total_seconds() < 0:
            self.complete = True
            return True
        else:
            return False

    def write_md5(self):
        """
        Append to the checksum file. 0000_SAMPLE_XYZ.tar.gz  md5sumhashcode
        Rsync needs to not include this file in transfer.
        """
        with open(self.checksum, 'w') as check_h:
            for subfolder in self.subfolders:
                if subfolder.md5sum is not None:
                    check_h.write("%s  fast5/%s\n" % (subfolder.md5sum, subfolder.tar_file))

    def get_bulk_metadata(self):
        """
        Merge all of the subfolder data frames into one for plotting purposes.
        We will delete this to save memory after each plot iteration
        A plot shouldn't take more than 2 minutes to be created.
        """
        self.df = pd.concat([subfolder.pd for subfolder in self.subfolders
                             if subfolder.pd is not None],
                            axis=0, ignore_index=True)
        self.df['RunDurationTime'] = self.df["EndTime"].apply(lambda x: x - self.start_time)
        self.df['RunDurationFloat'] = self.df["RunDurationTime"].apply(lambda x: x.total_seconds())
        self.df['EstLength'] = self.df[["StartTime", "EndTime"]].apply(lambda x: self.estimate_read_length(x), axis=1)

    def estimate_read_length(self, row):
        """
        Takes in a dataframe of starttime and endtime
        Retuns a series using the same index of the estimated length of a given read
        """
        speed = 450  # Bases per second
        # Return a timedelta object converted into a float
        return speed * (row.EndTime - row.StartTime).seconds

    def plot_yield(self):
        """
        Plot the estimated yield based on the metadata
        """
        plt.close('all')
        fig, ax = plt.subplots(1, figsize=(10, 10))
        # remove any unforeseen NA values
        self.df.dropna(inplace=True, how='any')
        self.df.sort_values(by="RunDurationTime", inplace=True)
        self.df.set_index("RunDurationFloat", inplace=True)
        # Plot cumulative yield over time
        self.df['EstLength'].cumsum().plot(ax=ax)
        # Define axis formatters
        ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
        ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
        # Set x and y labels
        ax.set_xlabel("Duration of run (HH:MM)")
        ax.set_ylabel("Yield")
        ax.set_title("Theoretical Yield for sample %s over time" % self.name)
        # Set x and y limits
        ax.set_ylim(ymin=0)
        # Add legend
        ax.legend([self.name])
        # Ensure labels are not missed'
        fig.tight_layout()
        savefig(os.path.join(self.plots_dir, "%s.theoretical_yield.png" % self.name))
        # Reset index for next plots
        self.df.reset_index(inplace=True)

    def plot_hist(self):
        plt.close('all')
        fig, ax = plt.subplots(1, figsize=(10, 10))
        num_bins = 50
        # Trim histogram,
        trimmed = self.df.query("EstLength < %d" % self.df.EstLength.quantile(0.995))["EstLength"]
        bins = np.linspace(start=0, stop=trimmed.max(), num=num_bins)
        bin_width = bins[1] - bins[0] 
        # Plot histogram
        trimmed.plot(kind="hist", ax=ax, density=1, bins=bins, alpha=0.6, weights=trimmed)

        # Local definition of formatter, requires input of bin_width and entire data frame
        def y_hist_to_human_readable(y, position):
            if y == 0:
                return 0
            s = humanfriendly.format_size(bin_width * trimmed.sum() * y, binary=False)
            return reformat_human_friendly(s)
        # Now use this format to show the base pairs per bin
        ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable))
        ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
        # Set titles
        ax.set_title("Read Distribution Graph for %s" % self.name)
        ax.grid(color='black', linestyle=':', linewidth=0.5)
        ax.set_ylabel("Bases per bin")
        ax.set_xlabel("Read length")
        # Create legend
        ax.legend([self.name])
        # Ensure labels are not missed
        fig.tight_layout()
        savefig(os.path.join(self.plots_dir, "%s.theoretical_hist.png" % self.name))
        del trimmed

    def plot_flowcell(self):
        poremap_df = get_poremap_from_yield_df(self.df)
        plt.close('all')
        # Use the df by channel to plot the flowcell, we want this one to be longer than it is wide
        fig, ax = plt.subplots(1, figsize=(20, 10))
        # Use the formatter we used for the yield plots.
        formatter_y = FuncFormatter(y_yield_to_human_readable)
        sns.heatmap(poremap_df,
                    # Remove labels from side, not appropriate
                    xticklabels=False,
                    yticklabels=False,
                    ax=ax,
                    # Prevent extreme values from over-scaling the bar
                    robust=True,
                    # Use the greens scale but in reverse, similar to MinKNOW
                    cmap="Greens_r",
                    # Format the keyword args for the side bar.
                    cbar_kws={"format": formatter_y,
                              "label": "Bases per channel"})
        # Create a small line down the middle of the graph as shown in MinKNOW
        ax.axvline([8], color='white', lw=15)
        # Add a nice big title!
        ax.set_title("Map of Yield by Channel", fontsize=25)
        # Ensure labels are not missed
        fig.tight_layout()
        savefig(os.path.join(self.plots_dir, "%s.theoretical_yield_map_by_pore.png" % self.name))
        del poremap_df

    def print_theoretical_stats(self):
        """
        Print the total yield, nx values and the run duration
        """
        percentiles = [0.1, 0.25, 0.5, 0.75, 0.9]

        # Get total yield metrics
        total_bp = self.df.EstLength.sum()
        total_bp_h = reformat_human_friendly(humanfriendly.format_size(total_bp, binary=False))
        # Describe the seq length histogram:
        total_bp_describe = self.df.EstLength.describe(percentiles=percentiles).to_string()
        total_bp_describe = '\n'.join([line.split()[0].ljust(8) + "\t" +
                                       "{:21.2f}".format(float(line.split()[1]))
                                       for line in total_bp_describe.split("\n")])

        # Calculate the NX of the read lengths where X is 0.1, 0.25, 0.5, 0.75 and 0.9
        nx = []
        seq_length_sorted_as_series = self.df['EstLength'].sort_values().reset_index(drop=True)
        seq_length_cumsum_as_series = seq_length_sorted_as_series.cumsum()
        # Iterate through series, add value to nx_list if crosses a cumsum percentile
        for index, seq_value in seq_length_sorted_as_series.iteritems():
            if (seq_length_cumsum_as_series[index] <= total_bp*percentiles[len(nx)] <=
              seq_length_cumsum_as_series[index+1]):
                nx.append(seq_value)
            if len(nx) == len(percentiles):
                # Found all the percentiles, no need to continue
                break
        nx_h = [reformat_human_friendly(humanfriendly.format_size(nX_value, binary=False))
                for nX_value in nx]
        # Print these stats to file
        with open(os.path.join(self.plots_dir, "%s.stats.txt" % self.name), 'w') as output_handle:
            # Print total basepairs
            output_handle.write("# Stats for sample '%s' #\n" % self.name)
            output_handle.write("Total basepairs:\n")
            output_handle.write(f"\t{total_bp:16,.0f}\t|\t{total_bp_h.rjust(9)}\n")
            output_handle.write("Description of Read Lengths:\n")
            # Tab indent each of the descriptor lines
            output_handle.writelines(f"\t{len_line}\n"
                                     for len_line in total_bp_describe.split("\n"))
            # Write output for each of the nX values
            output_handle.write("NX values:\n")
            output_handle.writelines(f"\tN{100*percentile:02.0f}:\t{nx_value:8,.0f}\t|\t{nx_h_value.rjust(9)}\n"
                                     for percentile, nx_value, nx_h_value in zip(percentiles, nx, nx_h))
            duration = self.df["RunDurationFloat"].max()  # In seconds
            hours, remainder = divmod(duration, 3600)
            minutes, seconds = divmod(remainder, 60)
            run_duration_h = f"{hours} hours, {minutes} minutes, {seconds:2,.0f} seconds"
            output_handle.write("Run Duration Time:\n")
            output_handle.write(f"\t{duration:8,.1f} seconds\t|\t{run_duration_h}\n")


class Sample:
    def __init__(self, sample_name, samplesheet, reads_path):
        self.pd = samplesheet.query("SampleName=='%s'" % sample_name)
        # Get the active runs for this sample
        self.runs = []
        self.is_running = True
        for index, run in self.pd.iterrows():
            if run.SampleName.startswith("%"):  # Sample name not specified in MinKNOW
                mux_path = os.path.join(reads_path, '_'.join([run.UTCMuxStartDate, run.UTCMuxStartTime]))
                seq_path = os.path.join(reads_path, '_'.join([run.UTCSeqStartDate, run.UTCSeqStartTime]))
                run.SampleName = re.sub("^%", "", run.SampleName)
            else:
                mux_path = os.path.join(reads_path, '_'.join([run.UTCMuxStartDate, run.UTCMuxStartTime, run.SampleName]))
                seq_path = os.path.join(reads_path, '_'.join([run.UTCSeqStartDate, run.UTCSeqStartTime, run.SampleName]))
            self.runs.append(Run(mux_path, run.SampleName, run.UTCMuxStartDate, run.UTCMuxStartTime, is_mux=True))
            self.runs.append(Run(seq_path, run.SampleName, run.UTCSeqStartDate, run.UTCSeqStartTime, is_mux=False))
    
    def is_run_complete(self):
        # All samples must be complete to return true.
        for run in self.runs:
            if not run.is_run_complete(): 
                return False
        return True


# Plotting scaling formatters
def reformat_human_friendly(s):
    """
    humanfriendly module returns with a few quirks
    1 = 1 byte ==> 1 b
    2 = 2 bytes ==> 2 b
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


def get_poremap_from_yield_df(run_df):
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
    channels_by_yield_df = pd.DataFrame(run_df.groupby("Channel")["EstLength"].sum())
    # Reset the index and have channel as a column instead of the index.
    channels_by_yield_df.reset_index(level=0, inplace=True)
    # Iterate through each row of the yield by channel dataframe.
    for yield_row in channels_by_yield_df.itertuples():
        channel_index = [(ix, iy)
                         for ix, row in enumerate(channels_by_order_array)
                         for iy, i in enumerate(row)
                         if int(i) == int(yield_row.Channel)][0]
        # Assign channel yield to position in MinKNOW
        channels_by_yield_array[channel_index] = yield_row.EstLength
    return channels_by_yield_array


"""
General process:
1. Get runs.
2. Each run has a boolean attribute - running.
3. While each run is running, continue to tar up subfolders but not those that don't have 4000 reads.
"""


def get_args():
    """
    Two simple arguments.
    1. Path to PCA directory
    2. Samplesheet
    """
    parser = argparse.ArgumentParser(description="Tar up the run folders of the PromethION")
    parser.add_argument("--samplesheet", type=str, required=True,
                        help="Path to tab delimited samplesheet. "
                             "Columns are SampleName, UTCMuxStartDate, UTCMuxStartTime, "
                             "UTCSeqStartDate, UTCSeqStartTime")
    parser.add_argument("--reads_path", type=str, required=True,
                        help="path/to/reads. "
                             "Ubuntu: /var/lib/MinKNOW/data/reads"
                             "Mac: /Library/MinKNOW/data"
                             "Windows: C:\\data\\reads")
    args = parser.parse_args()
    return args                                
               
                                     
def is_still_running(samples):
    for sample in samples:
        if not sample.is_run_complete():
            return True
    return False
            

def samplesheet_to_pd(samplesheet):
    return pd.read_csv(samplesheet, header=0, sep="\t", dtype=str, comment='#')


def main(args):
    samplesheet = samplesheet_to_pd(args.samplesheet) 
    samples = [Sample(sample, samplesheet, args.reads_path)
               for sample in samplesheet.SampleName.unique().tolist()]
    running = True
    first_pass = True
    while running:
        if not first_pass:
            running = is_still_running(samples)
            for sample in samples:
                for run in sample.runs:
                    time.sleep(15)
                    # This is going to break if no subfolders
                    if len([subfolder.pd for subfolder in run.subfolders
                            if subfolder.pd is not None]) == 0:
                        continue
                    run.slim_tarred_subfolders()
                    run.get_bulk_metadata()
                    run.plot_yield()
                    run.plot_hist()
                    run.plot_flowcell()
                    run.print_theoretical_stats()
        else:
            first_pass = False
        for sample in samples:
            for run in sample.runs: 
                run.get_subfolders()
                run.tar_subfolders()
                run.write_md5()


if __name__ == "__main__":
    main()
