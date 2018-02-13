#!/usr/bin/env python3

from promethion_starter import samplesheet_to_pd, config_to_pd
import pandas as pd
import os
import argparse
import paramiko
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import humanfriendly
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig


series_columns = ["Name", "Channel", "Read", "RNumber", # From file name
                  "MuxID", "StartTime", "EndTime"]  # From fast5 file - requried for plots


class Sample:
    def __init__(self, sample_name, sample_df):
        self.name = sample_name
        self.df = sample_df
        self.runs = []
        self.concatenated_data = None

class Slave:
    """Use the slave node config to access the data from the master node."""
    def __init__(self, slurm_id, ip):
        self.slurm_id = slurm_id
        self.ssh_ip = ip
        self.ssh_client = None
        self.reads_path = "/tmp/output/reads"

    def connect(self):
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh_client.connect(self.ssh_ip, key_filename='/home/prom/.ssh/id_rsa.pub')


class Run:
    def __init__(self, sample_name,
                 grnwch_mux_start_date, grnch_mux_start_time,
                 grnwch_seq_start_date, grnch_seq_start_time,
                 slurm_id, flowcell_id, ip):
        self.sample_name = sample_name
        self.grnwch_mux_start_date = grnwch_mux_start_date
        self.grnch_mux_start_time = grnch_mux_start_time
        self.grnwch_seq_start_date = grnwch_seq_start_date
        self.grnch_seq_start_time = grnch_seq_start_time
        self.slurm_id = slurm_id
        self.flowcell_id = flowcell_id
        self.concatenated_data = None
        self.slave = Slave(self.slurm_id, ip)
        self.path = os.path.join("/media/data", '_'.join([self.grnwch_seq_start_date,
                                                          self.grnwch_seq_start_time,
                                                          self.sample_name]))

def get_args():
    """Get arguments from commandline"""
    """
    Two simple arguments.
    1. Path to PCA directory
    2. Samplesheet
    """
    parser = argparse.ArgumentParser(description="Plot a run as it is going")
    parser.add_argument("--samplesheet", type=str, required=True,
                        help="Path to tab delimited samplesheet. "
                             "Columns are SampleName, GrnwchMuxStartDate, GrnwchMuxStartTime, "
                             "GrnwchSeqStartDate, GrnwchSeqStartTime SlurmID FlowcellID")
    parser.add_argument("--ip_config", type=str, required=True,
                        help="path/to/tab-delimited-config file. "
                             "One column of IP addresses ==> one column of slave nodes.")
    args = parser.parse_args()
    return args


def get_samples(samplesheet_df):
    # Get samples
    samples = []
    for sample_name in samplesheet_df["SampleName"].unique().tolist():
        samples.append(Sample(sample_name, samplesheet_df.query("SampleName==@sample_name")))


def get_runs(sample_df, ip_config):
    """Return a list of class run, each with a list of associated dataframes"""
    runs = []
    for index, row in sample_df.iterrows():
        ip = ip_config.query("SlurmID==@row.SlurmID")["IP"].item()
        runs.append(Run(row.sample_name,
                        row.GrnwchMuxStartDate, row.GrnwchMuxStartTime,
                        row.GrnwchSeqStartDate, row.GrnwchSeqStartTime,
                        row.SlurmID, row.FlowcellID, ip))
    return runs


def get_metadata(path, open_sftp):
    """Return a list of pandas dataframes imported over paramiko"""
    metadata_path = os.path.join(path, "metadata")
    # Get list of tsv files
    tsv_list = [tsv
                for tsv in open_sftp.listdir(path=metadata_path)
                if tsv.endswith(".tsv")]
    # Read in dataframe via read_csv()
    metadatas = []
    for tsv in tsv_list:
        with open_sftp.file(os.path.join(metadata_path, tsv), 'r') as tsv_ftpfile:
            metadata = pd.read_csv(tsv_ftpfile, header=True, index=False, sep="\t",
                                   parse_dates=["StartTime", "EndTime"])
        # Convert the following columns to numeric
        metadata[["Channel", "Read", "Mux", "RNumber"]] = metadata[["Channel", "Read", "Mux", "RNumber"]].apply(pd.to_numeric)
        metadata["EstLength"] = metadata[["StartTime", "EndTime"]].apply(estimate_read_length)
        # Add in the duration time column
        metadata["RunDurationTime"] = metadata["StartTime"] - metadata["StartTime"].min()
        metadatas.append(metadata)
    # Return a merged version of the metadata
    return merge_metadatas(metadatas)


def merge_metadatas(metadatas):
    """Merge a list of dataframes for a given run"""
    return pd.concat(metadatas, ignore_index=True)


def plot_samples(names, samples_df):
    """Plot an estimated yield plot and a histogram plot comparing each flowcell"""
    # Need to configure the index properly for this one and incorporate just the required dataframes
    # Yield plot
    # Set up plotting structure
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Set and sort the index by time
    samples_df.set_index("RunDurationTime", inplace=True)
    samples_df.sort_index(inplace=True)
    # Generate a cumulative yield
    sample_df["CumYield"] = sample_df["EstLength"].cumsum()
    # Plot cumulative yield over time.
    sample_df["CumYield"].plot(ax=ax)
    # Define axis formatters
    ax.xaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM")
    ax.set_ylabel("Yield")
    ax.set_title("Yield for sample %s over time" % name)
    # Set x and y limits
    ax.set_xlim(sample_df["DurationTime"].min(), sample_df["DurationTime"].max())
    ax.set_ylim(ymin=0)
    # Ensure labels are not mised
    fig.tight_layout()
    savefig("%s.combined_yield.png" % name)


def plot_single(name, sample_df):
    """Plot an estimated yield plot and a histogram plot for the sample"""
    # Yield plot
    # Set up plotting structure
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Set and sort the index by time
    sample_df.set_index("RunDurationTime", inplace=True)
    sample_df.sort_index(inplace=True)
    # Generate a cumulative yield
    sample_df["CumYield"] = sample_df["EstYield"].cumsum()
    # Plot cumulative yield over time.
    sample_df["CumYield"].plot(ax=ax)
    # Define axis formatters
    ax.xaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM")
    ax.set_ylabel("Yield")
    ax.set_title("Yield for sample %s over time" % name)
    # Set x and y limits
    ax.set_xlim(sample_df["DurationTime"].min(), sample_df["DurationTime"].max())
    ax.set_ylim(ymin=0)
    # Ensure labels are not mised
    fig.tight_layout()
    savefig("%s.combined_yield.png" % name)

    # Histogram
    # Set up plotting structure
    plt.close('all')
    fig, ax = plt.subplots(1)
    num_bins = 50
    # Clip data and plot histogram
    lengths = sample_df['EstLength']
    lengths = lengths[lengths < lengths.quantile(0.995)]
    lengths.plot(type='hist', ax=ax, normed=1, facecolor='blue', alpha=0.75)
    # Set the axis formatters
    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        if y == 0:
            return 0
        s = humanfriendly.format_size(lengths.sum() * y, binary=False)
        return reformat_human_friendly(s)
    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable_seq))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
    # Set titles
    ax.set_title("Read Distribution Graph for %s" % name)
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_ylabel("Bases per bin")
    # Ensure labels are not mised.
    fig.tight_layout()
    savefig("%s.combined_hist.png" % name)


def plot_runs(runs_df):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""


### A few extra definitions ###

def estimate_read_length(dataframe):
    """
    Takes in a dataframe of starttime and end time.
    Returns a series using the same index of the estimated length of a given read.
    """
    speed = 0.450 #bases per millisecond
    # Return a timedelta object converted into a float
    return speed * (dataframe["EndTime"] - dataframe["TimeTime"]).milliseconds()


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

def x_hist_to_human_readable(x, position):
    # Convert distribution to base pairs
    if x == 0:
        return 0
    s = humanfriendly.format_size(x, binary=False)
    return reformat_human_friendly(s)

def main():
    # Get args and samplesheet
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    ip_config = config_to_pd(args.ip_config)
    samples = get_samples(samplesheet)
    # Get runs for each sample
    for sample in samples:
        sample.runs = get_runs(sample.df, ip_config)
    # Get metadata for each run
    for sample in samples:
        for run in samples.runs:
            run.concatenated_data = get_metadata(run.path, run.slave.connect())
            run.concatenated_data["FlowcellID"] = run.flowcellID
        sample.concatenated_data = merge_metadatas([run.concatenated_data for run in samples.runs])
        sample.concatenated_data["SampleName"] = sample.name
    # Plot for each individual sample:
    for sample in samples:
        plot_single(sample.name, sample.concatenated_data.copy())
        plot_runs(sample.name, sample.concatenated_data.copy())
        for run in sample.runs:
            plot_single(run.flowcellID, run.concatenated_data)
    # Plot collective
    plot_samples(merge_metadatas([sample.concatenated_data for sample in samples]))