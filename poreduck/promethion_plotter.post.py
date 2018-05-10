#!/usr/bin/env python3

from promethion_starter import samplesheet_to_pd, config_to_pd
import pandas as pd
import os
import argparse
import paramiko
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import humanfriendly
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
import numpy as np
import seaborn as sns


series_columns = ["Name", "Channel", "Read", "RNumber", # From file name
                  "MuxID", "StartTime", "EndTime"]  # From fast5 file - requried for plots


class Sample:
    def __init__(self, sample_name, sample_df):
        self.name = sample_name
        self.df = sample_df
        self.flowcells = []


class Flowcell:
    def __init__(self, sample_name, flowcell_df, flowcell_ID):
        self.sample_name = sample_name
        self.flowcell_df = flowcell_df
        self.flowcell_ID = flowcell_ID
        self.runs = []


class Run:
    def __init__(self, sample_name,
                 grnwch_mux_start_date, grnwch_mux_start_time,
                 grnwch_seq_start_date, grnwch_seq_start_time,
                 slurm_id, flowcell_id, pca_path):
        self.sample_name = sample_name
        self.grnwch_mux_start_date = grnwch_mux_start_date
        self.grnwch_mux_start_time = grnwch_mux_start_time
        self.grnwch_seq_start_date = grnwch_seq_start_date
        self.grnwch_seq_start_time = grnwch_seq_start_time
        self.slurm_id = slurm_id
        self.flowcell_id = flowcell_id
        self.concatenated_data = None
        self.path = os.path.join(pca_path, slurm_id, "reads", '_'.join([self.grnwch_seq_start_date,
                                                                        self.grnwch_seq_start_time,
                                                                        self.sample_name]))
        self.metadata_path = os.path.join(self.path, "metadata", "merged")
        print(self.path)


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
    parser.add_argument("--pca_dir", type=str, required=True,
                        help="/PATH/TO/PCAXXXX")
    args = parser.parse_args()
    return args


def get_samples(samplesheet_df, pca_dir):
    # Get samples
    samples = []
    for sample_name in samplesheet_df["SampleName"].unique().tolist(): 
        samples.append(Sample(sample_name, samplesheet_df.query("SampleName=='%s'" % sample_name)))
    return samples 


def get_flowcells(sample_df, pca_dir):
    # Get flowcells from a given sample.
    flowcells = []
    for flowcell_id in sample_df["FlowcellID"].unique().tolist():
        sample_name = sample_df["SampleName"].unique().item()
        flowcells.append(Flowcell(sample_name, sample_df.query("FlowcellID=='%s'" % flowcell_id), flowcell_id))
    return flowcells


def get_runs(flowcell_df, pca_dir):
    """Return a list of class run, each with a list of associated dataframes"""
    runs = []
    for index, row in flowcell_df.iterrows():
        runs.append(Run(row.SampleName,
                        row.GrnwchMuxStartDate, row.GrnwchMuxStartTime,
                        row.GrnwchSeqStartDate, row.GrnwchSeqStartTime,
                        row.SlurmID, row.FlowcellID, pca_dir))
    return runs


def get_metadata(run_path, flowcell_id):
    """Return a list of pandas dataframes imported over paramiko"""
    metadata_path = os.path.join(run_path, "metadata", "merged")
    # Get list of tsv files
    tsv_list = [tsv
                for tsv in os.listdir(path=metadata_path)
                if tsv.endswith(".tsv")]
    # Read in dataframe via read_csv()
    metadatas = []
    for tsv in tsv_list:
        metadata = pd.read_csv(os.path.join(metadata_path, tsv), header=0, sep="\t",
                               parse_dates=["StartTime", "EndTime"])
        # Convert the following columns to numeric
        metadata[["Channel", "Read", "MuxID", "RNumber"]] = metadata[["Channel", "Read", "MuxID", "RNumber"]].apply(pd.to_numeric)
        # Add to list of dataframes 
        metadatas.append(metadata)
        
        # Add the flowcell ID
        metadata['FlowcellID'] = flowcell_id
    # Return a merged version of the metadata 
    return merge_metadatas(metadatas)


def merge_metadatas(metadatas):
    """Merge a list of dataframes for a given run"""
    # Calculate the Run Duration Time and Floats by flowcell id prior to merging:
    metadatas_all = []
    for flowcell in list(set([metadata['FlowcellID'].unique().item() for metadata in metadatas])):
        metadatas_fc = pd.concat([metadata 
                                   for metadata in metadatas 
                                   if metadata['FlowcellID'].unique().item() == flowcell],
                                  ignore_index=True)
        run_start = min(metadatas_fc["StartTime"])
        metadatas_fc["RunDurationTime"] = metadatas_fc["EndTime"].apply(lambda x: x - run_start)
        metadatas_fc["RunDurationFloat"] = metadatas_fc["RunDurationTime"].apply(lambda x: x.total_seconds())
        metadatas_all.append(metadatas_fc)
    return pd.concat(metadatas_all, ignore_index=True)


def plot_samples(names, samples_df):
    """Plot an estimated yield plot and a histogram plot comparing each flowcell"""
    # Need to configure the index properly for this one and incorporate just the required dataframes
    # Yield plot
    # Set up plotting structure
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Set and sort the index by time
    samples_df.set_index("RunDurationFloat", inplace=True)
    samples_df.sort_index(inplace=True)
    # Generate a cumulative yield
    samples_df["CumYield"] = samples_df.groupby(["SampleName"])["SeqLength"].apply(lambda x: x.cumsum())
    # Pivot dataframe using SampleName as pivot, CumYield as column. Save as yield_df
    yield_df = samples_df.pivot(index="SampleName",
                                columns='CumYield',
                                values='CumYields').filter(regex=r'^CumYields\.', axis=1)
    # Plot cumulative yield over time.
    yield_df.plot(ax=ax)
    # Define axis formatters
    ax.xaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_title("Yield over time by sample")
    # Set x and y limits
    #ax.set_xlim(samples_df["DurationTime"].min(), samples_df["DurationTime"].max())
    ax.set_ylim(ymin=0)
    # Ensure labels are not missed
    fig.tight_layout()
    savefig("%s.combined_yield.png" % '.'.join([samples_df['SampleName'].tolist()]))


def plot_histograms(names, samples_df):
    # Generate a weighted histogram for each sample
    # Histogram
    # Set up plotting structure
    plt.close('all')
    fig, ax = plt.subplots(1)
    num_bins = 50
    # Reset the index of the dataframe
    samples_df.reset_index(drop=True, inplace=True)
    # Clip data and plot histogram
    lengths = samples_df['SeqLength']
    samples_df['SeqLength'] = samples_df[lengths < lengths.quantile(0.9999)]
    # Pivot dataframe using SampleName as pivot, save the SeqLengths
    hist_df = samples_df.pivot(index="SampleName",
                               columns='SeqLength',
                               values='SeqLengths').filter(regex=r'^SeqLengths\.', axis=1)
    hist_df.plot(type='hist', ax=ax, normed=1, facecolor='blue', alpha=0.75)
    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        if y == 0:
            return 0
        s = humanfriendly.format_size(bin_width * lengths.sum() * y, binary=False)
        return reformat_human_friendly(s)
    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable_seq))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
    # Set titles
    ax.set_title("Read Distribution Graph for %s" % names)
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_ylabel("Bases per bin")
    # Ensure labels are not mised.
    fig.tight_layout()
    savefig("%s.combined_hist.png" % names)


def plot_flowcell(name, sample_df):
    """Plot an estimated yield plot and a histogram plot for the sample""" 
    # Histogram
    plt.close('all')
    fig, ax = plt.subplots(1)
    num_bins = 50
    # Clip data and plot histogram
    lengths = sample_df['SeqLength']
    lengths = lengths[lengths < lengths.quantile(0.999)]
    bins = np.linspace(start=0, stop=lengths.max(), num=num_bins)
    bin_width = bins[1] - bins[0]
    lengths.plot(kind='hist', ax=ax, normed=1, bins=bins, color='blue', alpha=0.75, weights=lengths)
    # Set the axis formatters
    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        if y == 0:
            return 0
        s = humanfriendly.format_size(bin_width * lengths.sum() * y, binary=False)
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
    savefig("%s.combined_hist.png" % name)

    # Plot flowcell
    
    # Split into chunks of 64 (rows of 4)
    c_w = 10
    c_l = 25
    c_num = 12
    
    # Create the values that make up the numbers of the far-right column of the grid
    channels_by_order_array = np.array([[c_no*c_w*c_l + c_w*l_no + w_no + 1
                                         for c_no in np.arange(c_num)
                                         for w_no in np.arange(c_w)]
                                         for l_no in np.arange(c_l)])
    
    # Initialise the array with zeros
    channels_by_yield_array = np.zeros(channels_by_order_array.shape)

    # Sum the values for each channel
    channels_by_yield_df = pd.DataFrame(sample_df.reset_index().groupby("Channel")['SeqLength'].sum())
    
    # Reset the index and have channel as a column instead of the index.
    channels_by_yield_df.reset_index(level=0, inplace=True)
    
    # Iterate through each row of the yield by channel dataframe.
    for yield_row in channels_by_yield_df.itertuples():
        channel_index = [(ix, iy)
                         for ix, row in enumerate(channels_by_order_array)
                         for iy, i in enumerate(row)
                         if int(i) == int(yield_row.Channel)][0]
        # Assign channel yield to position in MinKNOW
        channels_by_yield_array[channel_index] = yield_row.SeqLength

    # Close any previous plots
    plt.close('all')
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
    # Create three lines down the middle as shown in PromethION MinKNOW.
    [ax.axvline([x], color='white', lw=5) for x in [30, 60, 90]]
    # Nice big title!
    ax.set_title("Map of Yield by Channel", fontsize=25)
    # Ensure labels are not missed.
    fig.tight_layout() 
    savefig("%s.flowcellmap.png" % name)


def plot_sample(sample_df):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""
    # Plot total yield for the sample
    # Yield plot
    # Set up plotting structure
    name = sample_df["SampleName"].unique().item()
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Set and sort the index by time
    sample_df.dropna(inplace=True, how='any')
    sample_df.sort_values(by="RunDurationTime", inplace=True)
    sample_df.set_index("RunDurationFloat", inplace=True)
    # Plot cumulative yield over time.
    sample_df["SeqLength"].cumsum().plot(ax=ax)
    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_title("Yield for sample %s over time (combined)" % name)
    # Set x and y limits
    ax.set_ylim(ymin=0)
    # Ensure labels are not missed
    fig.tight_layout()
    savefig("%s.combined_yield.png" % name)
    
    # Plot yield by quality
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Take just the Class and Seqlength columns (And Date/name as index)
    q_classes = {"pass": "Green", "fail": "Red"}
    for q_class, colour in q_classes.items():
        sample_df.query("Class=='%s'" % q_class)["SeqLength"].cumsum().plot(ax=ax, color=colour)
    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_title("Yield for sample %s over time (by quality)" % name)
    # Set x and y limits
    ax.set_ylim(ymin=0)
    # Add legend (rename legend)
    ax.legend(q_classes.keys())
    # Ensure labels are not missed
    fig.tight_layout()
    savefig("%s.combined_yield_by_flowcell.class.png" % name)

    # Plot yield per flowcell
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Get the flowcell list and iterate through the line styles
    flowcells = sample_df['FlowcellID'].unique().tolist()
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    
    if len(linestyles) < len(flowcells):
        sys.exit("Error, don't have enough linestyles yet to plot all flowcells")
    for flowcell, linestyle in zip(flowcells, linestyles):
        sample_df.query("FlowcellID=='%s'" % flowcell)["SeqLength"].cumsum().plot(ax=ax, linestyle=linestyle)

    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_title("Yield for sample %s over time (by flowcell)" % name)
    # Set x and y limits
    ax.set_ylim(ymin=0)
    # Add legend (rename legend)
    ax.legend(flowcells)
    # Ensure labels are not missed
    fig.tight_layout()
    savefig("%s.combined_yield.all_flowcells.png" % name)

    # Plot histogram
    # Histogram
    plt.close('all')
    fig, ax = plt.subplots(1)
    num_bins = 50
    # Get linspacing of histogram
    trimmed = sample_df['SeqLength'][sample_df.SeqLength < sample_df.SeqLength.quantile(0.999)]
    bins = np.linspace(start=0, stop=trimmed.max(), num=num_bins)
    bin_width = bins[1] - bins[0]
    # Clip data and plot histogram
    for color, flowcell in enumerate(flowcells):
        lengths = sample_df.query("FlowcellID=='%s'" % flowcell)['SeqLength']
        lengths = lengths[lengths < lengths.quantile(0.999)]
        lengths.plot(kind='hist', ax=ax, normed=1, color='C%d' % color, bins=bins, alpha=0.6, weights=lengths)
        # Set the axis formatters
    def y_hist_to_human_readable_seq(y, position):
        # Convert distribution to base pairs
        if y == 0:
            return 0
        s = humanfriendly.format_size(bin_width * sample_df["SeqLength"].sum() * y, binary=False)
        return reformat_human_friendly(s)
    ax.yaxis.set_major_formatter(FuncFormatter(y_hist_to_human_readable_seq))
    ax.xaxis.set_major_formatter(FuncFormatter(x_hist_to_human_readable))
    # Set titles
    ax.set_title("Read Distribution Graph for %s" % name)
    ax.grid(color='black', linestyle=':', linewidth=0.5)
    ax.set_ylabel("Bases per bin")
    ax.set_xlabel("Read length")
    # Create legend
    ax.legend(flowcells)
    # Ensure labels are not missed.
    fig.tight_layout()
    savefig("%s.combined_hist.png" % name)  
    # Plot yield per flowcell per quality
    # Same method as before, plot yield by linetype and then by quality.
    # Plot yield per flowcell
    plt.close('all')
    fig, ax = plt.subplots(1)
    # Get the flowcell list and iterate through the line styles
    flowcells = sample_df['FlowcellID'].unique().tolist()
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    q_classes = {"pass": "Green", "fail": "Red"}    
    # Flowcell by line type, quality by color
    if len(linestyles) < len(flowcells):
        sys.exit("Error, don't have enough linestyles yet to plot all flowcells")
    for flowcell, linestyle in zip(flowcells, linestyles):
        for key, value in q_classes.items():
            sample_df.query("FlowcellID=='%s' & Class=='%s'" % (flowcell, key))["SeqLength"].cumsum().plot(ax=ax, linestyle=linestyle, color=value)

    # Define axis formatters
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
    # Set x and y labels
    ax.set_xlabel("Duration of run (HH:MM)")
    ax.set_ylabel("Yield")
    ax.set_title("Yield for sample %s over time (by flowcell/quality)" % name)
    # Set x and y limits
    ax.set_ylim(ymin=0)
    # Add legend (rename legend)
    ax.legend([flowcell+" "+ qclass 
               for flowcell in flowcells 
               for qclass in q_classes.keys()])
    # Ensure labels are not missed
    fig.tight_layout()
    savefig("%s.combined_yield.all_flowcells.quality.png" % name)


def print_stats(sample_df):
    percentiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    # Get total yield
    total_bp = sample_df['SeqLength'].sum()
    total_bp_h = reformat_human_friendly(humanfriendly.format_size(total_bp, binary=False))
    # Describe length
    length_describe = sample_df['SeqLength'].describe(percentiles=percentiles).to_string()
    # Describe quality
    qual_describe = sample_df['AvQual'].describe(percentiles=percentiles).to_string()

    # Reformat each of the methods such that they're rounded to two decimal places
    length_describe = '\n'.join([qual_line.split()[0].ljust(8) + "\t" + 
                                 "{:21.2f}".format(float(qual_line.split()[1]))
                                 for qual_line in length_describe.split("\n")])
    qual_describe = '\n'.join([qual_line.split()[0].ljust(8) + "\t" + 
                               "{:21.2f}".format(float(qual_line.split()[1]))
                               for qual_line in qual_describe.split("\n")])
    
    # Calculate the N50:
    nx = []
    seq_length_sorted_as_series = sample_df['SeqLength'].sort_values().reset_index(drop=True)
    seq_length_cumsum_as_series = seq_length_sorted_as_series.cumsum()
    for index, seq_value in seq_length_sorted_as_series.iteritems():
        if (seq_length_cumsum_as_series[index] <= total_bp*percentiles[len(nx)] <= 
          seq_length_cumsum_as_series[index+1]):
            nx.append(seq_value)
        if len(nx) == len(percentiles):
            # Found all the percentiles, no need to continue.
            break
    nx_h = [reformat_human_friendly(humanfriendly.format_size(nX_value, binary=False))
            for nX_value in nx]
 
    # Print these stats
    sample_name = sample_df["SampleName"].unique().item()
    with open(sample_name + ".stats.txt", 'a') as output_handle:
        # Print total basepairs
        output_handle.write("# Stats for sample '%s' #\n" % sample_name)
        output_handle.write("Total basepairs:\n")
        output_handle.write(f"\t{total_bp:16,.0f}\t|\t{total_bp_h.rjust(9)}\n")
        output_handle.write("Description of Read Lengths:\n")
        # Tab indent each of the descriptor lines
        output_handle.writelines(f"\t{len_line}\n"
                                 for len_line in length_describe.split("\n"))
        output_handle.write("Description of Read Qualities:\n")
        # Tab indent each of the descriptor lines
        output_handle.writelines(f"\t{qual_line}\n"
                                 for qual_line in qual_describe.split("\n"))
        output_handle.write("NX values:\n")
        output_handle.writelines(f"\tN{100*percentile:02.0f}:\t{nx_value:8,.0f}\t|\t{nx_h_value.rjust(9)}\n"
                                 for percentile, nx_value, nx_h_value in zip(percentiles, nx, nx_h)) 

    # Re-run this script for each flowcell
    flowcells = sample_df["FlowcellID"].unique().tolist()
    if len(flowcells) > 1:
        for flowcell in flowcells:
            with open(sample_name + ".stats.txt", 'a') as output_handle:
                output_handle.write("\n### Stats for flowcell '%s' ###\n" % flowcell)
            print_stats(sample_df.query("FlowcellID=='%s'" % flowcell))
    else:
        # Print the run-duration for each flowcell
        # Get run duration, from first read to last read.
        duration = sample_df["RunDurationFloat"].max()  # In seconds
        hours, remainder = divmod(duration, 3600)
        minutes, seconds = divmod(remainder, 60)
        run_duration_h = f"{hours} hours, {minutes} minutes, {seconds:2,.0f} seconds"
        with open(sample_name + ".stats.txt", 'a') as output_handle:
            output_handle.write(f"\t{duration:8,.1f} seconds\t|\t{run_duration_h}\n")
            
    
 
### A few extra definitions ###
def estimate_read_length(dataframe):
    """
    Takes in a dataframe of starttime and end time.
    Returns a series using the same index of the estimated length of a given read.
    """
    speed = 0.450  # bases per millisecond
    # Return a timedelta object converted into a float
    return speed * (dataframe["EndTime"] - dataframe["StartTime"]).microseconds


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


def x_yield_to_human_readable(x, position):
    # Convert time in seconds to hours or minutes 
    hours = int(x // 3600)
    minutes = int((x % 3600) // 60)
    seconds = int(x % 60)
    if x == 0:
        return 0
    s = f"{hours:02d}:{minutes:02d}"
    return s


def main():
    # Get args and samplesheet
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    ip_config = config_to_pd(args.ip_config)
    samples = get_samples(samplesheet, args.pca_dir)
    # Get runs for each sample
    for sample in samples: 
        sample.flowcells = get_flowcells(sample.df, args.pca_dir)
    for sample in samples:
        for flowcell in sample.flowcells:
            flowcell.runs = get_runs(flowcell.flowcell_df, args.pca_dir)
    # Get metadata for each run
    for sample in samples:
        for flowcell in sample.flowcells:
            for run in flowcell.runs:
                print(sample.name, run.flowcell_id)
                run.concatenated_data = get_metadata(run.path, run.flowcell_id)
                run.concatenated_data["FlowcellID"] = run.flowcell_id
                run.concatenated_data["SampleName"] = sample.name
            flowcell.concatenated_data = merge_metadatas([run.concatenated_data for run in flowcell.runs])
    # Plot for each individual sample:
    # Print stats for each individual sample
    for sample in samples:
        for flowcell in sample.flowcells:
            plot_flowcell(sample.name + "_" + flowcell.flowcell_ID, flowcell.concatenated_data)
        plot_sample(merge_metadatas([flowcell.concatenated_data 
                                     for flowcell in sample.flowcells]))
        if os.path.isfile(sample.name + ".stats.txt"):
            os.remove(sample.name + ".stats.txt")
        print_stats(merge_metadatas([flowcell.concatenated_data for flowcell in sample.flowcells]))
    # Plot collective
    plot_samples([sample.name 
                  for sample in samples],
                  merge_metadatas([flowcell.concatenated_data 
                                   for sample in samples 
                                   for flowcell in sample.flowcells]))

if __name__ == "__main__":
    main()
