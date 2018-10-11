import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
from scipy import stats
import matplotlib.pyplot as plt
import humanfriendly
from matplotlib.ticker import FuncFormatter
from matplotlib.pylab import savefig
from datetime import timedelta

import seaborn as sns


# Plot yield
def plot_yield(dataset, name, plots_dir):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""
    # Plot total yield for the sample
    # Yield plot
    # Set up plotting structure
    fig, ax = plt.subplots(1)

    # Plot setting start_time_float as axis index
    dataset.set_index("start_time_float_by_sample")["yield"].plot(ax=ax)

    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_yield_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set x and y labels
    ax.set_title("Yield over time for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative Yield")

    # Format nicely
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.yield.png" % name))
    plt.close('all')


# Plot reads
def plot_reads(dataset, name, plots_dir):
    """Plot an estimated yield plot and a histogram plot for each sample but by each flowcell"""
    # Plot total number of reads for the sample
    # Set up plotting structure
    fig, ax = plt.subplots(1)

    # Plot setting start_time_float as axis index
    dataset.set_index("start_time_float_by_sample")["read_count"].plot(ax=ax)

    # Set x and y ticks
    ax.yaxis.set_major_formatter(FuncFormatter(y_count_to_human_readable))
    ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set x and y labels
    ax.set_title("Read count over time for %s" % name)
    ax.set_xlabel("Time in (HH:MM)")
    ax.set_ylabel("Cumulative read count")

    # Format nicely
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.reads.png" % name))
    plt.close('all')
 


def plot_flowcell(dataset, name, plots_dir):
    # Set up plots 
    fig, ax = plt.subplots()
    fig.set_size_inches(15, 7)

    # Use the formatter we used for the yield plots.
    formatter_y = FuncFormatter(y_yield_to_human_readable)

    # Create the values that make up the numbers on the far-right column of the grid.
    c_w = 10
    c_l = 25
    c_num = 12

    # Create the array
    channels_by_order_array = np.array([[c_no * c_w * c_l + c_w * l_no + w_no + 1
                                         for c_no in np.arange(c_num)
                                         for w_no in np.arange(c_w)]
                                        for l_no in np.arange(c_l)])

    # Use the minknow_column_order function which reference the far-right column for a given row

    # to fill in the rest of the values for each row.
    channels_by_yield_array = np.zeros(channels_by_order_array.shape)

    # Sum the values for each channel.
    channels_by_yield_df = pd.DataFrame(dataset.groupby("channel")['channel_yield'].max())

    # Reset the index and have channel as a column instead of the index.
    channels_by_yield_df.reset_index(level='channel', inplace=True) 

    # Iterate through each row of the yield by channel dataframe.
    for yield_row in channels_by_yield_df.itertuples():
        channel_index = [(ix, iy)
                         for ix, row in enumerate(channels_by_order_array)
                         for iy, i in enumerate(row)
                         if int(i) == int(yield_row.channel)][0]

        # Assign channel yield to position in MinKNOW
        channels_by_yield_array[channel_index] = yield_row.channel_yield

    # Plot heatmap
    sns.heatmap(channels_by_yield_array,
                # Remove labels from side, they're not useful in this context.
                xticklabels=False,
                yticklabels=False,
                ax=ax,
                # Prevent extreme values from over-scaling the sidebar.
                robust=True,
                # Use the greens scale but in reverse, similar to MinKNOW.
                cmap=sns.diverging_palette(255, 133, l=60, as_cmap=True),
                # Format keyword args for the side bar.
                cbar_kws={"format": formatter_y,
                          "label": "Bases per channel"})
    
    # Create three lines down the middle as shown in PromethION MinKNOW.
    [ax.axvline([x], color='white', lw=5) for x in [30, 60, 90]]
    
    # Nice big title!
    ax.set_title("Map of Yield by Channel for %s" % name, fontsize=25)
    
    # Ensure labels are not missed.
    fig.tight_layout()

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.flowcellmap.png" % name))
    plt.close('all')


# Plot histogram
def plot_hist(dataset, name, plots_dir):
    # Set globals
    num_bins = 50
    max_quantile = 0.999

    # Open up plotting frame
    fig, ax = plt.subplots(1)

    # Get linspacing of histogram
    max_length = dataset['sequence_length_template'].quantile(max_quantile)
    trimmed = dataset.query("sequence_length_template < %s" % max_length)['sequence_length_template']

    # Develop the linspace
    bins = np.linspace(start=0, stop=trimmed.max(), num=num_bins)
    bin_width = bins[1] - bins[0]

    # Plot weighted histogram
    trimmed.plot(kind="hist", ax=ax, density=1, bins=bins, alpha=0.6, weights=trimmed)

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

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.hist.png" % name))
    plt.close('all')


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
    y = round(y, 3)
    s = humanfriendly.format_size(y, binary=False)
    return reformat_human_friendly(s)


def y_count_to_human_readable(y, position):
    # Use the same as y yield but strip off the last 'b'
    if y == 0:
        return 0
    y = round(y, 3)
    s = humanfriendly.format_size(y, binary=False)
    s = reformat_human_friendly(s).rstrip('b')
    return s


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


def plot_events_ratio(dataset, name, plots_dir):

    # Seaborn nomenclature for reg/lm plots are a little different 

    # Plot setting start_time_float as axis index
    max_quantile = 0.99

    # Trim the events ratio 
    max_ratio = dataset['events_ratio'].quantile(max_quantile)
    trimmed = dataset.query("events_ratio < %s" % max_ratio)
    
    # Set the background style for the plot
    sns.set_style('darkgrid')

    # Generate the plot 
    g = sns.lmplot(x='start_time_float_by_sample', y='events_ratio', data=trimmed,
                   hue='qualitative_pass', hue_order=['Passed', 'Failed'],
                   x_estimator=np.mean, truncate=True, x_bins=10, scatter_kws={'alpha': 0.1}, 
                   legend=False)

    # Create legend and rename
    leg_title = "Read Quality"
    leg = g.ax.legend(title=leg_title, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    # Zero base y-axis
    y_max = trimmed['events_ratio'].mean() * 2
    g.set(ylim=(0, y_max))

    # Set x and y labels
    g.set_axis_labels("Time in (HH:MM)", "Events ratio (events / base)")

    # Set title
    g.fig.suptitle("Events Ratio Graph for %s" % name)

    # Set x and y ticks:
    for ax in g.axes[0]:
        ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))
        ax.yaxis.set_major_formatter(FuncFormatter(y_count_to_human_readable))

    # Format nicely
    g.fig.tight_layout()

    # Reduce the plot size to make way for the suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save and close figure
    savefig(os.path.join(plots_dir, "%s.events_ratio.png" % name))
    plt.close("all")


def plot_quality_per_speed(dataset, name, plots_dir):

    # Seaborn nomenclature for joint plots are a little different
    sns.set_style("dark")
    g = sns.jointplot(x='pore_speed', y='mean_qscore_template',
                      data=dataset, kind='hex')

    # Add pearson stat
    g.annotate(stats.pearsonr)

    # Set axis labels
    g.set_axis_labels("Pore Speed (b/s)", "Mean q-score Template")

    # Set title
    g.fig.suptitle("Pore Speed against q-score template for %s" % name)

    # Format nicely.
    g.fig.tight_layout()

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)
    
    # Save and close the figure
    savefig(os.path.join(plots_dir, "%s.speed_vs_qscore.png" % name))
    plt.close('all')


def plot_pore_speed(dataset, name, plots_dir):
    # Plot setting start_time_float as axis index

    # Seaborn nomenclature for lmplots/regplots are a little different
    sns.set_style('darkgrid')

    g = sns.lmplot(x='start_time_float_by_sample', y='pore_speed', data=dataset,
                   hue='qualitative_pass', hue_order=["Passed", "Failed"],
                   x_estimator=np.mean, truncate=True, x_bins=10, scatter_kws={'alpha': 0.1},
                   legend=False)

    # Create legend with new alpha
    leg_title="Read Quality"
    leg = g.ax.legend(title=leg_title, framealpha=0.5)    
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    # Zero base y-axis
    y_max = dataset['pore_speed'].mean() * 2
    g.set(ylim=(0, y_max))

    # Set axis labels
    g.set_axis_labels("Time (HH:MM)", "Pore Speed (bases / second)")

    # Set axis formats
    for ax in g.axes[0]:
        ax.xaxis.set_major_formatter(FuncFormatter(x_yield_to_human_readable))

    # Set title
    g.fig.suptitle("Pore speed over time")

    # Reduce plot to make room for suptitle
    g.fig.subplots_adjust(top=0.95)

    # Save figure
    savefig(os.path.join(plots_dir, "%s.pore_speed.png" % name))
    plt.close('all')


def convert_sample_time_columns(dataset):
    # Use the utc in the fastq file to work around restarts
    min_start_time = dataset['start_time_utc'].min()
    dataset['start_time_timedelta_by_sample'] = dataset['start_time_utc'].apply(lambda x: x - min_start_time)

    # Convert to float because matplotlib doesn't seem to do timedelta on the x axis well.
    # Need to divide by another timedelta object in order to get float
    dataset['start_time_float_by_sample'] = dataset['start_time_timedelta_by_sample'].apply(lambda x: x / timedelta(seconds=1))

    # Sort values to start_time_float_by_sample (to assist yield plotting)
    dataset.sort_values(['start_time_float_by_sample'], inplace=True)

    # Reset index values to match
    dataset.reset_index(drop=True, inplace=True)

    # Return
    return dataset


def get_channel_yield(dataset):
    # Get the yield per channel
    return dataset.groupby(['channel'])['sequence_length_template'].cumsum()


def get_yield(dataset):
    # Get the yield dataset
    return dataset['sequence_length_template'].cumsum()


def get_read_count(dataset):
    # Get the read count dataset
    return dataset.reset_index()['index']


def plot_data(dataset, name, plots_dir):
    # Add in the start_time_float_by_sample (allows us to later iterate through plots by sample.
    dataset = convert_sample_time_columns(dataset)

    # Get read_count column
    dataset['read_count'] = get_read_count(dataset)

    # Get yield column
    dataset['yield'] = get_yield(dataset)

    # Get the cumulative channel yield
    dataset['channel_yield'] = get_channel_yield(dataset)

    #dataset.to_csv("dataset_test.csv", index=False, header=True)

    # Plot things
    # Matplotlib base plots
    plot_yield(dataset, name, plots_dir)
    plot_reads(dataset, name, plots_dir)
    plot_hist(dataset, name, plots_dir)

    # Seaborn plots
    plot_flowcell(dataset, name, plots_dir)
    plot_pore_speed(dataset, name, plots_dir)
    plot_quality_per_speed(dataset, name, plots_dir) 
    plot_events_ratio(dataset, name, plots_dir)

