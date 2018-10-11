import os
import pandas as pd
import numpy as np
import gzip
from Bio import SeqIO

def get_summary_files(summary_dirs):
    summary_files = [os.path.join(summary_dir, summary_file)
                     for summary_dir in summary_dirs
                     for summary_file in os.listdir(summary_dir)
                     if summary_file.endswith(".txt")
                     and "sequencing_summary" in summary_file]

    return summary_files


def get_fastq_files(fastq_dirs):
    fastq_files = [os.path.join(fastq_dir, fastq)
                   for fastq_dir in fastq_dirs
                   for fastq in os.listdir(fastq_dir)
                   if fastq.endswith(".fastq.gz")
                  ]

    return fastq_files


def get_series_from_seq(record):
    """Takes in a seq record and returns dataframe with index"""
    # Series index
    index = ["read_id", "run_id", "sample_id", "read", "channel", "start_time_utc"]
    # Get metadata
    fastq_id = record.id.split()[0].lstrip("@")
    row_as_dict = dict(x.split("=") for x in record.description.split()[1:])

    return pd.Series([fastq_id, row_as_dict['runid'],
                      row_as_dict['sampleid'], row_as_dict['read'],
                      row_as_dict['ch'], row_as_dict['start_time']],
                     index=index)


def get_fastq_dataframe(fastq_file, is_gzipped=True):
    """Use get_series_from_seq in list comprehension to generate dataframe then transpose"""
    # Open fastq file and return list comprehension
    # PD concat merges series as columns, we then transpose.
    try:
        if not is_gzipped:
            with open(fastq_file, "rt") as handle:
                fastq_df = pd.concat([get_series_from_seq(record)
                                      for record in SeqIO.parse(handle, "fastq")],
                                     sort=True,
                                     axis='columns').transpose()
        else:
            with gzip.open(fastq_file, "rt") as handle:
                fastq_df = pd.concat([get_series_from_seq(record)
                                      for record in SeqIO.parse(handle, "fastq")],
                                     sort=True,
                                     axis='columns').transpose()
        # Specify types for each line
        numeric_cols = ["read", "channel"]
        fastq_df[numeric_cols] = fastq_df[numeric_cols].apply(pd.to_numeric, axis='columns')
        # Convert StartTime to date
        fastq_df['start_time_utc'] = pd.to_datetime(fastq_df['start_time_utc'])
        return fastq_df
    except ValueError:
        print("Value error when generating dataframe for %s. Unknown cause of issue." % fastq_file)
        return pd.DataFrame(columns=["read_id", "run_id", "sample_id", "read", "channel", "start_time_utc"])


def read_fastq_datasets(fastq_files):
    # Read in each of the fastq files and retrieve header info
    dataset = pd.concat([get_fastq_dataframe(fastq_file, is_gzipped=True)
                         for fastq_file in fastq_files],
                        sort=True, ignore_index=True)
    # Return dataset
    return dataset


def read_summary_datasets(sequencing_summary_files):
    # Read in each of the csv files.
    dataset = [pd.read_csv(sequencing_summary_file, sep="\t", header=0)
               for sequencing_summary_file in sequencing_summary_files]

    # Merge the list of datasets
    dataset = merge_summary_dataset(dataset)

    # Reset the dtypes for the time columns
    dataset = set_summary_time_dtypes(dataset)

    # Get the cumulative channel yield
    dataset['channel_yield'] = get_channel_yield(dataset)

    # Get pass column
    dataset['pass'] = get_pass(dataset)

    # Get qualitative pass
    dataset['qualitative_pass'] = get_qualitative_pass(dataset)

    # Get duration ratio
    dataset['pore_speed'] = get_duration_ratio(dataset)

    # Get events ratio
    dataset['events_ratio'] = get_events_ratio(dataset)

    # Return the object
    return dataset


def merge_summary_dataset(dataset):
    """
    :rtype: pd.DataFrame
    """
    # Merge a list of datasets into a single pandas dataframe
    return pd.concat(dataset, sort=True, ignore_index=True)


def set_summary_time_dtypes(dataset):
    """
    :rtype: pd.DataFrame
    """
    # Return the time datasets as appropriate
    dataset.rename(columns={"start_time": "start_time_float",
                            "template_start": "template_time_float"}, inplace=True)

    # Add start_time_float template_time_float
    dataset['start_time_timedelta'] = pd.to_timedelta(dataset['start_time_float'], unit='s')
    dataset['template_start_time_timedelta'] = pd.to_timedelta(dataset['template_time_float'], unit='s')

    return dataset


def sort_summary_dataset(dataset):
    """
    :rtype: pd.DataFrame
    """
    # Sort the dataset by start_time (in seconds)
    return dataset.sort_values(['start_time_timedelta'])


def get_channel_yield(dataset):
    # Get the yield per channel 
    return dataset.groupby(['channel'])['sequence_length_template'].cumsum()


def get_yield(dataset):
    # Get the yield datset
    return dataset['sequence_length_template'].cumsum()


def get_pass(dataset):
    # Determine if sequence passed quality
    return dataset['mean_qscore_template'].apply(lambda x: True if x > 9 else False)


def get_qualitative_pass(dataset):
    # Describe the pass (Passed / Failed)
    return dataset['pass'].apply(lambda x: 'Passed' if x == True else "Failed")


def get_duration_ratio(dataset):
    # Return the length in bases over time in seconds.
    return dataset.apply(lambda x: np.nan if x.template_duration == 0 else x.sequence_length_template/x.template_duration, axis='columns')


def get_events_ratio(dataset):
    # Return the events-per-base ratio
    return dataset.apply(lambda x: np.nan if x.sequence_length_template == 0 else x.num_events / x.sequence_length_template, axis='columns')
