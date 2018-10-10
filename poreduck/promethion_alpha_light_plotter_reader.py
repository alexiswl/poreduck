import os
import pandas as pd
import numpy as np


def get_summary_files(summary_dir):
    summary_files = [os.path.join(summary_dir, summary_file)
                     for summary_file in os.listdir(summary_dir)
                     if summary_file.endswith("sequencing_summary.txt")
                     or summary_file.startswith("sequencing_summary_")
                    ] 
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

    # Get the cumulative channel yield
    #print(dataset)
    dataset['channel_yield'] = get_channel_yield(dataset)

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


def get_channel_yield(dataset):
    # Get the yield per channel 
    return dataset.groupby(['channel'])['sequence_length_template'].cumsum()


def get_yield(dataset):
    # Get the yield datset
    return dataset['sequence_length_template'].cumsum()


def get_pass(dataset):
    # Determine if sequence passed quality
    return dataset['mean_qscore_template'].apply(lambda x: True if x > 9 else False)


def get_duration_ratio(dataset):
    # Return the length in bases over time in seconds.
    return dataset.apply(lambda x: np.nan if x.template_duration == 0 else x.sequence_length_template/x.template_duration, axis='columns')


def get_events_ratio(dataset):
    # Return the events-per-base ratio
    return dataset.apply(lambda x: np.nan if x.sequence_length_template == 0 else x.num_events / x.sequence_length_template, axis='columns')
