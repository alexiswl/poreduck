#!/usr/bin/env python2

"""
This script is used to move fast5 files generate by MinKNOW 1.3 and below
into subfolders of 4000 as per MinKNOW 1.4+ which makes it easier data to handle.
It uses a parallelisation technique to move the files more quickly into each of the subfolders.
"""

import argparse  # For importing arguments
import os  # Get file lists, check directories
import itertools  # Iterate through fast5file and mv commands
import subprocess  # Run mv commands
import sys  # For exiting with errors
import pandas as pd
from itertools import izip

# Global variables
READS_DIR = ""


def main():
    args = get_args()
    set_directories(args)
    directory_list = move_fast5_files(args)
    archive_folders(args, directory_list)


def get_args():
    """Use the argparse command to retrieve the reads directory input,
    and the archive toggle and number of threads set by the user"""
    parser = argparse.ArgumentParser(description="Moves fast5 files from folder into subdirectories")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="The directory with all the fast5 files in them")
    parser.add_argument("--archive", default=False, dest='archive', action='store_true',
                        help="Zip up folder once complete.")
    parser.add_argument("--num_threads", default=5,
                        help="Number of mv commands we will allow to be running at any given time.")

    args = parser.parse_args()
    return args


def set_directories(args):
    """ Make sure reads directory exists.
        Then change to reads directory."""
    global READS_DIR
    READS_DIR = args.reads_dir
    if not os.path.isdir(READS_DIR):
        sys.exit("%s not a directory" % READS_DIR)
    READS_DIR = os.path.abspath(READS_DIR) + "/"
    os.chdir(READS_DIR)


def move_fast5_files(args):
    """ Move fast5 files to subfolders.
    fast5_files are sorted by their modification time
    We use the subprocess command to parallelise the set of mv commands
    """
    # Create pandas dataframe with x columns.
    fast5_df = pd.DataFrame(columns=['fast5_file', 'subfolder', 'mv_command'])

    fast5_df['fast5_file'] = [fast5_file for fast5_file in os.listdir(READS_DIR) if fast5_file.endswith(".fast5")]
    fast5_df['subfolder'] = [standardise_int_length(int(i / 4000)) for i in xrange(len(fast5_df))]
    fast5_df['mv_command'] = ["mv %s %s/" % (fast5_file, subfolder)
                              for fast5_file, subfolder in izip(fast5_df.fast5_file, fast5_df.subfolder)]

    subdirectories = fast5_df.subfolder.unique().tolist()
    print(subdirectories)
    for subdirectory in subdirectories:
        # Create directory
        if os.path.isdir(subdirectory):
            # If directory already exists, make sure nothing is inside
            if len(os.listdir(subdirectory)) > 0:
                sys.exit("Directory '%s' exists with files inside" % subdirectory)
            else:
                os.mkdir(subdirectory)

    processes = (subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 for cmd in fast5_df.mv_command.tolist())

    # We use the islice command to split our list of mv commands into five smaller lists.
    running_processes = list(itertools.islice(processes, args.num_threads))
    while running_processes:
        for i, process in enumerate(running_processes):
            if process.poll() is not None:  # Means that the process is complete!
                stdout, stderr = process.communicate()  # Get the output of the completed process
                if not stderr == "":  # Print stderr if it exists.
                    print stderr
                running_processes[i] = next(processes, None)
                # Run the next number in the list.
                if running_processes[i] is None:  # No more commands waiting to be processed.
                    del running_processes[i]  # Not a valid process.
                    break

    return subdirectories


def archive_folders(args, directory_list):
    """
    Use the tar command with pigz to compress each of the subfolders
    """
    # Archive each of the subfolders
    # If we haven't selected archive then we return immediately.
    if args.archive is None:
        return

    # Otherwise a simple tar command should do
    tar_commands = []
    for directory in directory_list:
        tar_commands.append("tar -cf - %s --remove-files | pigz -9 -p 8 > %s.tar.gz" %
                            (directory, directory))

    # Multi-thread our tar command
    processes = (subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 for cmd in tar_commands)

    # We use the islice command to split our commands into five smaller lists.
    running_processes = list(itertools.islice(processes, args.num_threads))

    while running_processes:
        for i, process in enumerate(running_processes):
            if process.poll() is not None:  # Means that the process is complete!
                stdout, stderr = process.communicate()  # Get the output of the completed process
                if not stderr == "":
                    print stderr
                running_processes[i] = next(processes, None)
                # Run the next number in the list.
                if running_processes[i] is None:  # No more commands waiting to be processed.
                    del running_processes[i]  # Not a valid process.
                    break


def standardise_int_length(my_integer):
    # Input of 15 returns 0015
    return "%04d" % int(my_integer)


# Run the main function.
main()

