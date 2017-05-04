#!/usr/bin/env python

import argparse
import os
import itertools
import subprocess
import sys

# Global variables
READS_DIR = ""
NUM_THREADS = 5


def main():
    args = get_args()
    set_directories(args)
    move_fast5_files(args)


def get_args():
    parser = argparse.ArgumentParser(description="Moves fast5 files from folder into subdirectories")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="The directory with all the fast5 files in them")
    parser.add_argument("--archive", default=False, dest='archive',
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
    os.path.abspath(READS_DIR) + "/"
    os.chdir(READS_DIR)


def move_fast5_files(args):
    """ Move fast5 files to subfolders.
    fast5_files are sorted by their modification time
    We use the subprocess command to parallelise the set of mv commands
    """

    fast5_files = [fast5_file for fast5_file in os.listdir(READS_DIR)
                   if fast5_file.endswith(".fast5")]
    fast5_files.sort(key=lambda x: os.path.getmtime(x))

    subdirectory_list = []
    subdirectory = 0
    mv_fast5_commands = []

    for f_index, fast5_file in itertools.izip(range(0, len(fast5_files)), fast5_files):
        subdirectory_list[subdirectory] = subdirectory_list[subdirectory] + "," + fast5_file
        if f_index + 1 % 4000 == 0:
            mv_fast5_commands.append(' '.join(subdirectory_list[subdirectory].split(","))
                                     + subdirectory)
            subdirectory += 1
        elif f_index + 1 == len(fast5_files):
            mv_fast5_commands.append(' '.join(subdirectory_list[subdirectory].split(","))
                                     + subdirectory)

    # Generate processing list
    # Talk to the shell. Note these commands won't run just yet.
    processes = (subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 for cmd in mv_fast5_commands)

    # We use the islice command to split our list of mv commands into five smaller lists.
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

main()
