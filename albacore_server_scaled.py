#!/usr/bin/env python

"""
This albacore script searches for tar zipped folders in a directory.
From here, it extracts the fast5 files and runs them through the albacore basecaller.

We can scale up the threads with qsub such that each of these tarred files can be run in parallel.

Each folder takes about 30 minutes to run on 10 threads. So with infinite threads and ram we could basecall everything
in half an hour!

The transfer_fast5_to_server.py script places a "TRANSFERRING" lock file on the folder.
This allows the transfer_fast5_to_server.py script to know when to stop looking.

"""

# Import necessary modules
import os  # General directory list stuff
import subprocess  # Executing each of the qsub commands.
import argparse  # For allowing users to configure the arguments used.
import sys  # For errors
import time  # For schnoozing!


# Set global and semi-global variables
TRANSFERRING = True
READS_DIR = ""
OUTPUT_DIR = ""
WORKING_DIR = ""
NUM_THREADS = 0
CHOSEN_CONFIG = ""
PARENT_DIRECTORY = ""
QSUB_LOG_DIR = ""

CONFIGS = {"FC106_RAD001": "FLO-MIN106_RAD001_linear.cfg",
           "FC106_LSK208_tc": "FLO-MIN106_LSK208_tc.cfg",
           "FC106_LSK108": "FLO-MIN106_LSK108_linear.cfg",
           "FC106_RAD002": "FLO-MIN106_RAD002_linear.cfg",
           "FC106_LSK208_2d": "FLO-MIN106_LSK208_2d.cfg"}


def main():
    # Basic house cleaning, get arguments, make sure they're legit.
    args = get_arguments()
    set_global_variables(args)
    check_directories()

    # While the transferring lock script exists.
    while TRANSFERRING:
        tarred_read_sets = get_tarred_files()

        for tarred_read_set in tarred_read_sets:
            extract_tarred_read_set(tarred_read_set)
            run_albacore(tarred_read_set)

        is_still_transferring()
        if len(tarred_read_sets) == 0:
            take_a_break()


def get_arguments():
    parser = argparse.ArgumentParser(
        description="The albacore-server scaled command incorporates qsub to spread the server load and rapidly " +
                    "generate data off the MinION.")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="/path/to/reads, should have a bunch of tar zipped files in it.")
    parser.add_argument("--config", type=str, choices=CONFIGS.keys(),
                        help="Pick a config")
    parser.add_argument("--output_dir", type=str, required=False, default=None,
                        help="Will be called 'basecalled' and sit adjacent to the reads folder if left blank.")
    parser.add_argument("--num_threads", type=int, required=False, default=5,
                        help="How many threads did you wish to use per parallel output, default is 5.")
    return parser.parse_args()


def set_global_variables(args):
    global READS_DIR, OUTPUT_DIR, WORKING_DIR, NUM_THREADS, CHOSEN_CONFIG
    READS_DIR = args.reads_dir
    if args.output_dir is not None:
        OUTPUT_DIR = args.output_dir
    CHOSEN_CONFIG = CONFIGS[args.config]
    NUM_THREADS = args.num_threads


def check_directories():
    global READS_DIR, OUTPUT_DIR, PARENT_DIRECTORY, QSUB_LOG_DIR
    if not os.path.isdir(READS_DIR):
        sys.exit("Error, %s does not exist" % READS_DIR)
    READS_DIR = os.path.abspath(READS_DIR) + "/"
    os.chdir(READS_DIR)

    PARENT_DIRECTORY = os.path.abspath(os.path.join(READS_DIR, os.pardir)) + "/"

    if OUTPUT_DIR is not None:
        OUTPUT_DIR = PARENT_DIRECTORY + "albacore/"
    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    QSUB_LOG_DIR = PARENT_DIRECTORY + "qsub_log/"

    if not os.path.isdir(QSUB_LOG_DIR):
        os.mkdir(QSUB_LOG_DIR)


def get_tarred_files():
    tarred_files = [READS_DIR + tarred_file for tarred_file in os.listdir(READS_DIR)
                    if tarred_file.endswith(".tar.gz")  # Is a zip file and
                    and tarred_file.replace(".tar.gz", "/")  # basecalled output folder
                    not in os.listdir(OUTPUT_DIR)]           # does not already exist.
    return tarred_files


def extract_tarred_read_set(tar_file):
    tar_command = "tar -xf %s" % tar_file
    tar_proc = subprocess.Popen(tar_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()
    print("Output of tar command", stdout, stderr)


def run_albacore(tarred_read_set):
    folder = tarred_read_set.replace(".tar.gz", "/")
    qsub_log_file = QSUB_LOG_DIR + folder.split("/")[-2] + ".o.log"
    qsub_error_file = QSUB_LOG_DIR + folder.split("/")[-2] + ".e.log"

    # This command is one giant line that we're going to break down into components.
    # The docker command opens up the docker, of Ubuntu 16.04
    docker_command = "sudo docker run --rm -v $(pwd):$(pwd) " \
                     "-u $(id -u):$(id -g) -w $(pwd) " \
                     "genomicpariscentre/albacore"
    # The read_fast5_basecaller is the algorithm that does the actual basecalling,
    # what would be run if we just had Ubuntu.
    basecaller_command = "read_fast5_basecaller.py " \
                         "--input %s " \
                         "--worker_threads %s " \
                         "--save_path %s " \
                         "--config %s" \
                         % (folder, NUM_THREADS, OUTPUT_DIR, CHOSEN_CONFIG)

    # These are both parsed into qsub which then determines what to do with it all.
    qsub_command = "qsub -o %s -e %s -S /bin/bash" % (qsub_log_file, qsub_error_file)

    # Put these all together into one grand command
    albacore_command = "echo \"%s %s\" | %s " % (docker_command, basecaller_command, qsub_command)

    # Execute via subprocess.
    albacore_proc = subprocess.Popen(albacore_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = albacore_proc.communicate()
    print("Output of albacore_proc", stdout, stderr)


def is_still_transferring():
    # Get list of files in dest directory.
    transferring_file = [t_file for t_file in os.listdir(PARENT_DIRECTORY)
                         if t_file == "TRANSFERRING"]  # Lock name is called transferring.
    # Can the file be found.
    if len(transferring_file) == 0:
        return False
    else:
        return True


def take_a_break():
    time.sleep(60)


main()
