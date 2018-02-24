#!/usr/bin/env python3

"""
This albacore script searches for tar zipped folders in a directory.
From here, it extracts the fast5 files and runs them through the albacore basecaller.

We can scale up the threads with qsub such that each of these tarred files can be run in
parallel.

Each folder takes about 30 minutes to run on 10 threads. So with infinite threads and ram
we could basecall everything
in half an hour!

The transfer_fast5_to_server.py script places a "TRANSFERRING" lock file on the folder.
This allows the transfer_fast5_to_server.py script to know when to stop looking.

"""

# Import necessary modules
import os  # General directory list stuff
import subprocess  # Executing each of the qsub commands.
import sys  # For errors
import time  # For schnoozing!
import pandas as pd  # For the status.csv file
import logging  # For logging
import fileinput  # Manipulating sge files
from datetime import datetime  # Logging times of actions
from pathlib import Path  # Creating lock files
import shutil  # Deleting opened directories
from Bio import SeqIO
import gzip
import statistics
import argparse

# Before we begin, are we using python 3.6 or greater?
try:
    assert sys.version_info >= (3, 6)
    my_python_3p6_string = f"If you see a syntax error here you need to update your python version"
except AssertionError:
    sys.exit("Error: Python version out of date. Require 3.6 or higher.")

STATUS_STANDARD_COLUMNS = ['name',
                           'albacore_submitted', 'albacore_jobid',
                           'albacore_commenced', 'albacore_complete',
                           ]
STATUS_DF = pd.DataFrame(columns=STATUS_STANDARD_COLUMNS)

CONFIGS = {"FC106_RAD001": "r94_250bps_linear.cfg",  # Rapid sequencing
           "FC106_LSK208_2d": "r94_250bps_2d.cfg",   # 2D unsure which one
           "FC106_LSK108": "r94_450bps_linear.cfg",  # For 1D ligation sequencing.
           "FC106_RAD002": "r94_450bps_linear.cfg"}  # Second Rapid sequencing kit.

FLOWCELLS = ["FLO-MIN107", "FLO-MIN106"]

KITS = ["SQK-LWP001", "SQK-NSK007", "VSK-VBK001", "SQK-RAS201", "SQK-RBK001", "SQK-LWB001",
        "SQK-RNA001", "SQK-RLI001", "SQK-RAD002", "SQK-RLB001", "SQK-RAB201", "SQK-LSK208", 
        "SQK-LSK108", "SQK-RAD003", "SQK-DCS108", "SQK-PCS108", "SQK-LSK308"]  # More kits to come


class Subfolder:
    """
    Equivalent of the subfolder class from the starter.
    This one finds the tar gzipped fast5 files in the fast5 folder
    """
    def __init__(self, name, main_dir):
        self.name = name
        # Get the relative name of the tar file
        self.tar_filename = name + ".fast5.tar.gz"
        # Set personal directories up
        self.reads_dir = os.path.join(main_dir, "fast5", name)
        self.albacore_summary_file = self.name + ".sequencing_summary.txt"
        self.albacore_log_file = self.name + ".pipeline.log"
        self.fastq_file = name + ".fastq.gz"
        self.albacore_submission_file = self.name + ".albacore.batch.sh"
        self.albacore_qsub_output_log = self.name + ".albacore.o.log"
        self.albacore_qsub_error_log = self.name + ".albacore.e.log"
        self.albacore_jobid = -1
        # Set dataframes up
        self.raw_dataframe = name + ".tsv"
        self.merged_dataframe = name + ".merged.tsv"
        # Status related attributes
        self.albacore_submitted = False
        self.albacore_commenced = False
        self.albacore_complete = False
        self.metadata_merged = False

    def to_series(self):
        return pd.Series(data=[self.name,
                               str(self.albacore_submitted), self.albacore_jobid,
                               str(self.albacore_commenced), str(self.albacore_complete)
                               ],
                         index=STATUS_STANDARD_COLUMNS
                         )

    def check_albacore_job_status(self, check_failed=True):
        if has_commenced(self.albacore_jobid):
            self.albacore_commenced = True
        if has_completed(self.albacore_jobid, check_failed):
            self.albacore_complete = True


"""
Main functions:
1. main
2. run_pipeline()
"""


def main(args):
    # Basic house cleaning, get arguments, make sure they're legit.
    args = get_args()
    dir_dict, configurations = check_directories(args)
    # Make sure there's no lockfile
    if os.path.isfile(configurations['lockfile']):
        print("Lock file exists. Exiting")
        return
    # Create directories
    dir_dict = create_directories(dir_dict)
    set_logger(configurations['logger_path'])
    create_basecalling_lock_file(configurations['lockfile'])

    # Get all the subfolders currently in the fast5 directory
    subfolders = get_subfolders(dir_dict=dir_dict)

    # Make sure STATUS_DF has all of the subfolders present
    status_df = generate_dataframe(subfolders)

    # Otherwise run through the pipeline
    run_pipeline(subfolders, configurations["max_processes"], dir_dict["metadata.merged"],
                 dir_dict, configurations)

    # Now wait for albacore to finish
    while is_still_basecalling(subfolders):

        # Otherwise run through the pipeline
        run_pipeline(subfolders, configurations["max_processes"], dir_dict["metadata.merged"],
                     dir_dict, configurations)

        # Generate dataframe throughout each iteration of the pipeline
        status_df = generate_dataframe(subfolders, status_df)

        # Have a 15 second break between iterations
        take_a_break()

    # Merge fastq files at the end of the run.
    generate_dataframe(subfolders, status_df)
    remove_basecalling_lock_file(configurations['lockfile'])


def run_pipeline(subfolders, processes, merged_dir, dir_dict, configurations):
    # Now run albacore
    index = sum([subfolder
                 for subfolder in subfolders
                    if subfolder.albacore_commenced
                    and not subfolder.albacore_complete])

    for subfolder in subfolders:
        if index > processes:
            break
        if not subfolder.albacore_submitted:
            run_albacore(subfolder, dir_dict, configurations)
            index += 1

    # Check for complete albacore jobs
    [subfolder.check_albacore_job_status()
     for subfolder in subfolders
        if subfolder.albacore_submitted
        and not subfolder.albacore_complete]

    # If albacore is complete. Merge the fastq output onto the metadata frame
    [merge_dataframes(subfolder, merged_dir)
     for subfolder in subfolders
        if subfolder.albacore_complete and
        not subfolder.metadata_merged]


"""
General script initialisation functions
2. check_directories
3. set_logger
4. pick_up_from_previous_run
"""


def get_args():
    # Albacore server arguments:
    albacore_parser = argparse.ArgumentParser(dscription="The albacore-server scaled command incorporates qsub to spread " +
                                                         "the server load and rapidly " + "generate data off the PromethION.")
    albacore_parser.add_argument("--fast5_dir", type=str, required=True,
                                 help="/path/to/reads/fast5, " +
                                      "should have a bunch of tar zipped files in it.")
    albacore_parser.add_argument("--flowcell", type=str, choices=FLOWCELLS,
                                 help="Pick a flowcell version")
    albacore_parser.add_argument("--kit", type=str, choices=KITS,
                                 help="Pick a kit type")
    albacore_parser.add_argument("--output_dir", type=str, required=False, default=None,
                                 help="Will be called 'albacore'" +
                                      " and sit adjacent to the reads folder if left blank.")
    albacore_parser.add_argument("--num_threads", type=int, required=False, default=5,
                                 help="How many threads did you wish to use " +
                                      "per parallel output default is 5.")
    albacore_parser.add_argument("--num_processes", type=int, required=False, default=10,
                                 help="How many processes did you wish to process at once?")
    albacore_parser.add_argument("--fastq_dir", type=str, required=False, default=None,
                                 help="Where should the fastq data be placed." +
                                      "Will be called 'fastq'" +
                                      "and sit adjacent to reads folder if left blank.")
    albacore_parser.add_argument("--qsub_directory", type=str, required=False, default=None,
                                 help="Where would you like to place the qsub files?")
    albacore_parser.add_argument("--qsub_host", type=str, required=False, default=None,
                                 help="Where would you like the qsub jobs to be run?")
    albacore_parser.add_argument("--max_processes", type=int, required=False, default=None,
                                 help="Limit the number of jobs that can be processed at any given moment." +
                                      "This command will prevent extraction/albacore jobs from being submitted while" +
                                      "there exists 'max_processes' jobs running / in the queue.")
    albacore_parser.add_argument("--log_directory", type=str, required=False, default=None,
                                 help="Where do you wish to store the poreduck logs? " +
                                      "Will be called 'poreduck_logs' and sit adjacent to reads folder if left blank")
    albacore_parser.add_argument("--albacore_version", type=str, required=True,
                                 help="Albacore 2 comes with the pass and fail folders." +
                                      "Albacore 1 does not. Also used in %VER% command in templates." +
                                      "Expressed as <majore_version>.<minor_version>.<patch>")
    albacore_parser.add_argument("--qsub_albacore_template", type=str, required=True,
                                 help="The qsub albacore template for your qsub command. " +
                                      "Check the poreduck examples for more information.")
    return albacore_parser.parse_args()


def check_directories(args):
    # Make sure the directories exist, change to reads directory,
    # Create any other necessary directories for the script to run.
    dir_dict = {}
    configurations = {}
    # Get arguments used throughout the script
    configurations['threads'] = args.num_threads
    configurations['processes'] = args.max_processes
    configurations['flowcell'] = args.flowcell
    configurations['kit'] = args.kit
    configurations['template'] = args.qsub_albacore_template
    configurations["albacore.version"] = args.albacore_version

    dir_dict["fast5"] = os.path.abspath(args.fast5_dir)
    if not os.path.isdir(args.reads_dir):
        sys.exit(f"Error, {args.reads_dir} does not exist")
    os.chdir(dir_dict["reads"])

    dir_dict["main"] = os.path.abspath(os.path.join(dir_dict["fast5"], os.pardir))
    dir_dict["metadata"] = os.path.join(dir_dict["main"], "metadata")

    if args.qsub_dir is None:
        dir_dict["qsub"] = os.path.join(dir_dict["main"], "qsub")
    else:
        dir_dict["qsub"] = os.path.abspath(args.qsub_dir)

    if args.fastq_dir is None:
        dir_dict["fastq"] = os.path.join(dir_dict["main"], "fastq")
    else:
        dir_dict["fastq"] = os.path.abspath(dir_dict["fastq"])


    if dir_dict["log"] == "":
        dir_dict["log"] = os.path.join(dir_dict["main"], "poreduck_logs")
    else:
        dir_dict["log"] = os.path.abspath(dir_dict["log"])

    configurations['logger_path'] = os.path.join(dir_dict["log"], '.'.join([str(datetime.now().date()),
                                                 str(datetime.now().time()).replace(":", "-"),
                                                 "albacore_pipeline.log"]))
    configurations['lock_file'] = os.path.join(dir_dict['main_dir'], "BASECALLING")

    return dir_dict, configurations


def set_logger(logger_path):
    global logger
    logging.basicConfig(filename=logger_path, level=logging.INFO,
                        format='%(asctime)s::\t%(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
    logger = logging.getLogger()


def create_directories(dir_dict):
    for subdir in dir_dict:
        if not os.path.isdir(subdir):
            os.mkdir(subdir)
    dir_dict["metadata.merged"] = os.path.join(dir_dict['metadata'], "merged")
    return dir_dict

"""
Pipeline functions:
1. extract_tarred_read_set
2. run_albacore
3. remove_folder
4. move_fastq_file
5. tar_albacore_folder
6. merge_fastq_files_wrapper
"""


def run_albacore(subfolder, dir_dict, configurations):
    """ Generate read_fast5_basecaller command
        And qsub command and pipe former into latter
        This function also sets the jobid of a given folder.
    """
    if subfolder.albacore_submitted:
        return   # Basecalling has already queued for this directory

    # Commence basecalling on this directory
    subfolder.albacore_submitted = True

    # The sbatch is all primed.
    # The following variables need to be set for sbatch.
    slurm_options = {"mem":"%d" % configurations['threads']*3,
                     "error":"albacore.{}.{}".format(subfolder.name, "%j"),
                     "chdir": dir_dict["main_dir"]
                     "job-name": "albacore"}
    # The following variables need to be exported
    export_options = {"NUM_THREADS": configurations["threads"],
                      "FLOWCELL": configurations["flowcell"],
                      "KIT": configurations["kit"],
                      "SUBFOLDER_NAME": subfolder.name,
                      "ALBACORE_DIR": dir_dict["albacore_dir"],
                      "ALBACORE_VER": configurations["albacore.version"]}

    job_submission_command = ["sbatch", configurations['template']]
    job_submission_command.extend(["--%s=%s" % (k,v) for k,v in slurm_options])
    job_submission_command.append("--export=%s" % ','.join(["%s=%s" % (k,v) for k,v in export_options]))

    job_submission_proc = subprocess.run(job_submission_command,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = job_submission_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of th   e array.
    logger.info(f"Output of albacore sge submission \nStdout:\"{stdout.rstrip()}\"\nStderr:\"{stderr.rstrip()}\"")
    # Get job ID:
    job_id = stdout.rstrip().split(";")[0]
    if job_id.isdigit():
        subfolder.albacore_jobid = int(job_id)
    else:
        logger.error(f"Error: Could not assign job id from {stdout.rstrip()}")
        sys.exit(f"Error: Could not assign job id from {stdout.rstrip()}")


def get_fastq_dataframe(fastq_file):
    fastq_df = pd.DataFrame(columns=["fastq_id", "read", "channel", "seq_length", "av_qual"])
    with gzip.open(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Get metadata
            fastq_id = record.id.split()[0]
            row_as_dict = dict(x.split("=") for x in record.description.split()[1:])
            # Get length
            fastq_length = len(record.seq)
            fastq_quality = statistics.mean(record.letter_annotations["phred_quality"])
            fastq_df = fastq_df.append(pd.Series([fastq_id, row_as_dict["read"], row_as_dict["ch"],
                                                  fastq_length, fastq_quality]),
                                       ignore_index=True)
    # Specify types for each line
    numeric_cols = ["read", "channel", "seq_length", "av_qual"]
    fastq_df[numeric_cols] = fastq_df[numeric_cols].apply(pd.to_numeric)
    fastq_df.set_index(["channel", "read"], inplace=True)
    return fastq_df


def merge_dataframes(subfolder, merged_dir):
    # Use the SeqIO library to extract the data info from the fastq file
    pass_df = get_fastq_dataframe(subfolder.pass_fastq_file)
    pass_df["Class"] = "pass"
    fail_df = get_fastq_dataframe(subfolder.fail_fastq_file)
    pass_df["Class"] = "fail"
    fastq_df = pd.concat(fail_df, pass_df, axis='rows')
    # Now merge with the raw metadata file
    raw_df = pd.read_csv(subfolder.metadata_raw, header=True)
    raw_df.set_index(["channel", "read"], inplace=True)
    pd.concat(raw_df, fastq_df, join='right').to_csv(merged_dir, subfolder.merged_dataframe))
    

"""
Miscellaneous pipeline non-core functions
1. get_subfolders
2. is_still_transferring
3. take_a_break
4. generate_dataframe
5. new_directories
6. is_still_basecalling
7. get_current_subfolder
8. num_current_jobs
9. create_basecalling_lock_file
10. remove_basecalling_lock_file
"""


def get_subfolders(subfolders=[], dir_dict={}):
    """ New subfolders will emerge to the directory as tarballs.
        Get a list of new tarballs then append them to the SUBFOLDERS
        list
    """
    # Don't reset existing Subfolder classes
    existing_tars = [subfolder.tar_filename
                     for subfolder in subfolders]

    # Or those that have been basecalled in a previous script.
    metadata_tar = [metadata.replace(".merged.tsv")
                    for metadata in os.listdir(dir_dict["metadata.merged"])]

    # Get a list of tarballs
    new_subfolders = sorted([Subfolder(tarred_folder.replace(".fast5.tar.gz", ""))
                             for tarred_folder in os.listdir(dir_dict["reads"])
                             if tarred_folder.endswith(".fast5.tar.gz")
                             and not tarred_folder in existing_tars
                             and not tarred_folder.replace(".fast5.tar.gz")
                             in metadata_tar])
    subfolders.extend(new_subfolders)
    return subfolders


def take_a_break():
    """
    Just a 15 second sleep ;)
    """
    time.sleep(15)


def generate_dataframe(subfolders, status_df=None):

    for subfolder in subfolders:
        if status_df is None or subfolder.name not in status_df.index.tolist():
            # Append row
            status_df = status_df.append(subfolder.to_series(), ignore_index=True)
            # Reset index
            status_df.set_index('name', drop=False, inplace=True)

    status_df.to_csv(status_df, index=False)
    return status_df


def is_still_basecalling(subfolders):
    """Any subfolders left to be basecalled?"""
    if len([subfolder for subfolder in subfolders
            if not subfolder.metadata_merged]) == 0:
        return False  # basecalling finished
    else:  # Basecalling still going
        return True

"""
qsub specific functions
1. has_commenced(job_id)
2. has_completed(job_id)
3. has_failed(job_id)
"""


def has_commenced(job_id):
    """
    Use the qacct command to see if there exists a start time
    No start_time is represented as -/-
    """
    qacct_command = f"sacct -j {job_id} | tail -n 1 | grep RUNNING | wc -l"
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    qacct_stdout = qacct_stdout.decode()
    qacct_stderr = qacct_stderr.decode()
    if int(qacct_stdout) > 1:
        # Start-time is non-existent
        return False
    if not qacct_stderr == "":
        logger.info(f"qacct stderr of {qacct_stderr}")
    return True


def has_completed(job_id, check_failed):
    """
    Use the qacct command to see if the job has finished.
    No end_time is represented as -/-
    """
    qacct_command = f"sacct -j {job_id} | grep 'COMPLETED\|FAILED'"
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    qacct_stdout = qacct_stdout.decode()
    qacct_stderr = qacct_stderr.decode()
    if qacct_stdout == "":
        # end_time is non-existent
        return False
    if not qacct_stderr == "":
        logger.info(f"qacct stderr of {qacct_stderr}")
    
    if check_failed and has_failed(job_id):
        logger.info(f"It appears that the job {job_id} has failed")
        sys.exit(f"Failing because {job_id} failed.")

    return True


def has_failed(job_id):
    """
    Use the qacct command to see if the job has failed.
    We pass the exit_status parameter into awk and sum it.
    If it's any greater than zero then the command has failed
    """

    qacct_command = f"sacct {job_id} | tail -n1 | grep FAILED | wc -l"
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    qacct_stdout = qacct_stdout.decode()
    qacct_stderr = qacct_stderr.decode()
    if qacct_stdout.strip() == "":
        return True  
    if int(qacct_stdout) > 0:
        logger.info(qacct_stdout)
        # Job has failed
        return True
    if not qacct_stderr == "":
        logger.info(f"qacct stderr of {qacct_stderr}")
    return False


def create_basecalling_lock_file(lockfile):
    # Create a file called BASECALLING in the working directory
    Path(lockfile).touch()


def remove_basecalling_lock_file(lockfile):
    # Remove the basecalling lock file.
    os.remove(lockfile)
