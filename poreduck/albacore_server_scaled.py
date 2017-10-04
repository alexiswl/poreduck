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
import argparse  # For allowing users to configure the arguments used.
import sys  # For errors
import time  # For schnoozing!
import pandas as pd  # For the status.csv file
import logging
import fileinput
from datetime import datetime
from pathlib import Path

# Before we begin, are we using python 3.6 or greater?
try:
    assert sys.version_info >= (3, 6)
    my_python_3p6_string = f"If you see a syntax error here you need to update your python version"
except AssertionError:
    sys.exit("Error: Python version out of date. Require 3.6 or higher.")

# Set global and semi-global variables
TRANSFERRING = True
BASECALLING_LOCK_FILE = None
READS_DIR = ""
CW_DIR = ""
ALBACORE_DIR = ""
WORKING_DIR = ""
NUM_THREADS = 0
MAX_PROCESSES = 0
CHOSEN_FLOWCELL = ""
CHOSEN_KIT = ""
PARENT_DIRECTORY = ""
QSUB_LOG_DIR = ""
FASTQ_DIR = ""
SUBFOLDERS = []
STATUS_CSV = ""
QSUB_HOST = None
QSUB_TYPE = ""
QSUB_EXTRACTION_TEMPLATE = ""
QSUB_ALBACORE_TEMPLATE = ""
STATUS_STANDARD_COLUMNS = ['name', 'extracted_submitted', 'extracted_jobid',
                           'extracted_commenced', 'extracted_complete',
                           'albacore_submitted', 'albacore_jobid',
                           'albacore_commenced', 'albacore_complete',
                           'folder_removed', 'fastq_moved']
STATUS_DF = pd.DataFrame(columns=STATUS_STANDARD_COLUMNS)

CONFIGS = {"FC106_RAD001": "r94_250bps_linear.cfg",  # Rapid sequencing
           "FC106_LSK208_2d": "r94_250bps_2d.cfg",   # 2D unsure which one
           "FC106_LSK108": "r94_450bps_linear.cfg",  # For 1D ligation sequencing.
           "FC106_RAD002": "r94_450bps_linear.cfg"}  # Second Rapid sequencing kit.
FLOWCELLS = ["FLO-MIN107", "FLO-MIN106"]
KITS = ["SQK-LWP001", "SQK-NSK007", "VSK-VBK001", "SQK-RAS201", "SQK-RBK001", "SQK-LWB001",
        "SQK-RNA001", "SQK-RLI001", "SQK-RAD002", "SQK-RLB001", "SQK-RAB201", "SQK-LSK208", 
        "SQK-LSK108", "SQK-RAD003", "SQK-DCS108", "SQK-PCS108", "SQK-LSK308"]  # More kits to come
QSUB_TYPES = ["SGE", "TORQUE"]
LOGGER = None
LOGGER_DIR = ""
LOGGER_PATH = ""
BARCODING = False


class Subfolder:
    def __init__(self, name):
        self.name = name
        # Get the relative name of the tar file
        self.tar_filename = name + ".tar.gz"
        # Set personal directories up
        self.reads_dir = os.path.join(READS_DIR, name)
        self.albacore_dir = os.path.join(ALBACORE_DIR, name)
        # Albacore related files
        if CHOSEN_KIT == "SQK-LSK308":
            self.workspace_dir = os.path.join(self.albacore_dir, "1dsq_analysis", "workspace")
        else:
            self.workspace_dir = os.path.join(self.albacore_dir, "workspace")
        self.albacore_summary_file = os.path.join(self.albacore_dir, "sequencing_summary.txt")
        self.albacore_log_file = os.path.join(self.albacore_dir, "pipeline.log")
        self.fastq_file = name + ".fastq"
        self.albacore_zip_file = self.albacore_dir + ".albacore.tar.gz"
        # Qsub related files / ids
        self.extracted_submission_file = os.path.join(QSUB_LOG_DIR, name + ".extract.batch.sh")
        self.extracted_qsub_output_log = os.path.join(QSUB_LOG_DIR, name + ".extract.o.log")
        self.extracted_qsub_error_log = os.path.join(QSUB_LOG_DIR, name + ".extract.e.log")
        self.extracted_jobid = -1
        self.albacore_submission_file = os.path.join(QSUB_LOG_DIR, name + ".albacore.batch.sh")
        self.albacore_qsub_output_log = os.path.join(QSUB_LOG_DIR, name + ".albacore.o.log")
        self.albacore_qsub_error_log = os.path.join(QSUB_LOG_DIR, name + ".albacore.e.log")
        self.albacore_jobid = -1
        # Status related attributes
        self.extracted_submitted = False
        self.extracted_commenced = False
        self.extracted_complete = False
        self.albacore_submitted = False
        self.albacore_commenced = False
        self.albacore_complete = False
        self.albacore_tarred = False
        self.folder_removed = False
        self.fastq_moved = False

    def to_series(self):
        return pd.Series(data=[self.name,
                               str(self.extracted_submitted), self.extracted_jobid,
                               str(self.extracted_commenced), str(self.extracted_complete),
                               str(self.albacore_submitted), self.albacore_jobid,
                               str(self.albacore_commenced), str(self.albacore_complete),
                               str(self.folder_removed), str(self.fastq_moved)],
                         index=STATUS_STANDARD_COLUMNS
                         )

    def check_albacore_job_status(self, check_failed=True):
        if has_commenced(self.albacore_jobid):
            self.albacore_commenced = True
            update_dataframe(get_current_subfolder(self.name))
        if has_completed(self.albacore_jobid, check_failed):
            self.albacore_complete = True
            update_dataframe(get_current_subfolder(self.name))

    def check_extraction_job_status(self, check_failed=True):
        if has_commenced(self.extracted_jobid):
            self.extracted_commenced = True
            update_dataframe(get_current_subfolder(self.name))
        if has_completed(self.extracted_jobid, check_failed):
            self.extracted_complete = True
            update_dataframe(get_current_subfolder(self.name))


"""
Main functions:
1. main
2. run_pipeline()
"""


def main(args):
    global TRANSFERRING
    # Basic house cleaning, get arguments, make sure they're legit.
    set_global_variables(args)
    check_directories()
    set_logger()
    create_basecalling_lock_file()

    # Pick up previous basecalling
    pick_up_from_previous_run()

    # While the transferring lock script exists.
    while TRANSFERRING:
        # Check that we're still transferring
        if not is_still_transferring():
            TRANSFERRING = False

        # Get all the subfolders currently in the fast5 directory
        get_subfolders()

        # Make sure STATUS_DF has all of the subfolders present
        generate_dataframe()

        # Otherwise run through the pipeline
        run_pipeline()

        # Have a break
        take_a_break()

    # Transferring from laptop complete just wait for albacore to finish.
    processing = True  # Set while command
    while processing:
        # Can we make this the last time around?
        if not is_still_basecalling():
            processing = False

        # Otherwise run through the pipeline
        run_pipeline()

        # Generate dataframe throughout each iteration of the pipeline
        generate_dataframe()

        # Have a 15 second break between iterations
        take_a_break()

    # Merge fastq files at the end of the run.
    generate_dataframe()    
    #merge_fastq_files_wrapper()
    remove_basecalling_lock_file()

def run_pipeline():
    # Extract tarballs
    [extract_tarred_read_set(subfolder) for subfolder in SUBFOLDERS
     if not subfolder.extracted_submitted]

    # Check for complete extraction jobs
    [subfolder.check_extraction_job_status() for subfolder in SUBFOLDERS
     if subfolder.extracted_submitted and
     not subfolder.extracted_complete]

    # Now run albacore
    [run_albacore(subfolder) for subfolder in SUBFOLDERS
     if subfolder.extracted_complete and
     not subfolder.albacore_submitted]

    # Check for complete albacore jobs
    [subfolder.check_albacore_job_status() for subfolder in SUBFOLDERS
     if subfolder.albacore_submitted and
     not subfolder.albacore_complete]

    # Once albacore is complete on a subfolder
    # Remove folder from existence, leaving just the tarball again.
    [remove_folder(subfolder) for subfolder in SUBFOLDERS
     if subfolder.albacore_complete and
     not subfolder.folder_removed]

    # Move the fastq file to the FASTQ_DIR
    [move_fastq_file(subfolder) for subfolder in SUBFOLDERS
     if subfolder.albacore_complete and
     not subfolder.fastq_moved]

    # Tar albacore directory
    [tar_albacore_folder(subfolder) for subfolder in SUBFOLDERS
     if subfolder.albacore_complete and
     not subfolder.albacore_tarred]


"""
General script initialisation functions
1. get_arguments
2. set_global_variables
3. check_directories
4. set_logger
4. pick_up_from_previous_run
"""


def get_arguments():
    parser = argparse.ArgumentParser(
        description="The albacore-server scaled command incorporates qsub to spread " +
                    "the server load and rapidly " + "generate data off the MinION.")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="/path/to/reads, " +
                             "should have a bunch of tar zipped files in it.")
    parser.add_argument("--flowcell", type=str, choices=FLOWCELLS,
                        help="Pick a flowcell version")
    parser.add_argument("--kit", type=str, choices=KITS,
                        help="Pick a kit type")
    parser.add_argument("--output_dir", type=str, required=False, default=None,
                        help="Will be called 'albacore'" +
                             " and sit adjacent to the reads folder if left blank.")
    parser.add_argument("--num_threads", type=int, required=False, default=5,
                        help="How many threads did you wish to use " +
                             "per parallel output default is 5.")
    parser.add_argument("--fastq_dir", type=str, required=False, default=None,
                        help="Where should the fastq data be placed." +
                             "Will be called 'fastq'" +
                             "and sit adjacent to reads folder if left blank.")
    parser.add_argument("--resume", type=str, required=False, default=None,
                        help="Resume the albacore run, need a csv file")
    parser.add_argument("--qsub_directory", type=str, required=False, default=None,
                        help="Where would you like to place the qsub files?")
    parser.add_argument("--qsub_host", type=str, required=False, default=None,
                        help="Where would you like the qsub jobs to be run?")
    parser.add_argument("--max_processes", type=int, required=False, default=None,
                        help="Limit the number of jobs that can be processed at any given moment." +
                             "This command will prevent extraction/albacore jobs from being submitted while" +
                             "there exists 'max_processes' jobs running / in the queue.")
    parser.add_argument("--log_directory", type=str, required=False, default=None,
                        help="Where do you wish to store the poreduck logs? " +
                             "Will be called 'poreduck_logs' and sit adjacent to reads folder if left blank")
    parser.add_argument("--barcoding", default=False, dest='barcoding', action='store_true',
                        help="Use this option to demultiplex library?")
    parser.add_argument("--qsub_type", choices=QSUB_TYPES, default="SGE",
                        help="What qsub system are you using?")
    parser.add_argument("--qsub_extraction_template", type=str, required=True,
                        help="The qsub extraction template for your qsub command. " +
                             "Check the poreduck examples for more information.")
    parser.add_argument("--qsub_albacore_template", type=str, required=True,
                        help="The qsub albacore template for your qsub command. " +
                             "Check the poreduck examples for more information.")

    return parser.parse_args()


def set_global_variables(args):
    # Global variables
    global READS_DIR, ALBACORE_DIR, WORKING_DIR, NUM_THREADS, CHOSEN_KIT, FASTQ_DIR
    global STATUS_CSV, QSUB_LOG_DIR, CW_DIR, QSUB_HOST, CHOSEN_FLOWCELL, MAX_PROCESSES
    global LOGGER_DIR, BARCODING, QSUB_TYPE, QSUB_ALBACORE_TEMPLATE, QSUB_EXTRACTION_TEMPLATE
    READS_DIR = args.reads_dir
    if args.output_dir is not None:
        ALBACORE_DIR = args.output_dir
    CHOSEN_KIT = args.kit
    CHOSEN_FLOWCELL = args.flowcell
    NUM_THREADS = args.num_threads
    if args.fastq_dir is not None:
        FASTQ_DIR = args.fastq_dir
    if args.resume is not None:
        STATUS_CSV = args.resume
    if args.qsub_directory is not None:
        QSUB_LOG_DIR = args.qsub_directory
    CW_DIR = os.getcwd()
    if args.qsub_host is not None:
        QSUB_HOST = args.qsub_host
    if args.max_processes is not None:
        MAX_PROCESSES = args.max_processes
    if args.log_directory is not None:
        LOGGER_DIR = args.log_directory
    BARCODING = args.barcoding
    QSUB_EXTRACTION_TEMPLATE = os.path.abspath(args.qsub_extraction_template)
    QSUB_ALBACORE_TEMPLATE = os.path.abspath(args.qsub_albacore_template)
    QSUB_TYPE = args.qsub_type


def check_directories():
    # Make sure the directories exist, change to reads directory,
    # Create any other necessary directories for the script to run.
    global READS_DIR, ALBACORE_DIR, PARENT_DIRECTORY, QSUB_LOG_DIR, FASTQ_DIR
    global STATUS_CSV, LOGGER_DIR, LOGGER_PATH, QSUB_TYPE, BASECALLING_LOCK_FILE
    if not os.path.isdir(READS_DIR):
        sys.exit(f"Error, {READS_DIR} does not exist")

    READS_DIR = os.path.abspath(READS_DIR)
    os.chdir(READS_DIR)

    PARENT_DIRECTORY = os.path.abspath(os.path.join(READS_DIR, os.pardir)) + "/"

    if ALBACORE_DIR == "":
        ALBACORE_DIR = os.path.join(PARENT_DIRECTORY, "albacore")
    else:
        ALBACORE_DIR = os.path.abspath(ALBACORE_DIR)

    if not os.path.isdir(ALBACORE_DIR):
        os.mkdir(ALBACORE_DIR)

    QSUB_LOG_DIR = os.path.join(PARENT_DIRECTORY, "qsub_log")
    if not os.path.isdir(QSUB_LOG_DIR):
        os.mkdir(QSUB_LOG_DIR)

    if FASTQ_DIR == "":
        FASTQ_DIR = os.path.join(PARENT_DIRECTORY, "fastq")
    else:
        FASTQ_DIR = os.path.abspath(FASTQ_DIR)

    if not os.path.isdir(FASTQ_DIR):
        os.mkdir(FASTQ_DIR)
    if LOGGER_DIR == "":
        LOGGER_DIR = os.path.join(PARENT_DIRECTORY, "poreduck_logs")
    else:
        LOGGER_DIR = os.path.abspath(LOGGER_DIR)
    if not os.path.isdir(LOGGER_DIR):
        os.mkdir(LOGGER_DIR)
    LOGGER_PATH = os.path.join(LOGGER_DIR, '.'.join([str(datetime.now().date()),
                                                     str(datetime.now().time()),
                                                     "albacore_pipeline.log"]))

    if not STATUS_CSV == "":
        if not os.path.isfile(os.path.join(CW_DIR, STATUS_CSV)):
            sys.exit("Resume csv specified but does not exist")
    else:
        STATUS_CSV = os.path.join(PARENT_DIRECTORY, "status.csv")
        if os.path.isfile(STATUS_CSV):
            print(f"Warning, STATUS_CSV not defined but still exists. Removing file {STATUS_CSV}")
            os.remove(STATUS_CSV)
    BASECALLING_LOCK_FILE = os.path.join(PARENT_DIRECTORY, "BASECALLING")


def set_logger():
    global LOGGER
    logging.basicConfig(filename=LOGGER_PATH, level=logging.INFO,
                        format='%(asctime)s::\t%(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
    LOGGER = logging.getLogger() 


def pick_up_from_previous_run():
    global SUBFOLDERS
    if os.path.isfile(STATUS_CSV):
        previous_dataframe = pd.read_csv(STATUS_CSV)
        LOGGER.info(f"Found {STATUS_CSV}. Putting run back on track.")

        # Set subfolders and all the statuses
        for index, row in previous_dataframe.iterrows():
            subfolder = row['name']
            SUBFOLDERS.append(Subfolder(subfolder)) 
            SUBFOLDERS[-1].extracted_submitted = row['extracted_submitted']
            SUBFOLDERS[-1].extracted_commenced = row['extracted_commenced']
            SUBFOLDERS[-1].extracted_jobid = row['extracted_jobid']
            SUBFOLDERS[-1].extracted_complete = row['extracted_complete']
            SUBFOLDERS[-1].albacore_submitted = row['albacore_submitted']
            SUBFOLDERS[-1].albacore_commenced = row['albacore_commenced']
            SUBFOLDERS[-1].albacore_jobid = row['albacore_jobid']
            SUBFOLDERS[-1].albacore_complete = row['albacore_complete']
            SUBFOLDERS[-1].folder_removed = row['folder_removed']
            SUBFOLDERS[-1].fastq_moved = row['fastq_moved']
            
            # Check to see if the jobs may have finished after the run crashed
            SUBFOLDERS[-1].check_extraction_job_status(check_failed=False) 
            SUBFOLDERS[-1].check_albacore_job_status(check_failed=False)
              
            # Now make sure none of the jobs have failed and reset if so
            if has_failed(SUBFOLDERS[-1].albacore_jobid):
                LOGGER.info(f"It appears {SUBFOLDERS[-1].name} has failed albacore processing. " +
                            "Reverting albacore submitted, commenced to false. " +
                            f"Removing job id. {SUBFOLDERS[-1].albacore_jobid}")
                SUBFOLDERS[-1].albacore_submitted = False
                SUBFOLDERS[-1].albacore_commenced = False
                SUBFOLDERS[-1].albacore_jobid = -1
                SUBFOLDERS[-1].albacore_complete = False
            if has_failed(SUBFOLDERS[-1].extracted_jobid):
                LOGGER.info(f"It appears {SUBFOLDERS[-1].name} has failed to extract. " +
                            "Reverting albacore submitted, commenced to false. " +
                            f"Removing job id. {SUBFOLDERS[-1].extracted_jobid}")
                SUBFOLDERS[-1].extracted_submitted = False
                SUBFOLDERS[-1].extracted_commenced = False
                SUBFOLDERS[-1].extracted_jobid = -1
                SUBFOLDERS[-1].extracted_complete = False


"""
Pipeline functions:
1. extract_tarred_read_set
2. run_albacore
3. remove_folder
4. move_fastq_file
5. tar_albacore_folder
6. merge_fastq_files_wrapper
"""


def extract_tarred_read_set(subfolder):
    # Extract the tarred reads, we currently do not delete after extraction...
    # maybe a good idea?
    # Has the dataset already been set for extraction?
    if subfolder.extracted_submitted:
        return
    elif num_current_jobs() >= MAX_PROCESSES and not MAX_PROCESSES == 0:
        return
    else:
        subfolder.extracted_submitted = True

    tar_command = f"pigz -dc {os.path.join(READS_DIR, subfolder.tar_filename)} | tar -xf -"

    qsub_replacement_dict = {"STDOUT": subfolder.extracted_qsub_output_log,
                             "STDERR": subfolder.extracted_qsub_error_log,
                             "HOSTNAME": QSUB_HOST, 
                             "COMMAND": tar_command,
                             "WORKING_DIRECTORY": READS_DIR}

    # Copy the standard qsub file from the main folder to the qsub folder
    cp_command = f"cp {QSUB_EXTRACTION_TEMPLATE} {subfolder.extracted_submission_file}"
    cp_proc = subprocess.Popen(cp_command, shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = cp_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    if not stdout == "" or not stderr == "":
        print("Copy command:", stdout, stderr)

    with fileinput.FileInput(subfolder.extracted_submission_file, inplace=True) as file:
        # Now edit this file based on our inputs.
        for line in file:
            for key, replacement in qsub_replacement_dict.items():
                # If the replacement is none we need to remove this line.
                if replacement is None:
                    if key in line:
                        line = ""
                    # Or continue if the key is not in this line.
                    else:
                        continue
                # Substitute in the replacement for this key in this line.
                else:
                    line = line.replace(f"%{key}%", str(replacement))
            # Print the substituted line which is captured by the fileinput.
            print(line.rstrip())

    # Submit job
    job_submission_command = f"qsub {subfolder.extracted_submission_file}"
    job_submission_proc = subprocess.Popen(job_submission_command, shell=True,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = job_submission_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of the array.
    LOGGER.info(f"Output of pigz sge submission \nStdout:\"{stdout.rstrip()}\"\nStderr:\"{stderr.rstrip()}\"")
    # Get job ID:
    if QSUB_TYPE == "SLURM":
        job_id = stdout.rstrip().split(";")[0]
    else:
        job_id = stdout.rstrip().split(".")[0]
    if job_id.isdigit():
        subfolder.extracted_jobid = int(job_id)
    else:
        LOGGER.error(f"Error: Could not assign job id from {stdout.rstrip()}")
        sys.exit(f"Error: Could not assign job id from {stdout.rstrip()}")
    update_dataframe(subfolder)


def run_albacore(subfolder):
    """ Generate read_fast5_basecaller command
        And qsub command and pipe former into latter
        This function also sets the jobid of a given folder.
    """
    if subfolder.albacore_submitted:
        return   # Basecalling has already queued for this directory
    elif num_current_jobs() >= MAX_PROCESSES and not MAX_PROCESSES == 0:
        return
    else:
        # Commence basecalling on this directory
        subfolder.albacore_submitted = True

    # Need to ramp up memory requirements if using the 1D^2 protocol
    if CHOSEN_KIT == "SQK-LSK308":
        memory_allocation = 4 + 4*NUM_THREADS
        basecaller_command_options = ["full_1dsq_basecaller.py"]
    else:
        memory_allocation = 4 + 2*NUM_THREADS
        basecaller_command_options = ["read_fast5_basecaller.py"]
    # The read_fast5_basecaller is the algorithm that does the actual base calling,
    # what would be run if we just had Ubuntu.
    basecaller_command_options.append(f"--input {subfolder.reads_dir}")
    basecaller_command_options.append(f"--worker_threads {NUM_THREADS}")
    basecaller_command_options.append(f"--save_path {subfolder.albacore_dir}")
    basecaller_command_options.append(f"--flowcell {CHOSEN_FLOWCELL}")
    basecaller_command_options.append(f"--kit {CHOSEN_KIT}")
    if BARCODING:
        basecaller_command_options.append("--barcoding")
    basecaller_command = ' '.join(basecaller_command_options)

    # These are both parsed into qsub which then determines what to do with it all.
    qsub_replacement_dict = {"STDOUT": subfolder.albacore_qsub_output_log,
                             "STDERR": subfolder.albacore_qsub_error_log,
                             "HOSTNAME": QSUB_HOST,
                             "MEM": memory_allocation, 
                             "COMMAND": basecaller_command,
                             "WORKING_DIRECTORY": ALBACORE_DIR}

    # Copy the standard qsub file from the main folder to the qsub folder
    cp_command = f"cp {QSUB_ALBACORE_TEMPLATE} {subfolder.albacore_submission_file}"
    cp_proc = subprocess.Popen(cp_command, shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = cp_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    if not stdout == "" or not stderr == "":
        print("Copy command:", stdout, stderr)

    with fileinput.FileInput(subfolder.albacore_submission_file, inplace=True) as file:
        # Now edit this file based on our inputs.
        for line in file:
            for key, replacement in qsub_replacement_dict.items():
                # If the replacement is none we need to remove this line.
                if replacement is None:
                    if key in line:
                        line = ""
                    # Or continue if the variable is not in this line.
                    else:
                        continue
                # Substitute in the replacement for this key in this line.
                else:
                    line = line.replace(f"%{key}%", str(replacement))
            # Print the substituted line which is captured by the fileinput.
            print(line.rstrip())

    # Submit job
    if QSUB_TYPE == "SGE" or QSUB_TYPE == "TORQUE":
        job_submission_command = f"qsub {subfolder.albacore_submission_file}"
    elif QSUB_TYPE == "SLURM":
        job_submission_command = f"sbatch {subfolder.albacore_submission_file}"

    job_submission_proc = subprocess.Popen(job_submission_command, shell=True,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = job_submission_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of th   e array.
    LOGGER.info(f"Output of albacore sge submission \nStdout:\"{stdout.rstrip()}\"\nStderr:\"{stderr.rstrip()}\"")
    # Get job ID:
    if QSUB_TYPE == "SLURM":
        job_id = stdout.rstrip().split(";")[0]
    else:
        job_id = stdout.rstrip().split(".")[0]
    if job_id.isdigit():
        subfolder.albacore_jobid = int(job_id)
    else:
        LOGGER.error(f"Error: Could not assign job id from {stdout.rstrip()}")
        sys.exit(f"Error: Could not assign job id from {stdout.rstrip()}")
    update_dataframe(subfolder)


def remove_folder(subfolder):
    """ Remove subfolder, leave us just with the '.tar.gz' file"""

    # Make sure that albacore has finished processing directory
    if not subfolder.albacore_complete:
        return

    # Check the folder hasn't already been removed.
    if subfolder.folder_removed:
        return
    else:  # Set this to true so that we don't try do it again!
        subfolder.folder_removed = True

    # Ensure that the tarred read set still exists
    if not os.path.isfile(os.path.join(READS_DIR, subfolder.tar_filename)):
        LOGGER.info("Um... the archive doesn't exist, I'm not going to delete the folder")
        return

    # Generate remove command
    remove_command = f"rm -rf {os.path.join(READS_DIR, subfolder.name)}"

    # Run remove command through subprocess
    remove_proc = subprocess.Popen(remove_command, shell=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = remove_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    if not stdout == "" or stderr == "":
        LOGGER.info(f"Output of deleting folder command is\nStdout:\"{stdout.rstrip()}\"\nStderr:\"{stderr.rstrip()}\"")
    update_dataframe(subfolder)


def move_fastq_file(subfolder):
    """ Get the fastq file from the albacore workspace and merge into generic fastq file
        As the fastq_file contains some weird run id info we can find using os.listdir
        to first find the fastq file in the workspace directory.
        Then we can move to the fastq folder and change name to subfolder.name + '.fastq'
        We will combine all fastq files at the end.
    """
    # Check if this job hasn't already been completed
    if subfolder.fastq_moved:
        return  # Job has already been

    # Dict where index is barcode or just "unbarcoded" and value is a list of paths
    # Generally just one item in the list of paths
    fastq_files_dict = {}

    if BARCODING:
        barcodes = [barcode for barcode in os.listdir(subfolder.workspace_dir)
                    if os.path.isdir(os.path.join(subfolder.workspace_dir, barcode)) and barcode.startswith("barcode")
                    or os.path.isdir(os.path.join(subfolder.workspace_dir, barcode)) and barcode == "unclassified"]
        for barcode in barcodes:
            # Get fastq files in the workspace directory
            fastq_files_list = [os.path.join(subfolder.workspace_dir, "pass", barcode, fastq)
                                for fastq in os.listdir(os.path.join(subfolder.workspace_dir, "pass", barcode))
                                if fastq.endswith(".fastq")]
            fastq_files_dict[barcode] = fastq_files_list
    else:
        fastq_files_list = [os.path.join(subfolder.workspace_dir, "pass", fastq)
                            for fastq in os.listdir(os.path.join(subfolder.workspace_dir, "pass"))
                            if fastq.endswith(".fastq")]
        if CHOSEN_FLOWCELL == "SQK-LSK308":
            fastq_files_dict["1d"] = fastq_files_list
        else:
            fastq_files_dict["unbarcoded"] = fastq_files_list

    if CHOSEN_FLOWCELL == "SQK-LSK308":
        # Get the 1dsq fastq file as well
        fastq_folder_1dsq = os.path.join(subfolder.albacore_dir, "1dsq_analysis", "1dsq_analysis", "workspace")
        fastq_files_list = [os.path.join(fastq_folder_1dsq, "pass", fastq)
                            for fastq in os.listdir(os.path.join(fastq_folder_1dsq, "pass"))
                            if fastq.endswith(".fastq")]
        fastq_files_dict["1dsq"] = fastq_files_list

    for barcode, fastq_files in fastq_files_dict.items():
        if len(fastq_files) > 1:
            LOGGER.info("It seems like we have more than 1 file here", fastq_files)

        fastq_file_index = -1

        for fastq_file in fastq_files:
            fastq_file_index += 1
            if BARCODING or CHOSEN_FLOWCELL == "SQK-LSK308":
                new_fastq_file = subfolder.fastq_file.replace(".fastq", '.'.join(["",
                                                                                  barcode,
                                                                                  str(fastq_file_index),
                                                                                  "fastq"]))
            else:
                new_fastq_file = subfolder.fastq_file.replace(".fastq", '.'.join(["",
                                                                                  str(fastq_file_index),
                                                                                  "fastq"]))

            # Create move command and run through subprocess.
            move_command = f"mv {fastq_file} {os.path.join(FASTQ_DIR, new_fastq_file)}"
            LOGGER.info(f"Moving fastq files in {subfolder.workspace_dir} to {FASTQ_DIR}")
            move_proc = subprocess.Popen(move_command, shell=True,
                                         stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = move_proc.communicate()
            stdout = stdout.decode()
            stderr = stderr.decode()
            if not stdout == "" or not stderr == "":
                LOGGER.info(f"Output of mv command is:\nStdout:\"{stdout.rstrip()}\"\nStderr:\"{stderr.rstrip()}\"")

    # Set as complete so we don't have to do it again.
    subfolder.fastq_moved = True
    update_dataframe(subfolder)


def tar_albacore_folder(subfolder):
    """
    Tar up the albacore folder into .tar.gz files.
    It's really just a bunch of metadata anyway
    """
    os.chdir(ALBACORE_DIR)
    tar_command = f"tar -cf - {os.path.basename(os.path.normpath(subfolder.albacore_dir))} --remove-files | " + \
                  f"pigz -p 16 > {subfolder.albacore_zip_file}"

    # Run tar_command through subprocess.
    tar_proc = subprocess.Popen(tar_command, shell=True,
                                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    if not stdout == "" or not stderr == "":
        LOGGER.info(f"Output of tar albacore folder command is:\n\"{stdout.rstrip()}\"\n\"{stderr.rstrip()}\"")
    subfolder.albacore_tarred = True
    os.chdir(READS_DIR)


def merge_fastq_files_wrapper():
    """ Merge all of the fastq files in the fastq directory to all.fastq"""

    if BARCODING:
        # fastq_file: 0003_49335_plasmids.barcode07.0.fastq
        # Get list of final fastq files
        fastq_files = [os.path.join(FASTQ_DIR, fastq_file) for fastq_file in os.listdir(FASTQ_DIR)
                       if len(fastq_file.split(".")) >= 3  
                       and (fastq_file.split(".")[-3].startswith("barcode")
                            or fastq_file.split(".")[-3] == "unclassified")
                       and "all" not in fastq_file]
        # Group fastq files by barcode
        fastq_by_barcodes = {}
        for fastq_file in fastq_files:
            barcode = fastq_file.split(".")[-3]
            if barcode not in fastq_by_barcodes:
                fastq_by_barcodes[barcode] = [fastq_file]
            else:
                fastq_by_barcodes[barcode].append(fastq_file)
        for barcode, fastq_files in fastq_by_barcodes.items():
            concatenated_fastq_file = os.path.join(FASTQ_DIR, barcode + ".all.fastq")
            merge_fastq_files(fastq_files, concatenated_fastq_file)
    else:
        concatenated_fastq_file = os.path.join(FASTQ_DIR, "all.fastq")
        fastq_files = [os.path.join(FASTQ_DIR, fastq) for fastq in os.listdir(FASTQ_DIR)
                       if fastq.endswith(".fastq") and not fastq == "all.fastq"]
        merge_fastq_files(fastq_files, concatenated_fastq_file)


def merge_fastq_files(fastq_files, concatenated_fastq_file):
    LOGGER.info(f"Concatenating all of the fastq files in {FASTQ_DIR} to {concatenated_fastq_file}")
    for fastq_file in fastq_files:
        concat_command = f"cat {fastq_file} >> {concatenated_fastq_file}"
        concat_proc = subprocess.Popen(concat_command, shell=True,
                                       stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = concat_proc.communicate()
        if not stdout == "" or not stderr == "":
            LOGGER.info(f"Output of cat command is:\nStdout:\"{stdout.rstrip()}\"\nStderr:\"{stderr.rstrip()}\"")

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


def get_subfolders():
    """ New subfolders will emerge to the directory as tarballs.
        Get a list of new tarballs then append them to the SUBFOLDERS
        list
    """
    global SUBFOLDERS
    # Get a list of tarballs
    subfolders = [tarred_folder.replace(".tar.gz", "") for tarred_folder
                  in os.listdir(READS_DIR) if tarred_folder.endswith(".tar.gz")]
    for subfolder in subfolders:
        is_initialised = False
        for initialised_subfolder in SUBFOLDERS:
            if subfolder == initialised_subfolder.name:
                is_initialised = True
                break
        if not is_initialised:
            SUBFOLDERS.append(Subfolder(subfolder))


def is_still_transferring():
    """ Is data still coming from the laptop.
        Get list of files in the destination directory
        and if any of them are called TRANSFERRING then the answer is yes!
    """
    transferring_file = [t_file for t_file in os.listdir(PARENT_DIRECTORY)
                         # Lock name is called TRANSFERRING.
                         if t_file == "TRANSFERRING"]
    # Can the lock file file be found?
    if len(transferring_file) == 0:
        return False
    else:
        return True


def take_a_break():
    """
    Just a 15 second sleep ;)
    """
    time.sleep(15)


def generate_dataframe():
    global STATUS_DF

    for subfolder in SUBFOLDERS:
        if subfolder.name not in STATUS_DF.index.tolist():
            # Append row
            STATUS_DF = STATUS_DF.append(subfolder.to_series(), ignore_index=True)
            # Reset index
            STATUS_DF.set_index('name', drop=False, inplace=True)

    STATUS_DF.to_csv(STATUS_CSV, index=False)


def update_dataframe(subfolder):
    global STATUS_DF
    """
    Edits the STATUS_DF row of a given subfolder
    """
    STATUS_DF.loc[STATUS_DF.index == subfolder.name, ] = subfolder.to_series().tolist()


def is_still_basecalling():
    """Any subfolders left to be basecalled?"""
    if len([subfolder for subfolder in SUBFOLDERS
            if not subfolder.albacore_complete]) == 0:
        return False  # basecalling finished
    else:  # Basecalling still going
        return True


def get_current_subfolder(subfolder_name):
    return [subfolder for subfolder in SUBFOLDERS
            if subfolder.name == subfolder_name][0]


def num_current_jobs():
    """Returns the number of subfolders that have either an extraction job submitted and not complete
    or an albacore job submitted and not complete"""
    return len([subfolder for subfolder in SUBFOLDERS
                if subfolder.extracted_submitted
                and not subfolder.extracted_complete
                or subfolder.albacore_submitted
                and not subfolder.albacore_complete])

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
    if QSUB_TYPE == "SGE":
        qacct_command = f"qacct -j {job_id} | grep start_time | grep -v '\-\/\-'"
    elif QSUB_TYPE == "TORQUE":
        qacct_command = f"tracejob {job_id} | grep QUEUED | wc -l"
    elif QSUB_TYPE == "SLURM":
        qacct_command = f"sacct -j {job_id} | tail -n 1 | grep RUNNING | wc -l"
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    qacct_stdout = qacct_stdout.decode()
    qacct_stderr = qacct_stderr.decode()
    if QSUB_TYPE == "SGE" and qacct_stdout == "" \
            or QSUB_TYPE == "TORQUE" and int(qacct_stdout) > 1\
            or QSUB_TYPE == "SLURM" and int(qacct_stdout) > 1:
        # Start-time is non-existent
        return False
    if not qacct_stderr == "":
        LOGGER.info(f"qacct stderr of {qacct_stderr}") 
    return True


def has_completed(job_id, check_failed):
    """
    Use the qacct command to see if the job has finished.
    No end_time is represented as -/-
    """
    if QSUB_TYPE == "SGE":
        qacct_command = f"qacct -j {job_id} | grep end_time | grep -v '\-\/\-'"
    elif QSUB_TYPE == "TORQUE":
        qacct_command = f"tracejob {job_id} 2> /dev/null | grep Exit_status "
    elif QSUB_TYPE == "SLURM":
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
        LOGGER.info(f"qacct stderr of {qacct_stderr}")
    
    if check_failed and has_failed(job_id):
        LOGGER.info(f"It appears that the job {job_id} has failed")
        sys.exit(f"Failing because {job_id} failed.")

    return True


def has_failed(job_id):
    """
    Use the qacct command to see if the job has failed.
    We pass the exit_status parameter into awk and sum it.
    If it's any greater than zero then the command has failed
    """
    if QSUB_TYPE == "SGE":
        qacct_command = f"qacct -j {job_id} | grep exit_status | awk '{{sum+=$2}} END {{print sum}}'"
    elif QSUB_TYPE == "TORQUE":
        qacct_command = f"tracejob {job_id} 2> /dev/null | grep Exit_status | grep - v preparing |" + \
                        f" cut - d' ' - f7 | cut - d'=' -f2"
    elif QSUB_TYPE == "SLURM":
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
        LOGGER.info(qacct_stdout)
        # Job has failed
        return True
    if not qacct_stderr == "":
        LOGGER.info(f"qacct stderr of {qacct_stderr}")
    return False


def create_basecalling_lock_file():
    # Create a file called BASECALLING in the working directory
    Path(BASECALLING_LOCK_FILE).touch()


def remove_basecalling_lock_file():
    # Remove the basecalling lock file.
    os.remove(BASECALLING_LOCK_FILE)