#!/usr/bin/env python

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
from datetime import datetime

# Set global and semi-global variables
TRANSFERRING = True
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
QSUB_HOST = ""
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
LOGGER = None
LOGGER_DIR = ""
LOGGER_PATH = ""
BARCODING = False
QSUB_SOURCE_STRING = ""


class Subfolder:
    def __init__(self, name):
        self.name = name
        # Get the relative name of the tar file
        self.tar_filename = name + ".tar.gz"
        # Set personal directories up
        self.reads_dir = os.path.join(READS_DIR, name)
        self.albacore_dir = os.path.join(ALBACORE_DIR, name)
        # Albacore related files
        self.workspace_dir = os.path.join(self.albacore_dir, "workspace")
        self.albacore_summary_file = os.path.join(self.albacore_dir, "sequencing_summary.txt")
        self.albacore_log_file = os.path.join(self.albacore_dir, "pipeline.log")
        self.fastq_file = name + ".fastq"
        self.albacore_zip_file = self.albacore_dir + ".albacore.tar.gz"
        # Qsub related files / ids
        self.extracted_qsub_output_log = os.path.join(QSUB_LOG_DIR, name + ".extract.o.log")
        self.extracted_qsub_error_log = os.path.join(QSUB_LOG_DIR, name + ".extract.e.log")
        self.extracted_jobid = -1
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


def main():
    global TRANSFERRING
    # Basic house cleaning, get arguments, make sure they're legit.
    args = get_arguments()
    set_global_variables(args) 
    check_directories()
    set_logger()
    
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

        # Have a break if no new subfolders
        if not new_subfolders():
            take_a_break()
            continue

        # Otherwise run through the pipeline
        run_pipeline()

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
    merge_fastq_files_wrapper()


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
    parser.add_argument("--qsub_source_string", type=str, required=False, default=None,
                        help="Anything that is required to be sourced prior to read_fast5_basecaller on qsub host?")

    return parser.parse_args()


def set_global_variables(args):
    # Global variables
    global READS_DIR, ALBACORE_DIR, WORKING_DIR, NUM_THREADS, CHOSEN_KIT, FASTQ_DIR
    global STATUS_CSV, QSUB_LOG_DIR, CW_DIR, QSUB_HOST, CHOSEN_FLOWCELL, MAX_PROCESSES
    global LOGGER_DIR, BARCODING, QSUB_SOURCE_STRING
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
    if args.qsub_source_string is not None:
        QSUB_SOURCE_STRING = args.qsub_source_string


def check_directories():
    # Make sure the directories exist, change to reads directory,
    # Create any other necessary directories for the script to run.
    global READS_DIR, ALBACORE_DIR, PARENT_DIRECTORY, QSUB_LOG_DIR, FASTQ_DIR
    global STATUS_CSV, LOGGER_DIR, LOGGER_PATH
    if not os.path.isdir(READS_DIR):
        sys.exit("Error, %s does not exist" % READS_DIR)

    READS_DIR = os.path.abspath(READS_DIR) + "/"
    os.chdir(READS_DIR)

    PARENT_DIRECTORY = os.path.abspath(os.path.join(READS_DIR, os.pardir)) + "/"

    if ALBACORE_DIR == "":
        ALBACORE_DIR = PARENT_DIRECTORY + "albacore/"
    else:
        ALBACORE_DIR = os.path.abspath(ALBACORE_DIR)

    if not os.path.isdir(ALBACORE_DIR):
        os.mkdir(ALBACORE_DIR)

    QSUB_LOG_DIR = PARENT_DIRECTORY + "qsub_log/"
    if not os.path.isdir(QSUB_LOG_DIR):
        os.mkdir(QSUB_LOG_DIR)

    if FASTQ_DIR == "":
        FASTQ_DIR = PARENT_DIRECTORY + "fastq/"
    else:
        FASTQ_DIR = os.path.abspath(FASTQ_DIR)

    if not os.path.isdir(FASTQ_DIR):
        os.mkdir(FASTQ_DIR)
    if LOGGER_DIR == "":
        LOGGER_DIR = PARENT_DIRECTORY + "poreduck_logs/"
    else:
        LOGGER_DIR = os.path.abspath(LOGGER_DIR)
    if not os.path.isdir(LOGGER_DIR):
        os.mkdir(LOGGER_DIR)
    LOGGER_PATH = os.path.join(LOGGER_DIR, '.'.join([str(datetime.now().date()), str(datetime.now().time()), "albacore_pipeline.log"]))

    if not STATUS_CSV == "":
        if not os.path.isfile(os.path.join(CW_DIR, STATUS_CSV)):
            sys.exit("Resume csv specified but does not exist")
    else:
        STATUS_CSV = os.path.join(PARENT_DIRECTORY, "status.csv")
        if os.path.isfile(STATUS_CSV):
            print("Warning, STATUS_CSV not defined but still exists. Removing file %s" % STATUS_CSV)
            os.remove(STATUS_CSV)


def set_logger():
    global LOGGER
    logging.basicConfig(filename=LOGGER_PATH, level=logging.INFO, format='%(asctime)s::\t%(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
    LOGGER = logging.getLogger() 


def pick_up_from_previous_run():
    global SUBFOLDERS
    if os.path.isfile(STATUS_CSV):
        previous_dataframe = pd.read_csv(STATUS_CSV)
        LOGGER.info('Found %s. Putting run back on track.' % STATUS_CSV)

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
                LOGGER.info("It appears %s has failed albacore processing. " \
                            "Reverting albacore submitted, commenced to false. " \
                            "Removing job id. %s" % (SUBFOLDERS[-1].name, SUBFOLDERS[-1].albacore_jobid))
                SUBFOLDERS[-1].albacore_submitted = False
                SUBFOLDERS[-1].albacore_commenced = False
                SUBFOLDERS[-1].albacore_jobid = -1
                SUBFOLDERS[-1].albacore_complete = False
            if has_failed(SUBFOLDERS[-1].extracted_jobid):
                LOGGER.info("It appears %s has failed to extract. " \
                            "Reverting albacore submitted, commenced to false. " \
                            "Removing job id. %s" % (SUBFOLDERS[-1].name, SUBFOLDERS[-1].extracted_jobid))
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
    elif num_current_jobs() > MAX_PROCESSES and not MAX_PROCESSES == 0:
        return
    else:
        subfolder.extracted_submitted = True

    tar_command = "pigz -dc %s | tar -xf -" % os.path.join(READS_DIR, subfolder.tar_filename)
    qsub_command = "qsub -o %s -e %s -N PIGZ -S /bin/bash -wd %s -l hostname=%s" % \
                   (subfolder.extracted_qsub_output_log,
                    subfolder.extracted_qsub_error_log,
                    READS_DIR, QSUB_HOST)
    tar2qsub_command = "echo \"%s\" | %s" % (tar_command, qsub_command)
    LOGGER.info("Submitting the following extraction job to sge:\n \"%s\"" % tar2qsub_command)
    tar_proc = subprocess.Popen(tar2qsub_command, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()

    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of the array.
    LOGGER.info("Output of pigz sge submission \nStdout:\"%s\"\nStderr:\"%s\"" % (stdout.rstrip(), stderr.rstrip()))
    subfolder.extracted_jobid = int(stdout.rstrip().split()[2])
    update_dataframe(subfolder)


def run_albacore(subfolder):
    """ Generate read_fast5_basecaller command
        And qsub command and pipe former into latter
        This function also sets the jobid of a given folder.
    """
    if subfolder.albacore_submitted:
        return   # Basecalling has already queued for this directory
    elif num_current_jobs() > MAX_PROCESSES and not MAX_PROCESSES == 0:
        return
    else:
        # Commence basecalling on this directory
        subfolder.albacore_submitted = True

    # Number of gigabytes required for a given qsub command
    memory_allocation = 4 + 2*NUM_THREADS
    # The read_fast5_basecaller is the algorithm that does the actual base calling,
    # what would be run if we just had Ubuntu.
    basecaller_command = "read_fast5_basecaller.py " \
                         "--input %s " \
                         "--worker_threads %s " \
                         "--save_path %s " \
                         "--flowcell %s " \
                         "--kit %s" \
                         % (subfolder.reads_dir, NUM_THREADS, subfolder.albacore_dir, CHOSEN_FLOWCELL, CHOSEN_KIT)
    if BARCODING:
        basecaller_command += " --barcoding"

    if not QSUB_SOURCE_STRING == "":
        basecaller_command = QSUB_SOURCE_STRING + ' && ' + basecaller_command

    # These are both parsed into qsub which then determines what to do with it all.
    qsub_command = "qsub " \
                   "-o %s " \
                   "-e %s " \
                   "-S /bin/bash -l hostname=%s " \
                   "-l h_vmem=%dG " \
                   "-N ALBACORE " \
                   "-wd %s -v OMP_NUM_THREADS=1" \
                   % (subfolder.albacore_qsub_output_log, subfolder.albacore_qsub_error_log,
                      QSUB_HOST, memory_allocation, PARENT_DIRECTORY)

    # Put these all together into one grand command
    albacore_command = "echo \"%s\" | %s" % (basecaller_command, qsub_command) 
    LOGGER.info("Submitting the following job to SGE:\n\"%s\"" % albacore_command)    

    # Execute qsub command via subprocess.
    albacore_proc = subprocess.Popen(albacore_command, shell=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = albacore_proc.communicate()
    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of the array.
    LOGGER.info("Output of albacore sge submission \nStdout:\"%s\"\nStderr:\"%s\"" % (stdout.rstrip(), stderr.rstrip()))
    subfolder.albacore_jobid = int(stdout.rstrip().split()[2])
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
    if not os.path.join(READS_DIR, subfolder.tar_filename + ".tar.gz"):
        LOGGER.info("Um... the archive doesn't exist, I'm not going to delete the folder")
        return

    # Generate remove command
    remove_command = "rm -rf %s" % os.path.join(READS_DIR, subfolder.name)

    # Run remove command through subprocess
    remove_proc = subprocess.Popen(remove_command, shell=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = remove_proc.communicate()
    if not stdout=="" or stderr=="":
        LOGGER.info("Output of deleting folder command is\nStdout:\"%s\"\nStderr:\"%s\"" % (stdout.rstrip(), stderr.rstrip()))
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

    fastq_files_dict = {}  # Dict where index is barcode or just "unbarcoded" and value is a list of paths
                           # Generally just one item in the list of paths
    if BARCODING:
        barcodes = [barcode for barcode in os.listdir(subfolder.workspace_dir)
                    if os.path.isdir(os.path.join(subfolder.workspace_dir, barcode)) and barcode.startswith("barcode")
                    or os.path.isdir(os.path.join(subfolder.workspace_dir, barcode)) and barcode == "unclassified"]
        for barcode in barcodes:
            # Get fastq files in the workspace directory
            fastq_files_list = [os.path.join(subfolder.workspace_dir, barcode, fastq) for fastq in os.listdir(os.path.join(subfolder.workspace_dir, barcode))
                                if fastq.endswith(".fastq")]
            fastq_files_dict[barcode] = fastq_files_list
    else:
        fastq_files_list = [os.path.join(subfolder.workspace_dir, fastq) for fastq in os.listdir(subfolder.workspace_dir)
                            if fastq.endswith(".fastq")]
        fastq_files_dict["unbarcoded"] = fastq_files_list

    for barcode, fastq_files in fastq_files_dict.iteritems():
        if len(fastq_files) > 1:
            LOGGER.info("It seems like we have more than 1 file here", fastq_files)

        fastq_file_index = -1

        for fastq_file in fastq_files:
            fastq_file_index += 1
            if BARCODING:
                new_fastq_file = subfolder.fastq_file.replace(".fastq", '.'.join(["",
                                                                                  barcode,
                                                                                  str(fastq_file_index),
                                                                                  "fastq"]))
            else:
                new_fastq_file = subfolder.fastq_file.replace(".fastq", '.'.join(["",
                                                                                  str(fastq_file_index),
                                                                                  "fastq"]))

            # Create move command and run through subprocess.
            move_command = "mv %s %s" % (fastq_file,
                                         os.path.join(FASTQ_DIR, new_fastq_file))
            LOGGER.info("Moving fastq files in %s to %s" % (subfolder.workspace_dir, FASTQ_DIR))
            move_proc = subprocess.Popen(move_command, shell=True,
                                         stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = move_proc.communicate()
            if not stdout == "" or not stderr == "":
                LOGGER.info("Output of mv command is:\nStdout:\"%s\"\nStderr:\"%s\"" % (stdout.rstrip(), stderr.rstrip()))

    # Set as complete so we don't have to do it again.
    subfolder.fastq_moved = True
    update_dataframe(subfolder)


def tar_albacore_folder(subfolder):
    """
    Tar up the albacore folder into .tar.gz files.
    It's really just a bunch of metadata anyway
    """
    os.chdir(ALBACORE_DIR)
    tar_command = "tar -cf - %s --remove-files | pigz -p 16 > %s" % \
                  (os.path.basename(os.path.normpath(subfolder.albacore_dir)), subfolder.albacore_zip_file)

    # Run tar_command through subprocess.
    tar_proc = subprocess.Popen(tar_command, shell=True,
                                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()
    if not stdout == "" or not stderr == "":
        LOGGER.info("Output of tar albacore folder command is:\n\"%s\"\n\"%s\"" % (stdout.rstrip(), stderr.rstrip()))
    subfolder.albacore_tarred = True
    os.chdir(READS_DIR)


def merge_fastq_files_wrapper():
    """ Merge all of the fastq files in the fastq directory to all.fastq"""

    if BARCODING:
        # fastq_file: 0003_49335_plasmids.barcode07.0.fastq
        # Get list of final fastq files
        fastq_files = [os.path.join(FASTQ_DIR, fastq_file) for fastq_file in os.listdir(FASTQ_DIR)
                       if (fastq_file.split(".")[-3].startswith("barcode")
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
        for barcode, fastq_files in fastq_by_barcodes.iteritems():
            concatenated_fastq_file = os.path.join(FASTQ_DIR, barcode + ".all.fastq")
            merge_fastq_files(fastq_files, concatenated_fastq_file)
    else:
        concatenated_fastq_file = os.path.join(FASTQ_DIR, "all.fastq")
        fastq_files = [os.path.join(FASTQ_DIR, fastq) for fastq in os.listdir(FASTQ_DIR)
                       if fastq.endswith(".fastq") and not fastq == "all.fastq"]
        merge_fastq_files(fastq_files, concatenated_fastq_file)


def merge_fastq_files(fastq_files, concatenated_fastq_file):
    LOGGER.info("Concatenating all of the fastq files in %s to %s" % (FASTQ_DIR, concatenated_fastq_file))
    for fastq_file in fastq_files:
        concat_command = "cat %s >> %s" % (fastq_file, concatenated_fastq_file)
        concat_proc = subprocess.Popen(concat_command, shell=True,
                                       stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = concat_proc.communicate()
        if not stdout == "" or not stderr == "":
            LOGGER.info("Output of cat command is:\nStdout:\"%s\"\nStderr:\"%s\"" % (stdout.rstrip(), stderr.rstrip()))

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


def new_subfolders():
    """Any new subfolders that haven't been added to the basecalling queue"""
    for subfolder in SUBFOLDERS:
        if not subfolder.albacore_submitted:
            return True
    return False


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
    qacct_command = "qacct -j %s | grep start_time | grep -v '\-\/\-'" % job_id
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    if qacct_stdout == "":
        # Start-time is non-existent
        return False
    if not qacct_stderr == "":
        LOGGER.info("qacct stderr of %s", qacct_stderr)
    return True


def has_completed(job_id, check_failed):
    """
    Use the qacct command to see if the job has finished.
    No end_time is represented as -/-
    """
    qacct_command = "qacct -j %s | grep end_time | grep -v '\-\/\-'" % job_id
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    if qacct_stdout == "":
        # end_time is non-existent
        return False
    if not qacct_stderr == "":
        LOGGER.info("qacct stderr of %s", qacct_stderr)
    
    if check_failed and has_failed(job_id):
        LOGGER.info("It appears that the job %d has failed" % job_id)
        sys.exit("Failing because %d failed. Good one Dave" % job_id)

    return True


def has_failed(job_id):
    """
    Use the qacct command to see if the job has failed.
    We pass the exit_status parameter into awk and sum it.
    If it's any greater than zero then the command has failed
    """
    qacct_command = "qacct -j {0} | grep exit_status | awk '{{sum+=$2}} END {{print sum}}'".format(job_id)
    
    qacct_proc = subprocess.Popen(qacct_command, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    qacct_stdout, qacct_stderr = qacct_proc.communicate()
    if qacct_stdout.strip() == "":
        return True  
    if int(qacct_stdout) > 0:
        LOGGER.info(qacct_stdout)
        # Job has failed
        return True
    if not qacct_stderr == "":
        LOGGER.info("qacct stderr of %s", qacct_stderr)
    return False


# Run the main function
main()
