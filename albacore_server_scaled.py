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
ALBACORE_DIR = ""
WORKING_DIR = ""
NUM_THREADS = 0
CHOSEN_CONFIG = ""
PARENT_DIRECTORY = ""
QSUB_LOG_DIR = ""
FASTQ_DIR = ""
ALBACORE_QSUBJOB_BY_FOLDER = {}
FASTQ_QSUBJOB_BY_FOLDER = {}

CONFIGS = {"FC106_RAD001": "r94_250bps_linear.cfg",  # Rapid sequencing
           "FC106_LSK208_2d": "r94_250bps_2d.cfg",   # 2D unsure which one
           "FC106_LSK108": "r94_450bps_linear.cfg",  # For 1D ligation sequencing.
           "FC106_RAD002": "r94_450bps_linear.cfg"}  # Second Rapid sequencing kit.


def main():
    global TRANSFERRING
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

        for folder, job in ALBACORE_QSUBJOB_BY_FOLDER.iteritems():
            if is_job_complete(job):
                perform_fastq_extraction(folder)

        if not is_still_transferring():
            TRANSFERRING = False

        if len(tarred_read_sets) == 0:
            take_a_break()

    # While albacore still running
    albacore_processing = True
    while albacore_processing:
        albacore_processing = False
        for folder, job in ALBACORE_QSUBJOB_BY_FOLDER.iteritems():
            if not is_job_complete(job):
                albacore_processing = True
                continue
            perform_fastq_extraction(folder)

    # While fastq still running
    fastq_processing = True
    while fastq_processing:
        fastq_processing = False
        for folder, job in FASTQ_QSUBJOB_BY_FOLDER.iteritems():
            if not is_job_complete(job):
                fastq_processing = True
                continue


def get_arguments():
    parser = argparse.ArgumentParser(
        description="The albacore-server scaled command incorporates qsub to spread the server load and rapidly " +
                    "generate data off the MinION.")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="/path/to/reads, should have a bunch of tar zipped files in it.")
    parser.add_argument("--config", type=str, choices=CONFIGS.keys(),
                        help="Pick a config")
    parser.add_argument("--output_dir", type=str, required=False, default=None,
                        help="Will be called 'albacore' and sit adjacent to the reads folder if left blank.")
    parser.add_argument("--num_threads", type=int, required=False, default=5,
                        help="How many threads did you wish to use per parallel output, default is 5.")
    parser.add_argument("--fastq_dir", type=str, required=False, default=None,
                        help="Where should the fastq data be placed. Will be called 'fastq' " +
                             "and sit adjacent to reads folder if left blank.")
    return parser.parse_args()


def set_global_variables(args):
    # Global variables
    global READS_DIR, ALBACORE_DIR, WORKING_DIR, NUM_THREADS, CHOSEN_CONFIG, FASTQ_DIR
    READS_DIR = args.reads_dir
    if args.output_dir is not None:
        ALBACORE_DIR = args.output_dir
    CHOSEN_CONFIG = CONFIGS[args.config]
    NUM_THREADS = args.num_threads
    if args.fastq_dir is not None:
        FASTQ_DIR = args.fastq_dir


def check_directories():
    # Make sure the directories exist, change to reads directory,
    # Create any other necessary directories for the script to run.
    global READS_DIR, ALBACORE_DIR, PARENT_DIRECTORY, QSUB_LOG_DIR, FASTQ_DIR
    if not os.path.isdir(READS_DIR):
        sys.exit("Error, %s does not exist" % READS_DIR)
    READS_DIR = os.path.abspath(READS_DIR) + "/"
    os.chdir(READS_DIR)

    PARENT_DIRECTORY = os.path.abspath(os.path.join(READS_DIR, os.pardir)) + "/"

    if ALBACORE_DIR == "":
        ALBACORE_DIR = PARENT_DIRECTORY + "albacore/"
    if not os.path.isdir(ALBACORE_DIR):
        os.mkdir(ALBACORE_DIR)

    QSUB_LOG_DIR = PARENT_DIRECTORY + "qsub_log/"

    if not os.path.isdir(QSUB_LOG_DIR):
        os.mkdir(QSUB_LOG_DIR)
   
    if FASTQ_DIR == "":
        FASTQ_DIR = PARENT_DIRECTORY + "fastq/"
    
    if not os.path.isdir(FASTQ_DIR):
        os.mkdir(FASTQ_DIR)


def get_tarred_files():
    # Get the tarred files, note, must not be already be being base called, vicious cycle!
    tarred_files = [READS_DIR + tarred_file for tarred_file in os.listdir(READS_DIR)
                    # If the file is a zip file and
                    if tarred_file.endswith(".tar.gz")
                    # albacore folder does not already exist
                    and not os.path.isdir(ALBACORE_DIR + tarred_file.replace(".tar.gz", ""))]
    return tarred_files


def extract_tarred_read_set(tar_file):
    # Extract the tarred reads, we currently do not delete after extraction, maybe a good idea?
    tar_command = "tar -xf %s" % tar_file
    tar_proc = subprocess.Popen(tar_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()
    print("Output of tar command", stdout, stderr)


def run_albacore(tarred_read_set):
    # Combine albacore and qsub commands together to be run on one line.
    folder = tarred_read_set.replace(".tar.gz", "/")
    qsub_log_file = QSUB_LOG_DIR + folder.split("/")[-2] + "albacore.o.log"
    qsub_error_file = QSUB_LOG_DIR + folder.split("/")[-2] + "albacore.e.log"
    output_folder = ALBACORE_DIR + folder.split("/")[-2]
    memory_allocation = 4 + NUM_THREADS  # Number of gigabytes required for a given qsub command
    # The read_fast5_basecaller is the algorithm that does the actual base calling,
    # what would be run if we just had Ubuntu.
    basecaller_command = "read_fast5_basecaller.py " \
                         "--input %s " \
                         "--worker_threads %s " \
                         "--save_path %s " \
                         "--config %s" \
                         % (folder, NUM_THREADS, output_folder, CHOSEN_CONFIG)

    # These are both parsed into qsub which then determines what to do with it all.
    qsub_command = "qsub -o %s -e %s -S /bin/bash -l h_vmem=%dG" % (qsub_log_file, qsub_error_file, memory_allocation)

    # Put these all together into one grand command
    albacore_command = "echo \"%s\" | %s " % (basecaller_command, qsub_command)
    print(albacore_command)

    # Execute qsub command via subprocess.
    albacore_proc = subprocess.Popen(albacore_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = albacore_proc.communicate()
    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of the array.
    print("Output of albacore command", stdout, stderr) 
    ALBACORE_QSUBJOB_BY_FOLDER[folder.split("/")[-2]] = stdout.rstrip().split()[2]
    print("Output of albacore_proc", stdout, stderr)


def get_albacore_subfolders():
    albacore_subfolders = [ALBACORE_DIR + subfolder + "/" for subfolder in os.listdir(ALBACORE_DIR)
                           if os.path.isdir(ALBACORE_DIR + subfolder)]
    return albacore_subfolders


def perform_fastq_extraction(albacore_folder):
    # If the albacore job has completed.
    # Use subprocess to run the fastq extraction on the folder.

    fastq_file = FASTQ_DIR + albacore_folder + ".fastq"
    qsub_log_file = QSUB_LOG_DIR + albacore_folder + ".fastq.o.log"
    qsub_error_file = QSUB_LOG_DIR + albacore_folder + ".fastq.e.log"

    poretools_command = "poretools fastq %s > %s" % (ALBACORE_DIR + albacore_folder, fastq_file)
    qsub_command = "echo \"%s\" | qsub -o %s -e %s -S /bin/bash" % \
                   (poretools_command, qsub_log_file, qsub_error_file)

    # Fastq command must not already be starting processing
    if albacore_folder not in FASTQ_QSUBJOB_BY_FOLDER:
        qsub_proc = subprocess.Popen(qsub_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = qsub_proc.communicate()
        # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
        # So job equal to third element of the array.
        FASTQ_QSUBJOB_BY_FOLDER[albacore_folder] = stdout.rstrip().split()[2]
        print("Output of fastq command for ", albacore_folder, stdout, stderr)


def is_still_transferring():
    # Get list of files in dest directory.
    transferring_file = [t_file for t_file in os.listdir(PARENT_DIRECTORY)
                         if t_file == "TRANSFERRING"]  # Lock name is called transferring.
    # Can the file be found.
    if len(transferring_file) == 0:
        return False
    else:
        return True


def is_job_complete(job):
    # Use qstat -u '*' to get all the jobs. We'll pass this into a pandas dataframe
    get_jobs_command = "qstat -u '*'"
    get_jobs_proc = subprocess.Popen(get_jobs_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    get_jobs_output, get_jobs_stderr = get_jobs_proc.communicate()

    # job id is the first column of each row
    #  120882 0.55500 supernova  user...
    for line in get_jobs_output.split("\n"):
        # Now split by space.
        first_column = line.split()[0]
        if first_column == job:
            return False  # Job still in queue

    # Otherwise we haven't found the job so it must be complete
    return True


def take_a_break():
    time.sleep(60)


main()
