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
SUBFOLDERS = []
STATUS_DF = pd.DataFrame(columns=['name', 'extracted', 'albacore_commenced',
                                  'albacore_complete', 'folder_removed', 'fastq_moved'])

CONFIGS = {"FC106_RAD001": "r94_250bps_linear.cfg",  # Rapid sequencing
           "FC106_LSK208_2d": "r94_250bps_2d.cfg",   # 2D unsure which one
           "FC106_LSK108": "r94_450bps_linear.cfg",  # For 1D ligation sequencing.
           "FC106_RAD002": "r94_450bps_linear.cfg"}  # Second Rapid sequencing kit.


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
        # Qsub related files / ids
        self.albacore_qsub_output_log = os.path.join(QSUB_LOG_DIR, name + ".albacore.o.log")
        self.albacore_qsub_error_log = os.path.join(QSUB_LOG_DIR, name + ".albacore.e.log")
        self.albacore_jobid = -1
        # Status related attributes
        self.extracted = False
        self.albacore_commenced = False
        self.albacore_complete = False
        self.folder_removed = False
        self.fastq_moved = False

    def to_series(self):
        return pd.Series(data=[self.name, self.extracted,
                               self.albacore_commenced, self.albacore_complete,
                               self.folder_removed, self.fastq_moved])


def main():
    global TRANSFERRING
    # Basic house cleaning, get arguments, make sure they're legit.
    args = get_arguments()
    set_global_variables(args)
    check_directories()

    # While the transferring lock script exists.
    while TRANSFERRING:
        get_subfolders()

        # Have a break if no new subfolders
        if not new_subfolders():
            take_a_break()
            continue

        # Otherwise, extract tarballs
        [extract_tarred_read_set(subfolder) for subfolder in SUBFOLDERS
         if not subfolder.extracted]

        # Now run albacore
        [run_albacore(subfolder) for subfolder in SUBFOLDERS
         if not subfolder.albacore_commenced]

        # Check for complete jobs
        [is_job_complete(subfolder) for subfolder in SUBFOLDERS
         if not subfolder.albacore_complete]

        # Once albacore is complete on a subfolder
        # Remove folder from existence, leaving just the tarball again.
        [remove_folder(subfolder) for subfolder in SUBFOLDERS
         if subfolder.albacore_complete and not subfolder.folder_removed]

        # Move the fastq file to the FASTQ_DIR
        [move_fastq_file(subfolder) for subfolder in SUBFOLDERS
         if subfolder.albacore_complete and not subfolder.fastq_moved]

        # Check that we're still transferring
        if not is_still_transferring():
            TRANSFERRING = False

    # Transferring from laptop complete just wait for albacore
    # to finish.
    albacore_processing = True
    while albacore_processing:
        albacore_processing = False
        have_break = True

        # Check for complete jobs
        [is_job_complete(subfolder) for subfolder in SUBFOLDERS
         if not subfolder.albacore_complete]

        # Still jobs in the works?
        # If so we need to go around the loop again
        if len([subfolder for subfolder in SUBFOLDERS
                if not subfolder.albacore_complete]) > 0:
            albacore_processing = True  # basecalling still going

        # Remove any folders of completed jobs.
        # Set have a break to be false if we find any subfolders to do this for
        if len([remove_folder(subfolder) for subfolder in SUBFOLDERS
                if subfolder.albacore_complete and not
                subfolder.folder_removed]) > 0:
            have_break = False

        # Move fastq files of any completed jobs.
        # Set have a break to be false if we find any subfolders to do this for
        if len([move_fastq_file(subfolder) for subfolder in SUBFOLDERS
                if subfolder.albacore_complete and not
                subfolder.fastq_moved]) > 0:
            have_break = False

        # If no subfolders recently completed basecalling
        # have a break otherwise check for some more.
        if have_break:
            take_a_break()
            continue

    # Merge fastq files at the end of the run.
    merge_fastq_files()


def generate_dataframe():
    global STATUS_DF
    STATUS_DF = pd.DataFrame(columns=['name', 'extracted',
                                      'albacore_commenced', 'albacore_complete',
                                      'folder_removed', 'fastq_moved'])
    for subfolder in SUBFOLDERS:
        STATUS_DF.append(subfolder.to_series())

    STATUS_DF.to_csv(file="status.csv", index=False)


def new_subfolders():
    """Any new subfolders that haven't commenced basecalling"""
    for subfolder in SUBFOLDERS:
        if not subfolder.albacore_commenced:
            return True
    return False


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
        print("Um... the archive doesn't exist, I'm not going to delete the folder")
        return

    # Generate remove command
    remove_command = "rm -rf %s" % os.path.join(READS_DIR, subfolder.name)

    # Run remove command through subprocess
    remove_proc = subprocess.Popen(remove_command, shell=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = remove_proc.communicate()
    print("Output of deleting folder is", stdout, stderr)
    generate_dataframe()


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

    # Get fastq files in the workspace directory
    fastq_files = [fastq for fastq in subfolder.workspace_dir
                   if fastq.endswith(".fastq")]
    if len(fastq_files) > 1:
        print("It seems like we have more than 1 file here", fastq_files)

    fastq_file_index = 0
    for fastq_file in fastq_files:
        # Rename fastq to have a .0.fastq at the start
        subfolder.fastq_file = subfolder.fastq_file.replace(".fastq", "") + str(fastq_file_index) + ".fastq"
        # We will move and rename but make sure we're not overwriting anything
        if os.path.isfile(subfolder.fastq_file):
            fastq_file_index += 1
            print("Hmmm, looks like we might accidentally overwrite something here, adding 1 to the index")
            subfolder.fastq_file = subfolder.fastq_file.replace(".fastq", "") + str(fastq_file_index) + ".fastq"
        # Create move command and run through subproces.
        move_command = "mv %s %s" % (os.path.join(subfolder.workspace_dir, fastq_file),
                                     os.path.join(FASTQ_DIR, subfolder.fastq_file))
        move_proc = subprocess.Popen(move_command, shell=True,
                                     stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = move_proc.communicate()
        if not stdout == "" or not stderr == "":
            print("Output of mv command is:", stdout, stderr)

    # Set as complete so we don't have to do it again.
    subfolder.fastq_moved = True
    generate_dataframe()


def merge_fastq_files():
    """ Merge all of the fastq files in the fastq directory to all.fastq"""
    concatenated_fastq_file = os.path.join(FASTQ_DIR, "all.fastq")
    fastq_files = [os.path.join(FASTQ_DIR, fastq) for fastq in os.listdir(FASTQ_DIR)
                   if fastq.endswith(".fastq") and not fastq == "all.fastq"]

    for fastq_file in fastq_files:
        concat_command = "cat %s >> %s" % (fastq_file, concatenated_fastq_file)
        concat_proc = subprocess.Popen(concat_command, shell=True,
                                       stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = concat_proc.communicate()
        if not stdout == "" or not stderr == "":
            print("Output of cat command is:", stdout, stderr)


def get_arguments():
    parser = argparse.ArgumentParser(
        description="The albacore-server scaled command incorporates qsub to spread " +
                    "the server load and rapidly " + "generate data off the MinION.")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="/path/to/reads, " +
                             "should have a bunch of tar zipped files in it.")
    parser.add_argument("--config", type=str, choices=CONFIGS.keys(),
                        help="Pick a config")
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


def extract_tarred_read_set(subfolder):
    # Extract the tarred reads, we currently do not delete after extraction...
    # maybe a good idea?

    # Has the dataset already been extracted?
    if subfolder.extracted:
        return
    else:
        subfolder.extracted = True  # It has been extracted now.

    tar_command = "tar -xf %s" % os.path.join(READS_DIR, subfolder.tar_filename)
    tar_proc = subprocess.Popen(tar_command, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()
    if not stdout == "" or not stderr == "":
        print("Output of tar command", stdout, stderr, subfolder.name)
    generate_dataframe()


def run_albacore(subfolder):
    """ Generate read_fast5_basecaller command
        And qsub command and pipe former into latter
        This function also sets the jobid of a given folder.
    """
    if subfolder.albacore_commenced:
        return   # Basecalling has already commenced on this directory
    else:
        # Commence basecalling on this directory
        subfolder.albacore_commenced = True

    # Number of gigabytes required for a given qsub command
    memory_allocation = 4 + NUM_THREADS
    # The read_fast5_basecaller is the algorithm that does the actual base calling,
    # what would be run if we just had Ubuntu.
    basecaller_command = "read_fast5_basecaller.py " \
                         "--input %s " \
                         "--worker_threads %s " \
                         "--save_path %s " \
                         "--config %s" \
                         % (subfolder.reads_dir, NUM_THREADS, subfolder.albacore_dir, CHOSEN_CONFIG)

    # These are both parsed into qsub which then determines what to do with it all.
    qsub_command = "qsub -o %s -e %s -S /bin/bash -l h_vmem=%dG" % (subfolder.albacore_qsub_output_log,
                                                                    subfolder.albacore_qsub_error_log,
                                                                    memory_allocation)

    # Put these all together into one grand command
    albacore_command = "echo \"%s\" | %s" % (basecaller_command, qsub_command)
    print(albacore_command)

    # Execute qsub command via subprocess.
    albacore_proc = subprocess.Popen(albacore_command, shell=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = albacore_proc.communicate()
    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of the array.
    print("Output of albacore command", stdout, stderr) 
    subfolder.albacore_jobid = stdout.rstrip().split()[2]
    generate_dataframe()


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


def is_job_complete(subfolder):
    """ Use qstat -u '*' to get all the jobs.
        We then split by \n, if any are equal to the job id assigned to the folder
        then set albacore_complete attribute to true"""
    if subfolder.albacore_complete:
        return  # Albacore has already been set to complete.

    # Run and get output of qstat -u command
    get_jobs_command = "qstat -u '*'"
    get_jobs_proc = subprocess.Popen(get_jobs_command, shell=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    get_jobs_output, get_jobs_stderr = get_jobs_proc.communicate()

    # job id is the first column of each row
    #  120882 0.55500 supernova  user...
    for line in get_jobs_output.split("\n")[:-1]:
        # Now split by space.
        first_column = line.split()[0]
        if first_column == subfolder.albacore_jobid:
            return  # Job still in queue

    # Otherwise we haven't found the job so it must be complete
    subfolder.albacore_complete = True


def take_a_break():
    time.sleep(60)


main()
