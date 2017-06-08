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
STATUS_CSV = ""
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

    def check_albacore_job_status(self):
        if has_commenced(self.albacore_jobid):
            self.albacore_commenced = True
            update_dataframe(get_current_subfolder(self.name))
        if has_completed(self.albacore_jobid):
            self.albacore_complete = True
            update_dataframe(get_current_subfolder(self.name))

    def check_extraction_job_status(self):
        if has_commenced(self.extracted_jobid):
            self.extracted_commenced = True
            update_dataframe(get_current_subfolder(self.name))
        if has_completed(self.extracted_jobid):
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
    merge_fastq_files()


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
4. pick_up_from_previous_run
"""


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
    parser.add_argument("--resume", type=str, required=False, default=None,
                        help="Resume the albacore run, need a csv file")
    parser.add_argument("--qsub_directory", type=str, required=False, default=None,
                        help="Where would you like to place the qsub files?")

    return parser.parse_args()


def set_global_variables(args):
    # Global variables
    global READS_DIR, ALBACORE_DIR, WORKING_DIR, NUM_THREADS, CHOSEN_CONFIG, FASTQ_DIR
    global STATUS_CSV, QSUB_LOG_DIR
    READS_DIR = args.reads_dir
    if args.output_dir is not None:
        ALBACORE_DIR = args.output_dir
    CHOSEN_CONFIG = CONFIGS[args.config]
    NUM_THREADS = args.num_threads
    if args.fastq_dir is not None:
        FASTQ_DIR = args.fastq_dir
    if args.resume is not None:
        STATUS_CSV = args.resume
    if args.qsub_directory is not None:
        QSUB_LOG_DIR = args.qsub_directory


def check_directories():
    # Make sure the directories exist, change to reads directory,
    # Create any other necessary directories for the script to run.
    global READS_DIR, ALBACORE_DIR, PARENT_DIRECTORY, QSUB_LOG_DIR, FASTQ_DIR
    global STATUS_CSV
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

    if not STATUS_CSV == "":
        if not os.path.isfile(STATUS_CSV):
            sys.exit("Resume csv specified but does not exist")
    else:
        STATUS_CSV = os.path.join(PARENT_DIRECTORY, "status.csv")
        if os.path.isfile(STATUS_CSV):
            print("Warning, STATUS_CSV not defined but still exists. Removing file %s" % STATUS_CSV)
            os.remove(STATUS_CSV)


def pick_up_from_previous_run():
    if os.path.isfile(STATUS_CSV):
        previous_dataframe = pd.read_csv(STATUS_CSV, header=True)

        # Set subfolders and all the statuses
        for index, row in previous_dataframe.iterrows():
            subfolder = row['name']
            SUBFOLDERS.append(Subfolder(subfolder))

        for subfolder in SUBFOLDERS:
            SUBFOLDERS[subfolder].extracted_submitted = row['extracted_submitted']
            SUBFOLDERS[subfolder].extracted_commenced = row['extracted_commenced']
            SUBFOLDERS[subfolder].extracted_jobid = row['extracted_jobid']
            SUBFOLDERS[subfolder].extracted_complete = row['extracted_completed']
            SUBFOLDERS[subfolder].albacore_submitted = row['albacore_submitted']
            SUBFOLDERS[subfolder].albacore_commenced = row['albacore_commenced']
            SUBFOLDERS[subfolder].albacore_jobid = row['albacore_jobid']
            SUBFOLDERS[subfolder].albacore_complete = row['albacore_complete']
            SUBFOLDERS[subfolder].folder_removed = row['folder_removed']
            SUBFOLDERS[subfolder].fastq_moved = row['fastq_moved']
            # Now make sure none of the jobs have failed and reset if so
            if has_failed(SUBFOLDERS[subfolder].albacore_jobid):
                SUBFOLDERS[subfolder].albacore_submitted = False
                SUBFOLDERS[subfolder].albacore_commenced = False
                SUBFOLDERS[subfolder].albacore_jobid = -1
                SUBFOLDERS[subfolder].albacore_complete = False
            if has_failed(SUBFOLDERS[subfolder].extracted_jobid):
                SUBFOLDERS[subfolder].extracted_submitted = False
                SUBFOLDERS[subfolder].extracted_commenced = False
                SUBFOLDERS[subfolder].extracted_jobid = -1
                SUBFOLDERS[subfolder].extracted_complete = False


"""
Pipeline functions:
1. extract_tarred_read_set
2. run_albacore
3. remove_folder
4. move_fastq_file
5. tar_albacore_folder
6. merge_fastq_files
"""


def extract_tarred_read_set(subfolder):
    # Extract the tarred reads, we currently do not delete after extraction...
    # maybe a good idea?

    # Has the dataset already been set for extraction?
    if subfolder.extracted_submitted:
        return
    else:
        subfolder.extracted_submitted = True

    tar_command = "pigz -dc %s | tar -xf -" % os.path.join(READS_DIR, subfolder.tar_filename)
    qsub_command = "qsub -o %s -e %s -S /bin/bash -wd %s -l hostname=melb-compute06" % \
                   (subfolder.extracted_qsub_output_log,
                    subfolder.extracted_qsub_error_log,
                    READS_DIR)
    tar_proc = subprocess.Popen("echo \"%s\" | %s" % (tar_command, qsub_command), shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()

    # Stdout equal to 'Your job 122079 ("STDIN") has been submitted\n'
    # So job equal to third element of the array.
    print("Output of extraction command", stdout, stderr)
    subfolder.extracted_jobid = int(stdout.rstrip().split()[2])
    update_dataframe(subfolder)


def run_albacore(subfolder):
    """ Generate read_fast5_basecaller command
        And qsub command and pipe former into latter
        This function also sets the jobid of a given folder.
    """
    if subfolder.albacore_submitted:
        return   # Basecalling has already queued for this directory
    else:
        # Commence basecalling on this directory
        subfolder.albacore_submitted = True

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
    qsub_command = "qsub " \
                   "-o %s " \
                   "-e %s " \
                   "-S /bin/bash -l hostname=melb-compute06 " \
                   "-l h_vmem=%dG " \
                   "-wd %s -v OMP_NUM_THREADS=1" \
                   % (subfolder.albacore_qsub_output_log, subfolder.albacore_qsub_error_log,
                      memory_allocation, PARENT_DIRECTORY)

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
        print("Um... the archive doesn't exist, I'm not going to delete the folder")
        return

    # Generate remove command
    remove_command = "rm -rf %s" % os.path.join(READS_DIR, subfolder.name)

    # Run remove command through subprocess
    remove_proc = subprocess.Popen(remove_command, shell=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = remove_proc.communicate()
    print("Output of deleting folder is", stdout, stderr)
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

    # Get fastq files in the workspace directory
    fastq_files = [fastq for fastq in os.listdir(subfolder.workspace_dir)
                   if fastq.endswith(".fastq")]
    if len(fastq_files) > 1:
        print("It seems like we have more than 1 file here", fastq_files)

    fastq_file_index = 0
    for fastq_file in fastq_files:
        # Rename fastq to have a .0.fastq at the start
        subfolder.fastq_file = subfolder.fastq_file.replace(".fastq", "") + "." + str(fastq_file_index) + ".fastq"
        # We will move and rename but make sure we're not overwriting anything
        if os.path.isfile(os.path.join(FASTQ_DIR, subfolder.fastq_file)):
            fastq_file_index += 1
            print("Hmmm, looks like we might accidentally overwrite something here, adding 1 to the index")
            subfolder.fastq_file = subfolder.fastq_file.replace(str(fastq_file_index-1) + ".fastq",
                                                                str(fastq_file_index) + ".fastq")

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
    update_dataframe(subfolder)


def tar_albacore_folder(subfolder):
    """
    Tar up the albacore folder into .tar.gz files.
    It's really just a bunch of metadata anyway
    """
    tar_command = "tar -cf - %s --remove_files | pigz -p 16 > %s" % \
                  (subfolder.albacore_dir, subfolder.albacore_zip_file)

    # Run tar_command through subprocess.
    tar_proc = subprocess.Popen(tar_command, shell=True,
                                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = tar_proc.communicate()
    if not stdout == "" or not stderr == "":
        print ("Output of tar albacore folder command is", stdout, stderr)
    subfolder.albacore_tarred = True


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


"""
Miscellaneous pipeline non-core functions
1. get_subfolders
2. is_still_transferring
3. take_a_break
4. generate_dataframe
5. new_directories
6. is_still_basecalling
7. get_current_subfolder
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
        if subfolder not in STATUS_DF.index.tolist():
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
        print("qacct stderr of %s", qacct_stderr)
    return True


def has_completed(job_id):
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
        print("qacct stderr of %s", qacct_stderr)
    
    if has_failed(job_id):
        print("It appears that the job %d has failed" % job_id)
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
   
    if int(qacct_stdout) > 0:
        # Job has failed
        return True
    if not qacct_stderr == "":
        print("qacct stderr of %s", qacct_stderr)
    return False


# Run the main function
main()
