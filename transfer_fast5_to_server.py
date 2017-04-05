#!/usr/bin/env python

"""

This script is designed to neatly handle large amounts of MinION data.
At first MinION fast5 file structure was quite cool and a novel experience,
however the yield from a set MinION run has sky-rocketed beyond even what Clive Brown expected.
We now have the predicament that one read for each fasta file is simply too much for a computer/server to handle.

ONT have taken the step in MinKNOW version 1.4+ to reduce some problems by restricting the number of files in
each folder to 4000. This is done by creating sub-folders in the 'reads' directory, 0, 1, 2, as needed.

However this comes with a couple of bugs and a couple more issues. Any scripts that relied on all reads
being in one folder now have implement a recursive stage, and the number of files isn't strictly 4000.
Why? Because 'mux-reads' don't seem to count (so folder 0 will often have around 6000-7000 reads), and
if a run is restarted then you will end up with at least 8000 reads in each folder.

This script is designed to:
1. Detect amy 'completed' folders, those that have 4000 reads in them.
   1a. Rename this folder specific to this run so it won't be accidentally overwritten.
2. Tar up, check integrity and then md5sum this folder. We will use tiam5.sh to do this from my boubs repo.
3. Rsync the tar.gz file over to the server, and remove the source file to save space on the computer.

Rsync script will eventually timeout, which will most likely mean that the run is complete.
"""

# Import necessary modules,
import os  # list directories and path checking
import sys  # For stopping in the event of errors.
import subprocess  # Running rsync and tar functions.
import argparse  # Allow users to set commandline arguments and show help
import getpass  # Prompts user for password, just a oneoff to runt he script.
from pexpect import pxssh  # Connecting via ssh to make sure that the parent of the destination folder is there.
import time  # For snoozing
import pandas as pd  # Create dataframe of list of files with attributes for each.

# Set global variables that aren't actually global,
# Just easier than piping them into everything.
READS_DIR = ""
SERVER_NAME = ""
SERVER_USERNAME = ""
PASSWORD = ""
DEST_DIRECTORY = ""
TIMEOUT = 1800
CSV_DIR = ""


def main():
    # Get argument list, password for server and set directories.
    args = get_arguments()
    set_global_variables(args)
    check_directories()

    # Get rsync up and running - shouldn't be anything to sync
    run_rsync_command()
    rsync_running = True

    # While loop to continue tarring up folders
    while rsync_running:
        # Get list of sub-directories
        subdirs = get_subdirs()

        # May not be any new folders
        if len(subdirs) == 0:
            have_a_break()
            pass

        # For any new folders.
        for subdir in subdirs:
            # Is folder finished?
            folder_status = check_folder_status(subdir)
            if folder_status == "still writing":
                continue

            # Tar up folder(s)
            tar_folders(subdir)

    # Now we need to tar up the last folder.


def get_arguments():
    parser = argparse.ArgumentParser(
        description="The transfer_fast5_to_server transfers MiNION data from a laptop in realtime." +
                    "The process will finish when rsync times out. Use the --timeout option to adjust this.")
    parser.add_argument("--reads_dir", type=str, required=True,
                        help="/path/to/reads, should have a bunch of subfolders from 0 to N")
    parser.add_argument("--server_name", type=str, required=True,
                        help="If you were to ssh username@server, please type in the server bit.")
    parser.add_argument("--user_name", type=str, required=True,
                        help="If you were to ssh username@server, please type in the username bit")
    parser.add_argument("--dest_directory", type=str, required=True,
                        help="Where abouts on the server do you wish to place these files?")
    parser.add_argument("--timeout", type=int, required=False, default=1800,
                        help="How long did you wish to wait before rsync times out. Default is 30 mins")
    return parser.parse_args()


def set_global_variables(args):
    global READS_DIR, SERVER_NAME, SERVER_USERNAME, PASSWORD, DEST_DIRECTORY, TIMEOUT
    READS_DIR = args.reads_dir
    SERVER_NAME = args.server_name
    SERVER_USERNAME = args.user_name
    PASSWORD = get_password()
    DEST_DIRECTORY = args.dest_directory
    TIMEOUT = args.timeout


def get_password():
    return getpass.getpass('password: ')


def check_directories():
    global READS_DIR, CSV_DIR
    # Check if reads directory exists
    if not os.path.isdir(READS_DIR):
        sys.exit("Error, reads directory, %s, does not exist" % READS_DIR)
    READS_DIR = os.path.abspath(READS_DIR) + "/"
    # We will now continue working from the reads directory.
    os.chdir(READS_DIR)

    # Create CSV directory if it does not exist
    CSV_DIR = READS_DIR + "csv"
    if not os.path.isdir(CSV_DIR):
        os.mkdir(CSV_DIR)

    # Check if server is active using the ping command.
    ping_command = subprocess.Popen("ping -c 1 %s" % SERVER_NAME, shell=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ping_command.wait()
    out, error = ping_command.communicate()
    print(out, error)

    # Check if folder on server is present.
    # Log into server, then check for folder.
    dest_parent = '/'.join(DEST_DIRECTORY.split("/")[:-1])
    print(dest_parent)

    s = pxssh.pxssh()  
    
    if not s.login(SERVER_NAME, SERVER_USERNAME, PASSWORD):
	print("SSH failed on login")
    else:
	print("SSH passed")
      
    s.sendline('if [ -d %s ]; then echo "PRESENT"; fi' % dest_parent)  # Command to check if folder is there.
    s.prompt()  # match the prompt
    output = s.before  # Gets the `output of the sendline command
    
    print(output.rstrip().split('\n')[-1]) 
    if not output.rstrip().split('\n')[-1] == "PRESENT":
        # Parent folder is not present. Exit.
        sys.exit("Error, parent directory of %s does not exist" % DEST_DIRECTORY)


def run_rsync_command():
    # Create password file:
    pwd_filename = "~/.minion_transfer_pwd.txt"
    subprocess.call("touch %s" % pwd_filename, shell=True)  # Create the password
    subprocess.call("chmod 600 %s" % pwd_filename, shell=True)  # Reduce the permissions so only you can read/write to it.
    subprocess.call("echo %s > %s" % (PASSWORD, pwd_filename), shell=True)

    # Generate list of rsync options to be used.
    rsync_command_options = []
    # Free space on the laptop by deleting files that have been transferred
    rsync_command_options.append("--remove-source-files")
    rsync_command_options.append("--timeout=%d" % TIMEOUT)
    rsync_command_options.append("--password-file=%s" % pwd_filename)

    # Using the 'rsync [OPTION]... SRC [SRC]... [USER@]HOST:DEST' permutation of the command
    rsync_command = "rsync %s %s*.tar.gz %s@%s:%s" % (' '.join(rsync_command_options),
                                                      READS_DIR,
                                                      SERVER_USERNAME,
                                                      SERVER_NAME,
                                                      DEST_DIRECTORY)
    subprocess.Popen(rsync_command, shell=True)


def get_subdirs():
    subdirs = [READS_DIR + directory + "/" for directory in os.listdir(READS_DIR)
               if os.path.isdir(directory)  # Make sure that the subdirectory is a directory
               and not directory == "tmp"   # And not the tmp directory
               and not directory == "csv"]  # And not our csv directory that we'll create.

    # Remove those that are not ints.
    for subdir in subdirs:
        try:
            int(subdir)
            return True
        except ValueError:
            subdirs.remove(subdir)

    # Sort subdirectories by writing time.
    return sorted(subdirs, key=os.path.getmtime)


def have_a_break():
    time.sleep(60)


def check_folder_status(subdir):
    # Returns the folder status, also generates a little csv file for each of the corresponding folders
    # for you take home.
    return_status = "still writing"

    os.chdir(subdir)
    fast5_files = [fast5_file for fast5_file in os.listdir(subdir)
                   if fast5_file.endswith(".fast5")]

    # Create pandas data frame with each fast5 file as a row.
    # Final columns will include:
    #       fast5 file name       - name of fast5 file
    #       Modification time - time of modification of file, useful for deciding if run has finished.
    #       rnumber tag       - unique to each run, so can pull out different runs in a file.
    #       mux scan          - Is the run part of the mux scan or sequencing run? (YES/NO)
    #       channel           - Channel ID of the run.
    #       read number       - What number read is this.

    fast5_pd = pd.DataFrame(index=fast5_files, columns=['mtime', 'rnumber', 'mux', 'channel', 'read_no'])
    fast5_pd['mtime'] = [os.path.getmtime(fast5_file) for fast5_file in fast5_files]
    fast5_pd['rnumber'] = [fast5_file.split('_')[-4] for fast5_file in fast5_files]
    fast5_pd['mux'] = ["YES" if "mux_scan" in fast5_file else "NO" for fast5_file in fast5_files]
    fast5_pd['channel'] = [fast5_file.split('_')[-3] for fast5_file in fast5_files]
    fast5_pd['read_no'] = [fast5_file.split('_')[-2] for fast5_file in fast5_files]

    # Get list of runs in the folder.
    runs = fast5_pd['rnumber'].unique().tolist()

    for run in runs:
        # Before moving the mux files we need to make sure that there is some sequencing run files in the folder
        if len(fast5_pd.loc[fast5_pd.rnumber == run and fast5_pd.mux == "NO"]) == 0:
            continue  # No sequencing run files in the folder, skipping folder.

        # Move mux scan files for a given run
        fast5_to_move = fast5_pd.loc[fast5_pd.rnumber == run and fast5_pd.mux == "YES"]
        if len(fast5_to_move) == 0:
            pass
        else:
            is_mux = True
            return_status = "moving files"
            move_fast5_files(subdir, fast5_to_move, run, is_mux)
            #fast5_pd.to_csv(CSV_DIR + subdir + rnumber)

        # Move standard sequencing run files for a given run
        fast5_to_move = fast5_pd.loc[fast5_pd.rnumber == run and fast5_pd.mux == "NO"]

        # Before moving the sequencing run files, we need to make sure that this is == 4000
        if len(fast5_to_move) != 4000:
            continue  # Folder is not full
        else:
            return_status = "moving_files"
            is_mux = False
            move_fast5_files(subdir, fast5_to_move, run, is_mux)

    return return_status


def move_fast5_files(subdir, fast5_files, run, is_mux):
    if len(fast5_files) == 0:
        return
    mux = ""
    if is_mux:
        mux = "_mux_scan"

    # Create a folder in the reads directory.
    subdir = '/'.join(subdir.split('/')[:-1])  # Take off the last slash
    new_dir = subdir + "_" + run + mux
    os.mkdir(new_dir)

    for fast5_file in fast5_files:
        subprocess.call("mv %s %s" % (fast5_file, new_dir))


def tar_folders(subdir_prefix):
    # Get the subdirectories that start with the initial subdirectory.
    # So 0 may now be 0_12345 where 12345 is the rnumber.
    subdirs = [subdir for subdir in get_subdirs()
               if subdir.startswith(subdir_prefix + "_")]

    # Now tar up each folder individually
    for subdir in subdirs:
        tar_command = "tar -cf - %s | pigz -9 -p 32 > %s.tar.gz" % (subdir, subdir)
        tar_proc = subprocess.Popen(tar_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = tar_proc.communicate()
        if stderr is not None:
            print(stderr)


main()
