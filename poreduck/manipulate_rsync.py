#!/usr/bin/env python3

"""
Python is probably overkill for this step.
However I believe it is important to document this code well as this is a rather serious script
"""

from fileinput import FileInput
import re
import sys
import pandas as pd
import argparse

"""
The file we want to edit is /usr/bin/sync-read-files.
This script will need to be activated everytime the PromethION software is updated.
The only goal of this script is to change
what rsync carries across. Instead of .fast5 we want .fast5.tar.gz.
We also want to add the --prune-dirs script so that we don't have any blank folders coming with us.
"""

"""
# Steps required.
1. Edit /etc/fstab through sudo and add extra line mounting /etc to /sda1 instead so it is writable.
2. Edit /usr/bin/sync-read-files and move to /tmp/new-sync-read-files
3. Kill all current rsync scripts that are running
4. Edit /etc/cron.d/sync-read-files to include .fast5.tar.gz files and .tsv files
"""

bin_file_path = "/usr/bin/sync-read-files"

print("Manipulating the file %s" % bin_file_path)

# Use FileInput to manipulate the path
try:
    with FileInput(files=(bin_file_path), inplace=True) as bin_file_h:
        for line in bin_file_h:
            line = re.sub("--include '*.fast5'", "--include '*.fast5.tar.gz'", line)
            line = re.sub("rsync -rzt", "rsync -rzt --prune-empty-dirs", line)
            print(line)
except IOError as e:
    if e[0] == errno.EPERM:
        sys.exit("Error, you need to run this command as root")
    else:
        sys.exit("Unknown IOError")

print("Completed manipulation on main node. Please check path")
print("We also need to finish the manipulation on the other slave nodes")

print("First for each slave we change group permissions on the file")
paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh_client.connect(self.ssh_ip, key_filename='/home/prom/.ssh/id_rsa.pub')
        ssh_client.exec_command("sudo chmod 775 %s" % bin_file_path)

print("Then we copy up this file")


