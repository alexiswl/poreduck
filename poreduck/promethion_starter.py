#!/usr/bin/env python3

"""
Given a sample sheet and a run directory, tar up a PromethION run.
"""

import argparse
import h5py
from datetime import datetime, timedelta
import pandas as pd
import os
import shutil
import subprocess
import paramiko
from tempfile import NamedTemporaryFile

"""
Class types

Run:
A run exists for sample name.
Each run has a subset of the sample sheet dataframes.
We imply the read paths from each flowcell being used.

Subfolder:
Folder specified as an int value.
By default each folder holds 4000 fast5 files.
Sets up a dataframe for each fastq file


Fastq File:
Grabs meta data about each fastq file.
Along with the expected finish time of the run.

"""


class Fast5file:
    def __init__(self, filename, input_folder, open_sftp, is_mux=False):
        self.filename = filename
        self.filepath = os.path.join(input_folder, self.filename)
        # To get the rest of the attributes from the filename,
        # Split the filename into that which is before 'sequencing_run'
        # And that which is after sequencing run
        if is_mux:
            pre_seq_pivot, post_seq_pivot = self.filename.split("mux_scan", 1)
        else:
            pre_seq_pivot, post_seq_pivot = self.filename.split("sequencing_run", 1)
        pre_seq_pivot = pre_seq_pivot.split("_")
        post_seq_pivot = post_seq_pivot.split("_")
        # Slurm node, date and f_id are all premux
        # Not needed as of yet
        # Channel, read, rnumber and sample_id are all post pivot
        self.channel = post_seq_pivot[-2]
        self.read = post_seq_pivot[-4]
        self.rnumber = post_seq_pivot[-6]
        self.sample_id = post_seq_pivot[0:-6]
        self.corrupted = False
        # Download the fast5file to /tmp.
        tmp_file = NamedTemporaryFile(delete=False)  # We will delete at the end
        open_sftp.get(self.filepath, tmp_file.name)
        # Now get inside the fast5 file
        with h5py.File(tmp_file.name) as f:
            # Get attributes from /Raw/Reads/
            try:
                read_attributes = dict(f['Raw/Reads/Read_%s' % self.read].attrs.items())
            except KeyError:
                self.corrupted = True
                print("%s is corrupted" % self.filename)
                return
            # Get values from inside the fast5 value
            self.mux_id = read_attributes["start_mux"]
            self.read_id = read_attributes["read_id"]
            # Get the Experiment duration from the context tags.
            try:
                context_tags = dict(f['UniqueGlobalKey/context_tags'].attrs.items())
                tracking_id = dict(f['UniqueGlobalKey/tracking_id'].attrs.items())
                channel_id = dict(f['UniqueGlobalKey/channel_id'].attrs.items())
            except KeyError:
                self.corrupted = True
                print("%s is corrupted" % self.filename)
                return
            # Time (even mux have start times)
            self.exp_start_time = datetime.strptime(tracking_id['exp_start_time'].decode(), 
                                                    "%Y-%m-%dT%H:%M:%SZ") 
            # Minute format and in byte strings (not found in mux)
            if not is_mux:
                self.exp_duration_set = timedelta(minutes=int(context_tags['experiment_duration_set'].decode()))
            else:
                self.exp_duration_set = timedelta(minutes=10)
            # Read start time = start_time / sampling rate + exp_start_time
            self.duration_time = int(read_attributes["duration"])/int(channel_id["sampling_rate"])
            self.pore_start_time = timedelta(seconds=int(read_attributes["start_time"])/int(channel_id["sampling_rate"])) + self.exp_start_time
            self.pore_end_time = self.pore_start_time + timedelta(seconds=self.duration_time)
        # Now remove the link to the file (essentially delete it!)
        os.unlink(tmp_file.name)

    def to_series(self):
        series_columns = ["Name", "Channel", "Read", "RNumber", # From file name
                          "MuxID", "StartTime", "EndTime"]  # From fast5 file - requried for plots

        return pd.Series(data=[self.filename, self.channel, self.read, self.rnumber,
                               self.mux_id, self.pore_start_time, self.pore_end_time],
                         index=series_columns)


class Subfolder:
    def __init__(self, reads_path, number, metadata_dir, slave, is_mux=False, threshold=4000):
        # Int name of the fast5 file
        self.number = number
        self.standard_int = self.number.zfill(4)
        self.is_mux = is_mux
        self.pardir = reads_path
        self.path = os.path.join(self.pardir, number) 
        self.metadata_dir = metadata_dir
        self.metadata_path = ""
        # Initialise the stage parameters. 
        self.is_full = False
        self.is_tarred = False
        self.run_complete = False
        # Initiliase fast5_file list and dataframe
        self.fast5_files = []
        self.pd = None
        self.slave = slave
        # Initiliase start and end times
        self.rnumber = ""
        self.start_time = None
        self.end_time = None
        self.threshold = threshold
        self.new_folder_name = ""
        self.new_folder_path = ""
        self.tar_file = ""
        self.tar_path = ""
        self.num_fast5_files = 0

    def get_new_folder_name(self): 
        # Create the new folder name
        if self.is_mux:
            mux_seq = "mux_scan"
        else:
            mux_seq = "sequencing_run"
        self.new_folder_name = '_'.join([self.standard_int, mux_seq, self.rnumber])
        self.new_folder_path = os.path.join(self.pardir, self.new_folder_name)
        self.tar_file = self.new_folder_name + ".tar.gz"
        self.tar_path = os.path.join(self.pardir, self.tar_file)
        self.metadata_path = os.path.join(self.metadata_dir, self.new_folder_name+".tsv")

    def get_fast5_files(self):
        # Reopen up open_sftp
        open_sftp = self.slave.ssh_client.open_sftp()
        # Fast5 class
        self.fast5_files = [Fast5file(fast5_file, self.path, open_sftp, is_mux=self.is_mux)
                            for fast5_file in open_sftp.listdir(path=self.path)
                            if fast5_file.endswith(".fast5")]
        self.num_fast5_files = len(self.fast5_files)
    
    def get_dataframe(self):
        # Generate dataframe for subfolder
        self.pd = None
        for fast5_file in self.fast5_files:
            if fast5_file.corrupted:
                continue
            if self.pd is None:
                first_series = fast5_file.to_series()
                
                self.pd = pd.DataFrame(data=[first_series])

            else:
                self.pd = self.pd.append(fast5_file.to_series(),
                                         ignore_index=True)

    def write_dataframe(self):
        open_sftp = self.slave.ssh_client.open_sftp()
        tsv_ftpfile = open_sftp.file(self.metadata_path, 'w')
        self.pd.to_csv(tsv_ftpfile, header=True, index=False, sep="\t")
        tsv_ftpfile.flush()
        tsv_ftpfile.close()
        open_sftp.close()
   
    def check_if_full(self):    
        if self.is_full:
            # We shouldn't be calling this twice.
            return
        # Use opensftp to get the current status of the directory.
        open_sftp = self.slave.ssh_client.open_sftp()
        raw_fast5_count = len([fast5
                               for fast5 in open_sftp.listdir(path=self.path)
                               if fast5.endswith(".fast5")])
        open_sftp.close()
        # Determine if this folder is still being written to
        if raw_fast5_count < self.threshold and not self.is_mux:
            self.is_full = False
            return
        self.is_full = True
        # Get fast5 files
        self.get_fast5_files()
        # Get dataframe
        self.get_dataframe()
        # May still be full if run has stopped and last bin.
        self.run_complete = self.is_run_complete()
        # Redetermine if bin is full
        if not self.is_full and self.run_complete:
            # Get fast5 files
            self.get_fast5_files()
            # Get data frame
            self.get_dataframe()

        # Set values of rnumber, start_time and end_time
        try:
            self.rnumber = self.pd.RNumber.unique().item()
        except ValueError:
            print(self.pd.RNumber.unique())
        self.start_time = self.pd.StartTime.min()
        self.end_time = self.pd.EndTime.max()
        # Get new folder name
        self.get_new_folder_name()
        self.write_dataframe()
    
    def is_run_complete(self):
        # We cannot complete this operation is the bin is not full
        if not self.is_full:
            return False
        # Get standard fast5 file (not that simple)
        while True:
            fast5_files_iter = iter(self.fast5_files)
            fast5_file = next(fast5_files_iter)
            if not fast5_file.corrupted:
                break
        # Get expected finish time from fast5 file
        exp_finished_time = self.get_run_finish_time(fast5_file)
        # Determine if run is complete.
        current_time = datetime.utcnow()
        # If difference is less than zero, run is finished
        diff = exp_finished_time - current_time 
        if diff.total_seconds() < 0:
            return True
        else:
            return False

    def get_run_finish_time(self, fast5_file):
        start_time = fast5_file.exp_start_time
        minutes = fast5_file.exp_duration_set
        return start_time + minutes

    def tar_folder(self):
        # Always ensure the previous stage has been completed
        if not self.is_full:
            return
        # Always ensure the this stage has been 
        # First move the folder to the new folder path location through opensftp
        open_sftp = self.slave.ssh_client.open_sftp()
        open_sftp.rename(self.path, self.new_folder_path)
        open_sftp.close()

        """Tar folder using pigz"""
        # When tarring we need to be in the directory, rather than use the absolute path
        # Now declare tar and pigz command
        tar_command = ' '.join(["tar", "-cf", '-', self.new_folder_name, "--remove-files"])
        gzip_command = ' '.join(["gzip", "-", '>', self.tar_file+".tmp"])
        tar_and_gzip_command = ' | '.join([tar_command, gzip_command])
        # Use openssh client to run tar gzip through ssh
        stdin, stdout, stderr = self.slave.ssh_client.exec_command('; '.join(["cd %s" % self.pardir,
                                                                              tar_and_gzip_command]))
        print(stdin.__dict__, stdout.__dict__, stderr.__dict__)
        # Move output from .tmp to tar file
        open_sftp = self.slave.ssh_client.open_sftp()
        open_sftp.rename(self.tar_path+".tmp", self.tar_path)
        open_sftp.close()
        self.is_tarred = True

        
class Run:
    def __init__(self, path, slave, is_mux=False):
        self.path = path 
        self.fast5_path = os.path.join(self.path, "fast5")
        self.is_mux = is_mux
        self.complete = False
        self.subfolders = []
        self.metadata_dir = os.path.join(self.path, "metadata")
        # Slave object for run. Use to connect to data.
        self.slave = slave
        self.slave.connect()
        # Create metadata directory through the slave clients opensftp
        opensftp = self.slave.ssh_client.open_sftp()
        try:  # Try change to the directory, create if does not exist
            opensftp.chdir(self.metadata_dir)
        except IOError:  # Directory does not exist.
            opensftp.mkdir(self.metadata_dir)
        opensftp.close()

    def get_subfolders(self):
        # Get all folders through opensftp
        opensftp = self.slave.ssh_client.open_sftp()
        folders = opensftp.listdir(path=self.fast5_path)
        folders = [folder 
                   for folder in folders
                   if folder.isdigit()
                  ]
        for folder in folders:
            # Don't read in folders
            if folder in [subfolder.number for subfolder in self.subfolders]:
                continue 
            # Append new folders
            self.subfolders.append(Subfolder(self.fast5_path, folder, self.metadata_dir, self.slave, is_mux=self.is_mux))

    def tar_subfolders(self):
        for folder in self.subfolders:
            if folder.is_tarred:
                continue
            if not folder.is_full:
                folder.check_if_full()
            if folder.is_full and not folder.is_tarred:
                folder.tar_folder()
     
    def is_run_complete(self):
        if not self.subfolders[0].is_full:
            return False
        if self.subfolders[0].is_run_complete():
            return True


class Sample:
    def __init__(self, sample_name, samplesheet, config_pd, pca_value):
        self.pd = samplesheet.query("SampleName=='%s'" % sample_name)
        # Get the active slaves for this run.
        self.slaves = [Slave(slurm_id, config_pd, pca_value)
                       for slurm_id in self.pd.SlurmID.tolist()]
        self.runs = []
        self.is_running = True
        for index, run in self.pd.iterrows():
            slave = [slave
                     for slave in self.slaves
                     if slave.slurm_id == run.SlurmID][0]
            mux_path = os.path.join(slave.reads_path, '_'.join([run.GrnwchMuxStartDate, run.GrnwchMuxStartTime, run.SampleName]))
            seq_path = os.path.join(slave.reads_path, '_'.join([run.GrnwchSeqStartDate, run.GrnwchSeqStartTime, run.SampleName]))
            self.runs.append(Run(mux_path, slave, is_mux=True))
            self.runs.append(Run(seq_path, slave, is_mux=False))
    
    def is_run_complete(self):
        # All ports must be complete to return true.
        for run in self.runs:
            if not run.is_run_complete(): 
                return False
        return True


class Slave:
    """Use the slave node config to access the data from the master node."""
    def __init__(self, slurm_id, config_pd, pca_value):
        self.slurm_id = slurm_id
        self.ssh_ip = config_pd.query("SlurmID=='%s'" % self.slurm_id)['IP'].item()
        print(self.ssh_ip)
        self.ssh_client = None
        self.reads_path = os.path.join("/media/data/", pca_value, self.slurm_id, "reads")

    def connect(self):
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh_client.connect(self.ssh_ip, key_filename='/home/prom/.ssh/id_rsa.pub')

"""
General process:
1. Get runs.
2. Each run has a boolean attribute - running.
3. While each run is running, continue to tar up subfolders but not those that don't have 4000 reads.
"""


def get_args():
    """
    Two simple arguments.
    1. Path to PCA directory
    2. Samplesheet
    """
    parser = argparse.ArgumentParser(description="Tar up the run folders of the PromethION")
    parser.add_argument("--samplesheet", type=str, required=True,
                        help="Path to tab delimited samplesheet. "
                             "Columns are SampleName, GrnwchMuxStartDate, GrnwchMuxStartTime, "
                             "GrnwchSeqStartDate, GrnwchSeqStartTime SlurmID")
    parser.add_argument("--pca_id", type=str, required=True,
                        help="/path/to/PCA00XX/")
    parser.add_argument("--ip_config", type=str, required=True,
                        help="path/to/tab-delimited-config file. "
                             "One column of IP addresses ==> one column of slave nodes.")
    args = parser.parse_args()
    return args                                
               
                                     
def is_still_running(samples):
    for sample in samples:
        if not sample.is_run_complete():
            return True
    return False
            

def samplesheet_to_pd(samplesheet):
    return pd.read_csv(samplesheet, header=0, sep="\t", dtype=str)


def config_to_pd(config):
    return pd.read_csv(config, header=0, sep="\t", dtype=str)


def main():
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    config_pd = config_to_pd(args.ip_config)
    samples = [Sample(sample, samplesheet, config_pd, args.pca_id)
               for sample in samplesheet.SampleName.unique().tolist()]
    running = True
    while running:
        for sample in samples:
            for run in sample.runs:
                run.get_subfolders()
                run.tar_subfolders()
        running = is_still_running(samples)

main()
