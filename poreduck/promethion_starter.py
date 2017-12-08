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
    def __init__(self, filename, input_folder, is_mux=False):
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
        # Now get inside the fast5 file
        with h5py.File(self.filepath) as f:
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
                self.corrupted=True
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
           
            
    def to_series(self):
        series_columns = ["Name", "Channel", "Read", "RNumber", # From file name
                          "MuxID", "StartTime", "EndTime"]  # From fast5 file - requried for plots

        return pd.Series(data=[self.filename, self.channel, self.read, self.rnumber,
                               self.mux_id, self.pore_start_time, self.pore_end_time],
                         index=series_columns)


class Subfolder:
    def __init__(self, reads_path, number, metadata_dir, is_mux=False, threshold=4000):
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
        # Initiliase start and end times
        self.rnumber = ""
        self.start_time = None
        self.end_time = None
        self.threshold = threshold
    
    def get_new_folder_name(self): 
        # Create the new folder name
        if self.is_mux:
            mux_seq = "mux_scan"
        else:
            mux_seq = "sequencing_run"
        self.new_folder_name = '_'.join([self.standard_int, mux_seq, self.rnumber])
        self.new_folder_path = os.path.join(self.pardir, self.new_folder_name)
        self.tar_file = self.new_folder_name + ".tar.gz"
        self.metadata_path = os.path.join(self.metadata_dir, self.new_folder_name+".tsv")

    def get_fast5_files(self):
        # Fast5 class
        self.fast5_files = [Fast5file(fast5_file, self.path, is_mux=self.is_mux)
                            for fast5_file in os.listdir(self.path)
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
        self.pd.to_csv(self.metadata_path, header=True, index=False, sep="\t")        
   
    def check_if_full(self):    
        if self.is_full:
            # We shouldn't be calling this twice.
            return
        raw_fast5_count = len([fast5 
                               for fast5 in os.listdir(self.path) 
                               if fast5.endswith(".fast5")])
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
        return(start_time + minutes)

    def tar_folder(self):
        # Always ensure the previous stage has been completed
        if not self.is_full:
            return
        # Always ensure the this stage has been 
        # First move the folder to the new folder path location.
        shutil.move(self.path, self.new_folder_path)
        
        """Tar folder using pigz"""
        # When tarring we need to be in the directory, rather than use the absolute path
        oldpwd = os.getcwd()
        os.chdir(self.pardir)
        # Now declare tar and pigz command
        tar_command = ["tar", "-cf", '-', self.new_folder_name, "--remove-files"]
        pigz_command = ["pigz", "-9", "-p", "16"]
        # Pipe output of tar command into pigz command. 
        with open(self.tar_file+".tmp", 'w') as output_h:
            tar_proc = subprocess.Popen(tar_command, stdout=subprocess.PIPE)
            pigz_proc = subprocess.Popen(pigz_command, stdin=tar_proc.stdout,
                                         stdout=output_h, stderr=subprocess.PIPE)
        tar_proc.stdout.close()
        pigz_output = pigz_proc.communicate()
        print(pigz_output)
        if not pigz_proc.returncode == 0:
            print(pigz_command.stderr.decode())
        shutil.move(self.tar_file+".tmp", self.tar_file)
        os.chdir(oldpwd)
        self.is_tarred = True

        
class Run:
    def __init__(self, path, is_mux = False):
        self.path = path 
        self.fast5_path = os.path.join(self.path, "fast5")
        self.is_mux = is_mux
        self.complete = False
        self.subfolders = []
        self.metadata_dir = os.path.join(self.path, "metadata")
        if not os.path.isdir(self.metadata_dir):
            os.mkdir(self.metadata_dir)

    def get_subfolders(self):

        # Get all folders
        folders = [folder 
                   for folder in os.listdir(self.fast5_path)
                   if os.path.isdir(os.path.join(self.fast5_path, folder))
                   and folder.isdigit()
                  ]
        for folder in folders:
            # Don't read folders
            if folder in [subfolder.number for subfolder in self.subfolders]:
                continue 
            # Append new folders
            self.subfolders.append(Subfolder(self.fast5_path, folder, self.metadata_dir, is_mux=self.is_mux))

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
    def __init__(self, sample_name, samplesheet, run_dir):
        self.pd = samplesheet.query("SampleName=='%s'" % sample_name)
        self.runs = []
        self.is_running = True
        for index, run in self.pd.iterrows():  
            mux_path = os.path.join(run_dir, run.SlurmID, "reads", '_'.join([run.GrnwchMuxStartDate, run.GrnwchMuxStartTime, run.SampleName]))
            seq_path = os.path.join(run_dir, run.SlurmID, "reads", '_'.join([run.GrnwchSeqStartDate, run.GrnwchSeqStartTime, run.SampleName]))
            self.runs.append(Run(mux_path, is_mux=True))
            self.runs.append(Run(seq_path, is_mux=False))
    
    def is_run_complete(self):
        # All ports must be complete to return true.
        for run in self.runs:
            if not run.is_run_complete(): 
                return False
        return True


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
                        help="Path to tab delimited samplesheet. Columns are SampleName, GrnwchMuxStartDate, GrnwchMuxStartTime, GrnwchSeqStartDate, GrnwchSeqStartTime SlurmID")
    parser.add_argument("--pca_dir", type=str, required=True,
                        help="/path/to/PCA00XX/")
    args = parser.parse_args()
    return args                                
               
                                     
def is_still_running(samples):
    for sample in samples:
        if not sample.is_run_complete():
            return True
    return False
            

def samplesheet_to_pd(samplesheet):
    return pd.read_csv(samplesheet, header=0, sep="\t", dtype=str) #{"GrnwchMuxStartDate": str,
                                                                    # "GrnwchMuxStartTime": str,
                                                                    # "GrnwchSeqStartDate": str,
                                                                     #     "GrnwchSeqStartTime": str})


def main():
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    samples = [Sample(sample, samplesheet, args.pca_dir)
               for sample in samplesheet.SampleName.unique().tolist()]
    running = True
    while running:
        for sample in samples:
            for run in sample.runs:
                run.get_subfolders()
                run.tar_subfolders()
        running = is_still_running(samples)

main()
