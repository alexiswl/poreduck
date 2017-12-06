#!/usr/bin/env python3

"""
Given a sample sheet and a run directory, tar up a PromethION run.
"""

import argparse
import hdf5
from datetime import datetime, timedelta

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
        # Now get inside the fast5 file
        with h5py.File(self.filepath) as f:
            # Get attributes from /Raw/Reads/
            read_attributes = dict(f['Raw/Reads/Read_%s' % self.read].attrs.items())
            # Get values from inside the fast5 value
            self.mux_id = read_attributes["start_mux"]
            self.read_id = read_attributes["read_id"]
            # Get the Experiment duration from the context tags.
            if is_mux:
                return
            context_tags = dict(f['UniqueGlobalKey/context_tags'].attrs.items())
            tracking_id = dict(f['UniqueGlobalKey/tracking_id'].attrs.items())
            channel_id = dict(f['UniqueGlobalKey/channel_id'].attrs.items())
            # Time format
            self.exp_start_time = datetime.datetime.strptime(tracking_id['exp_start_time'].decode(), 
                                                             "%Y-%m-%dT%H:%M:%SZ") 
            # Minute format and in byte strings.
            self.exp_duration_set = timedelta(minutes=context_tags['experiment_duration_set'].decode()))
            # Read start time = start_time / sampling rate + exp_start_time
            self.duration_time = read_attributes["duration"]/channel_id["sampling_rate"]
            self.pore_start_time = timedelta(seconds=read_attributes["start_time"]/channel_id["sampling_rate"]) + self.exp_start_time
            self.pore_end_time = self.pore_start_time + timedelta(seconds=self.duration_time)
           
            
    def to_series(self):
        series_columns = ["Name", "Channel", "Read", "RNumber" # From file name
                          "MuxID", "StartTime", "EndTime"]  # From fast5 file - requried for plots

        return pd.series(data=[self.filename, self.channel, self.read, self.rnumber,
                               self.mux_id, self.pore_start_time, self.pore_end_time],
                         index=series_columns)


class Subfolder:
    def __init__(self, reads_path, number, is_mux=False, threshold=4000):
        # Int name of the fast5 file
        self.number = int(number)
        self.standard_int = f"{int(my_integer):04}"
        self.ismux = is_mux
        self.pardir = reads_path
        self.path = os.path.join(self.pardir, number)
        self.fast5_path = os.path.join(self.path, "fast5")
        # Initialise the stage parameters. 
        self.is_full = False
        self.is_tarred = False
        self.run_complete = False
        # Initiliase fast5_file list and dataframe
        self.fast5_files = []
        self.pd = None
        # Initiliase start and end times
        self.rnumber = 0
        self.start_time = None
        self.end_time = None
        self.threshold = threshold
        
        # Create the new folder name
        if is_mux:
            mux_seq = "mux_scan"
        else:
            mux_seq = "sequencing_run"
        self.new_folder_name = '_'.join([self.standard_int, mux_seq, self.rnumber])
        self.new_folder_path = os.path.join(self.pardir, self.new_folder_name)
        self.tar_file = self.new_folder_name + ".tar.gz"
        
    def get_fast5_files(self):
        # Fast5 class
        self.fast5_files = [Fast5file(fast5_file, self.fast5_path, is_mux=is_mux)
                            for fast5_file in os.listdir(self.fast5_path)
                            if fast5_file.endswith(".fast5")]
        self.num_fast5_files = len(self.fast5_files)
    
    def get_dataframe(self):
        # Generate dataframe for subfolder
        self.pd = None
        for fast5_file in self.fast5_files:
            if self.pd is None:
                self.pd = pd.DataFrame(data=fast5_file.to_series())
            else:
                self.pd = self.pd.append(fast5_file.to_series(),
                                         ignore_index=True)
                
    def check_if_full(self):    
        if self.is_full:
            # We shouldn't be calling this twice.
            return
        self.run_complete = self.is_run_complete()
        # Determine if this folder is still being written to
        if self.num_fast5_files < self.threshold or not self.run_complete:
            self.is_full = False
            return
        self.is_full = True
        # Get fast5 files
        self.get_fast5_files()
        # Get dataframe
        self.get_dataframe()
        # Set values of rnumber, start_time and end_time
        self.rnumber = self.pd.Rnumber.unique().item()
        self.start_time = self.pd.StartTime.min()
        self.end_time = self.pd.EndTime.max()
    
    def is_run_complete(self):
        # Get expected finish time
        exp_finished_time = self.get_run_finish_time(self.fast5_files[0])
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
        tar_command = ["tar", "-cf", "-", self.path, "--remove-files"]
        pigz_command = ["pigz", "-9", "-p", "16", ">", self.tar_file]
        # Pipe output of tar command into pigz command.
        tar_proc = subprocess.run(tar_command, stdout=subprocess.PIPE)
        pigz_command = subprocess.run(pigz_command, stdin=tar_proc.stdout,
                                      stderr=subprocess.PIPE)
        tar_proc.stdout.close()
        if not pigz_command.returncode == 0:
            print(pigz_command.stderr.decode())
        os.chdir(oldpwd)
        self.is_tarred = True

        
class Run:
    def __init__(self, date, start_time, path, is_mux = False):
        self.start_date = date
        self.start_time = start_time
        self.path = path
        self.is_mux = is_mux
        self.complete = False
        self.subfolders = []

    def get_subfolders(self):
        for folder in os.listdir(self.path):
            if folder in [subfolder.number for subfolder in self.subfolders]:
                continue
            self.subfolders.append(Subfolder(self, self.path, self.number, is_mux = False))

    def tar_subfolders(self):
        for folder in self.subfolders:
            if folder.is_tarred:
                continue
            if not folder.is_full:
                folder.check_if_full()
            if folder.is_full and not folder.is_tarred:
                folder.tar_folder()
     
    def is_run_complete(self):
        if self.subfolders[0].is_run_complete():
            return True


class Sample:
    def __init__(self, sample_name, samplesheet, run_dir):
        self.pd = samplesheet.query("SampleName=='%s'" % sample_name)
        self.runs = []
        self.is_running = True
        for run in self.pd.iterrows():
            mux_path = os.path.join(run_dir, run.slurm_id, '_'.join([run.GrnwchMuxStartDate, run.GrnwchMuxStartTime, run.SampleName]))
            seq_path = os.path.join(run_dir, run.slurm_id, '_'.join([run.GrnwchSeqStartDate, run.GrnwchSeqStartTime, run.SampleName]))
            self.runs.append(Run(run.date, run.greenwich_mux_start_time, run.mux_path, is_mux=True))
            self.runs.append(Run(run.date, run.greenwich_seq_start_time, run.seq_path, is_mux=False))
    
    def is_still_running(self):
        # All ports must be complete to return true.
        for run in self.runs():
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
    parser = argparse.ArgumentParser(description="Tar up the run folders of the PromethION)
    parser.add_argument("--samplesheet", type=str, required=True,
                        help="Path to tab delimited samplesheet. Columns are SampleName, GrnwchMuxStartDate, GrnwchMuxStartTime, GrnwchSeqStartDate, GrnwchSeqStartTime")
    parser.add_argument("--pca_dir", type=str, required=True,
                        help="/path/to/PCA00XX/")
    args = parser.parse_args()
    return args                                
               
                                     
def is_still_running(samples):
    for sample in samples:
        if sample.is_still_running():
            return True
    return False
            

def samplesheet_to_pd(samplesheet):
    return pd.read_csv(samplesheet, header=0, sep="\t")


def main():
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    samples = [Sample(sample_name, samplesheet, args.data_dir)
               for sample in samplesheet.SampleName.unique().tolist()]
    running = True
    while running:
        for sample in samples:
            for run in sample.runs:
                run.get_subfolders()
                run.tar_subfolders()
        running = is_still_running(samples)

main()
