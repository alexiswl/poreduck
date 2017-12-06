#!/usr/bin/env python3

"""
Given a sample sheet and a run directory, tar up a PromethION run.
"""

import argparse
import hdf5


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


To do list:
1. Get values from inside the fast5 file

"""

class Fast5file:
    def __init__(self, filename, input_folder, is_mux=False):
        self.filename = filename
        self.filepath = os.path.join(filename, self.filename)
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
        # Get values from inside the fast5 value
        # mod time.
    def to_series(self):
        series_columns = ["Name", "Channel", "Read", "RNumber" # From file name
                          "MuxID", "DurationTime", "ModTime"]  # From fast5 file - requried for plots

        return pd.series(data=[self.filename, self.channel, self.read, self.rnumber,
                               self.mux_id, self.duration_time,
                               self.exp_start_time, self.exp_duration_set],
                         index=series_columns)


class Subfolder:
    def __init__(self, reads_path, number, is_mux=False, threshold=4000):
        # Int name of the fast5 file
        self.number = int(number)
        self.standard_int = f"{int(my_integer):04}"
        self.ismux = is_mux
        self.pardir = reads_path
        self.path = os.path.join(self.pardir, number)
        # Fast5 class
        self.fast5_files = [Fast5file(fast5_file, self.path, is_mux=is_mux)
                            for fast5_file in os.listdir(self.path)
                            if fast5_file.endswith(".fast5")]
        self.num_fast5_files = len(self.fast5_files)

        # Generate dataframe for subfolder
        self.pd = None
        for fast5_file in self.fast5_files:
            if self.pd is None:
                self.pd = pd.DataFrame(data=fast5_file.to_series())
            else:
                self.pd = self.pd.append(fast5_file.to_series(),
                                         ignore_index=True)
        # Get expected finish time
        exp_finished_time = get_expected_finish_time(self.fast5_files[0])

        get_current_time = ##
        is_run_finished = get_current_time - exp_finished_time
        # Determine if this folder is still being written to
        if self.num_fast5_files < threshold or not is_run_finished:
            self.complete = False
            return
        self.complete = True
        self.rnumber = self.pd.Rnumber.unique().item()
        self.start_time = self.pd.ModTime.min()
        self.end_time = self.pd.ModTime.max()
        # Create the new folder name
        if is_mux:
            mux_seq = "mux_scan"
        else:
            mux_seq = "sequencing_run"
        self.new_folder_name = '_'.join([self.standard_int, mux_seq, self.rnumber])
        self.new_folder_path = os.path.join(self.pardir, self.new_folder_name)
        self.tar_file = self.new_folder_name + ".tar.gz"


    def tar_folder(self):
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



class Run:
    def __init__(self, date, start_time, path, is_mux = False):
        self.start_date = date
        self.start_time = start_time
        self.path = path
        self.is_mux = is_mux
        self.subfolders = []

    def get_subfolders(self):
        for folder in os.listdir(self.path):
            if folder in [subfolder.number for subfolder in self.subfolders]:
                continue
            self.subfolders.append(Subfolder(self, self.path, self.number, is_mux = False))

    def tar_subfolders(self):
        for folder in self.subfolders:
            if folder.complete:
                folder.tar_folder()


class Sample:
    def __init__(self, sample_name, samplesheet, run_dir):
        self.pd = samplesheet.query("SampleName=='%s'" % sample_name)
        self.runs = []
        for run in self.pd.iterrows():
            mux_path = os.path.join(run_dir, run.slurm_id, '_'.join([run.Date, run.GrnwchMuxStartTime, run.SampleName]))
            seq_path = os.path.join(run_dir, run.slurm_id, '_'.join([run.Date, run.GrnwchSeqStartTime, run.SampleName]))
            self.runs.append(Run(run.date, run.greenwich_mux_start_time, run.path, is_mux=True))
            self.runs.append(Run(run.date, run.greenwich_seq_start_time, run.path, is_mux=False))



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


def samplesheet_to_pd():


def main():
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    samples = [Sample(sample_name, samplesheet, args.data_dir)
               for sample in samplesheet.SampleName.unique().tolist()]
    while running:
        for sample in samples:
            for run in sample.runs:
                run.get_subfolders()
                run.tar_subfolders()


main()