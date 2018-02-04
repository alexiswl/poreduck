#!/usr/bin/env python3

"""
Find .tar.gz files given a samplesheet 30 at a time.
Wait for 30 to be removed from aarnet cloud.
Repeat until all .tar.gz files are complete
"""

import pandas as pd
import paramiko
import numpy as np


class Subfolder:
    """
     Class subfolder, useful for tracking data
    """
    def __init__(self, name, run, log_series):
        self.name = name
        self.tar_source_path = os.path.join(run.source_path, "fast5", self.name+".tar.gz")
        self.metadata_source_path = os.path.join(run.source_path, "metadata", self.name+".tsv")
        self.tar_dest_path = os.path.join(run.dest_path, "fast5", self.name+".tar.gz")
        self.metadata_dest_path = os.path.join(run.dest_path, "metadata", self.name+".tsv")

        self.log_series = log_series
        if len(self.log_df) == 0:  # Previously didn't exist in old directory
            self.tar_file_exists = self.log_series["TarFileExists"] = True
            self.tar_file_transferring = self.log_series["TarFileTransferring"] = False
            self.tar_file_transferred = self.log_series["TarFileTransferred"] = False
            self.tar_file_removed = self.log_series["TarFileRemoved"] = False
        else:
            self.tar_file_exists = self.log_series["TarFileExists"].item()
            self.tar_file_transferring = self.log_series["TarFileTransferring"].item()
            self.tar_file_transferred = self.log_series["TarFileTransferred"].item()
            self.tar_file_removed = self.log_series["TarFileRemoved"].item()

    def to_series(self):
        return self.log_series


class Run:
    """
    A sample in a given flowcell for a given timestamp
    """
    def __init__(self, source_path, dest_path, slave, log_df):
        self.source_path = source_path
        self.dest_path = dest_path
        self.slave = slave
        self.log_df = log_df
        self.subfolders = []

    def to_df(self):
        return self.log_df

    def get_files(self):
        self.subfolders.extend(get_files(self, self.log_df))


class Sample:
    """
    Each sample class has a list of runs using the samplesheet input
    """
    def __init__(self, sample_name, samplesheet, config_pd, pca_value, output_dir, log_df):
        self.pd = samplesheet.query("SampleName=='%s'" % sample_name)
        # Get the active slaves for this run.
        self.slaves = [Slave(slurm_id, config_pd, pca_value)
                       for slurm_id in self.pd.SlurmID.tolist()]
        self.runs = []
        for index, run in self.pd.iterrows():
            slave = [slave
                     for slave in self.slaves
                     if slave.slurm_id == run.SlurmID][0]
            mux_pd = self.pd.query("GrnwchMuxStartDate=='%s' & GrnwchMuxStartTime=='%s'" %(run.GrnwchMuxStartDate, run.GrnwchMuxStartTime))
            seq_pd = self.pd.query("GrnwchSeqStartDate=='%s' & GrnwchSeqStartTime=='%s'" % (run.GrnwchSeqStartDate, run.GrnwchSeqStartTime))
            mux_log_df = pd.merge(mux_pd, log_df, how='left', on=['SampleName',
                                                                  'GrnwchMuxStartDate', 'GrnwchMuxStartTime',
                                                                  'GrnwchSeqStartDate', 'GrnwchSeqStartTime'])
            seq_log_df = pd.merge(seq_pd, log_df, how='left', on=['SampleName',
                                                                  'GrnwchMuxStartDate', 'GrnwchMuxStartTime',
                                                                  'GrnwchSeqStartDate', 'GrnwchSeqStartTime'])
            mux_source_path = os.path.join(slave.reads_path, '_'.join([run.GrnwchMuxStartDate, run.GrnwchMuxStartTime, run.SampleName]))
            mux_dest_path = os.path.join(output_dir, slave.slurm_id, '_'.join([run.GrnwchMuxStartDate, run.GrnwchMuxStartTime, run.SampleName]))
            seq_source_path = os.path.join(slave.reads_path, '_'.join([run.GrnwchSeqStartDate, run.GrnwchSeqStartTime, run.SampleName]))
            seq_dest_path = os.path.join(output_dir, slave.slurm_id, '_'.join([run.GrnwchSeqStartDate, run.GrnwchSeqStartTime, run.SampleName]))
            self.runs.append(Run(mux_path, mux_source_path, mux_dest_path, slave, mux_log_df))
            self.runs.append(Run(seq_path, seq_source_path, seq_dest_path, slave, seq_log_df))


def get_args():
    # Arguments include getting samples.
    parser = argparse.ArgumentParser(description="Transfer to webdav server")
    parser.add_argument("--samplesheet", type=str, required=True,
                        help="Path to tab delimited samplesheet. "
                             "Columns are SampleName, GrnwchMuxStartDate, GrnwchMuxStartTime, "
                             "GrnwchSeqStartDate, GrnwchSeqStartTime SlurmID")
    parser.add_argument("--pca_id", type=str, required=True,
                        help="/path/to/PCA00XX/")
    parser.add_argument("--ip_config", type=str, required=True,
                        help="path/to/tab-delimited-config file. "
                             "One column of IP addresses ==> one column of slave nodes.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Probably something like ~/aarnet/PromethION")
    args = parser.parse_args()
    return args


def samplesheet_to_pd(samplesheet):
    return pd.read_csv(samplesheet, header=0, sep="\t", dtype=str, comment='#')


def config_to_pd(config):
    return pd.read_csv(config, header=0, sep="\t", dtype=str, comment='#')


def get_prev_transfer_log(samplesheet, log_file):
    if not os.path.isfile(log_file):
        return pd.DataFrame(data=None, columns=["SampleName", "GrnwchMuxStartDate", "GrnwchMuxStartTime",
                                                "GrnwchSeqStartDate", "GrnwchSeqStartTime",
                                                "SlurmID", "TarName",
                                                "TarFileExists", "TarFileTransferring",
                                                "TarFileTransferred", "TarFileRemoved"],
                            dtype={'TarFileExists': np.bool, 'TarFileTransferring': np.bool,
                                   'TarFileTransferred': np.bool, 'TarFileRemoved': np.bool})

    return pd.read_csv(log_file, header=0, sep="\t", comment='#',
                       dtype={'TarFileExists': np.bool, 'TarFileTransferring': np.bool,
                              'TarFileTransferred': np.bool, 'TarFileRemoved': np.bool})


def get_files(run, log_df):
    """
    For each file we will return the list of type Subfolder.
    We will use to track if each Subfolder has been transferred.
    :param run:
    :return:
    """
    pm = paramiko.SSHClient()
    pm.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    pm.connect(self.ssh_ip, key_filename='/home/prom/.ssh/id_rsa.pub')
    open_sftp = pm.open_sftp()
    tarred = sorted([tar.replace(".tar.gz", "") for tar in open_sftp.listdir(path=os.path.join(run.path, "fast5"))
                     if tar.endswith(".tar.gz")])
    subfolders = []
    # Append to the class subfolders
    for tar in tarred:
        # Ensure not already in there
        if tar in [subfolder.name for subfolder in run.subfolders]:
            continue
        else:
            subfolders.append(Subfolder(tarred, run, log_df.query("TarName=='%s'" % tarred).squeeze()))
    return subfolders


def get_samples(samplesheet, config_pd, pca_id):
    """
    Grab the list of samples from the sample sheet
    :param samplesheet:
    :return: List of class Sample
    """
    samples = [Sample(sample, samplesheet, config_pd, pca_id)
               for sample in samplesheet.SampleName.unique().tolist()]
    return samples


def run_rsync(samples):
    """
    We only transfer 30 tar.gz and .tsv at a time.
    :param samples:
    :return: None
    """
    rsync_pre_command = ["rsync", "-zr", "--prune-empty-dirs", '--include="*/"']
    rsync_command = rsync_pre_command.copy()
    # We don't want to transfer more than 30 at a time
    counter = 0
    for sample in samples:
        for run in sample.runs:
            for subfolder in run.subfolders:
                if counter > 30:
                    pass
                elif subfolder.tar_file_exists and not subfolder.tar_file_transferring and not subfolder.tar_file_transferred:
                    counter += 1
                    rsync_command.append("--include='%s'" % subfolder.name+".tar.gz")
                    rsync_command.append("--include='%s'" % subfolder.name+".tsv")
                    subfolder.tar_file_transferring = True
    rsync_command.append("--exclude='*'")
    rsync_proc = subprocess.run(rsync_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if not rsync_proc.returncode == 0:
        print("Error code: %s in %s stderr is %s" % (rsync_proc.returncode, rsync_command, rsync_proc.stderr.decode()))
        sys.exit("Error code: %s in %s stderr is %s" % (rsync_proc.returncode, rsync_command, rsync_proc.stderr.decode()))
    for sample in samples:
        for run in sample.runs:
            for subfolder in run.subfolders:
                if subfolder.tar_file_exists and subfolder.tar_file_transferring:
                    subfolder.tar_file_transferring = False
                    subfolder.tar_file_transferred = True


def wait_for_file_download_to_complete():
    files_exist = True
    while files_exist:
        files_exist = False
        for sample in samples:
            for run in sample.runs:
                for subfolder in run.subfolders:
                    if subfolder.tar_file_transferring:
                        if os.path.isfile(subfolder.tar_dest_path) or os.path.isfile(subfolder.metadata_dest_path):
                            files_exist = True
                        else:
                            subfolder.tar_file_removed = True


def main():
    args = get_args()
    samplesheet = samplesheet_to_pd(args.samplesheet)
    config = config_to_pd(args.config)
    get_prev_transfer_log(samplesheet, os.path.join(args.outdir, "prom_transfers.log"))
    samples = get_samples(samplesheet, config, args.pca_id)
    # Generate the while loop:
    transferring = True
    while transferring:
        transferring = False
        run_rsync()
        wait_for_file_download_to_complete()
        # Check for more runs to do.
        for sample in samples:
            for run in sample.runs:
                for subfolder in run.subfolders:
                    if not subfolder.tar_file_transferred:
                        transferring = False



