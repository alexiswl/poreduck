#!/usr/bin/env python3

import argparse
import tarfile
import logging
import os
import sys
from datetime import datetime
import subprocess
import shutil
import time
import gzip

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

"""
Usage: Given a folder of fast5 data create a new folder matching the zero filled filename of that folder.
Also use the nomenclature 
Fast5: reads/0000_<FlowcellID>_RNUMBER.fast5.tar.gz
Fastq: fastq/0000_<FlowcellID>_RNUMBER.fastq.gz
Sequencing_summary: sequencing_summary/0000_<FlowcellID>_RNUMBER.sequencing_summary.txt
"""


def get_args():
    parser = argparse.ArgumentParser(description="Tar up a folder of nanopore data")
    parser.add_argument('--sequencing_summary_path',
                        help="Path to sequencing_summary_path", required=True)
    parser.add_argument('--fastq_path',
                        help="Path to fastq files", required=True)
    parser.add_argument("--fast5_path",
                        help="Path to folder of fast5 files", required=True)
    parser.add_argument("--flowcellID",
                        help="Flowcell ID to add to suffix of tar file", required=True)
    parser.add_argument("--rnumber",
                        help="Add random number to folder", required=True)
    parser.add_argument("--md5_fast5",
                        help="File to append to md5sum for fast5 files", required=True)
    parser.add_argument("--md5_fastq",
                        help="File to append to md5sym for fastq files", required=True)
    parser.add_argument("--inplace", action='store_true', default=False, help="Remove folders as well")
    parser.add_argument("--overwrite", action='store_true', default=False,
                        help="Overwrite output file rather than append to it")
    parser.add_argument("--dry-run", dest='dry_run', action='store_true', default=False,
                        help="Don't actually tar anything, just output the logs")
    args = parser.parse_args()
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)
    return args


def get_output_name(args):
    # Return the appropriate name for the tar folder
    folder_basename = os.path.basename(os.path.normpath(args.fast5_path))
    output_name = [folder_basename.zfill(4)]
    if args.flowcellID:
        output_name.append(args.flowcellID)
    if args.rnumber:
        output_name.append(args.rnumber)
    output_name = '_'.join(output_name)
    logging.info("Folder name is %s" % output_name)
    return output_name


def get_fast5_output_path(fast5_path, output_name):
    # Get the .tar.gz. output to the folder
    folder_dirpath = os.path.dirname(os.path.normpath(fast5_path))
    output_path = os.path.join(folder_dirpath, output_name + ".fast5.tar.gz")
    logging.info("Tar will be written to %s" % output_path)
    return output_path


def get_fastq_output_path(fastq_path, output_name):
    # Get the fastq output to the folder
    folder_dirpath = os.path.dirname(os.path.normpath(fastq_path))
    output_path = os.path.join(folder_dirpath, "fastq", output_name + ".fastq.gz")
    logging.info("Fastq will be written to %s" % output_path)
    return output_path


def get_sequencing_summary_output_path(sequencing_summary_path, output_name):
    # Get the .tar.gz. output to the folder
    folder_dirpath = os.path.dirname(os.path.normpath(sequencing_summary_path))
    output_path = os.path.join(folder_dirpath, "sequencing_summary", output_name + ".sequencing_summary.txt")
    logging.info("Sequencing summary will be written to %s" % output_path)
    return output_path


def move_sequencing_summary_file(summary_path, output_path, overwrite=False, dry_run=False, inplace=False):
    # Move to the sequencing_summary file to new output path
    if not dry_run:
        if os.path.isfile(output_path) and not overwrite:
            logging.info("Summary file %s already exists in destination and overwrite not set. "
                         "Skipping" % output_path)
        if not inplace:
            shutil.copy(summary_path, output_path)
        else:
            shutil.move(summary_path, output_path)
    else:
        logging.info("Would have moved summary from %s into %s" % (summary_path, output_path))


def zip_and_move_fastq_file(fastq_path, output_path, overwrite=False, inplace=False, dry_run=False):
    if not dry_run:
        if os.path.isfile(output_path) and not overwrite:
            logging.info("Fastq file %s already exists in destination and overwrite not set. "
                         "Skipping" % output_path)
        # Zip and move the summary file
        with open(fastq_path, 'rb') as f_in, gzip.open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        if inplace:
            # Wait for file system to catch up then remove
            time.sleep(1)
            os.remove(fastq_path)
    else:
        logging.info("Would have gzipped and moved fastq %s into %s" % (fastq_path, output_path))


def tar_up_folder(fast5_path, output_path, overwrite=False, inplace=False, dry_run=False):
    # Tar up the folder provided
    # Get the output path
    logging.info("Output path is %s" % output_path)

    # Normalise the output path, it's come straight from args
    fast5_path = os.path.normpath(fast5_path)

    # To overwrite or not to overwrite
    if overwrite:
        file_handler_setting = "w|gz"
    else:
        # Can't actually append a compressed tar yet we're going to log and return if the file exists
        if os.path.isfile(output_path):
            logging.info("Tar file %s exists, not overwriting" % output_path)
            return
        else:
            file_handler_setting = "w|gz"  # Default is append

    # Get number of files in the path
    fast5_files = [fast5_file
                   for fast5_file in os.listdir(fast5_path)
                   if fast5_file.endswith(".fast5")]

    if len(fast5_files) == 0:
        logger.error("No fast5 files in %s" % fast5_path)
        sys.exit(1)

    if not dry_run:
        # Get time
        start_time = datetime.now()
        logging.info("Starting tarring %s into %s" % (fast5_path, output_path))
        logging.info("%d files to tar" % len(fast5_files))

        # Open up the output_path file
        archive = tarfile.open(output_path, file_handler_setting)
        # Add each of the fast5 files to the archive
        for fast5_file in fast5_files:
            input_file = os.path.join(fast5_path, fast5_file)
            output_file = os.path.join(os.path.basename(fast5_path), fast5_file)
            # Add file to archive
            archive.add(input_file, arcname=output_file)
        # If inplace also remove the input file from the system
        if inplace:
            time.sleep(3)
            # Remove folder from filesystem
            shutil.rmtree(fast5_path)
        # Close the archive
        archive.close()
        
        # Log the time taken to write the archive file
        end_time = datetime.now()
        diff_time = end_time - start_time
        logging.info("Finished tarring %s in %s" % (fast5_path, output_path))
        logging.info("Added %d files to the tar archive" % len(fast5_files))
        logging.info("Process completed in %s" % round(diff_time.total_seconds(), 2))
    else:
        logging.info("Would have tarred %s into %s" % (fast5_path, output_path))


def get_md5sum(output_path):
    logging.info("Obtaining the md5sum for %s" % output_path)

    # Grab the md5sum of the file
    md5_command = ['md5sum', output_path]
    md5_proc = subprocess.run(md5_command, stdout=subprocess.PIPE)
    md5_output = md5_proc.stdout.decode().splitlines()[0]
    logging.info("Obtained %s as md5 for %s" % (md5_output, output_path))

    # Return the md5 of the file for writing to a checksum file
    return md5_output


def write_md5sum(md5_sum, output_md5_file):
    logging.info("Writing md5sum '%s' to %s" % (md5_sum, output_md5_file))
    with open(output_md5_file, 'a') as md5_h:
        md5_h.write(md5_sum + "\n")


def main():
    # Get args
    args = get_args()
    # Set output variables
    output_name = get_output_name(args)
    output_fast5_path = get_fast5_output_path(args.fast5_path, output_name)
    output_fastq_path = get_fastq_output_path(args.fastq_path, output_name)
    output_sequencing_summary_path = get_sequencing_summary_output_path(args.fastq_path, output_name)

    # Create folders for fastq and sequencing summary files if necessary
    fastq_dir = os.path.join(os.path.dirname(os.path.normpath(output_fastq_path)), "fastq")
    sequencing_summary_dir = os.path.join(os.path.dirname(os.path.normpath(output_sequencing_summary_path)),
                                          "sequencing_summary")

    # Create fastq directory
    if not os.path.isdir(fastq_dir):
        os.mkdir(fastq_dir)

    # Create sequencing summary directory
    if not os.path.isdir(sequencing_summary_dir):
        os.mkdir(sequencing_summary_dir)

    # Tar up folder
    tar_up_folder(args.fast5_path, output_fast5_path,
                  args.overwrite, args.inplace, args.dry_run)
    # Move fastq folder
    zip_and_move_fastq_file(args.fastq_path, output_fastq_path,
                            args.overwrite, args.inplace, args.dry_run)
    # Move sequencing summary file
    move_sequencing_summary_file(args.sequencing_summary_path, output_sequencing_summary_path,
                                 args.overwrite, args.inplace, args.dry_run)

    # Get md5 for fastq and fast5
    if not args.dry_run:
        md5sum_fast5 = get_md5sum(output_fast5_path)
        md5sum_fastq = get_md5sum(output_fastq_path)

        # Write md5
        write_md5sum(md5sum_fast5, args.md5_fast5)
        write_md5sum(md5sum_fastq, args.md5_fastq)


if __name__ == "__main__":
    main()
