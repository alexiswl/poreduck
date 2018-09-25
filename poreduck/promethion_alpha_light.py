#!/usr/bin/env python3

import argparse
import tarfile
import logging
import os
import sys
from datetime import datetime
import subprocess

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

"""
Usage: Given a folder of fast5 data create a new folder matching the zero filled filename of that folder.
Also use the nomenclature 0000_<FlowcellID>_RNUMBER.tar.gz 
"""


def get_args():
    parser = argparse.ArgumentParser(description="Tar up a folder of nanopore data")
    parser.add_argument("--folder_path", help="Path to folder of fast5 files")
    parser.add_argument("--flowcellID", help="Flowcell ID to add to suffix of tar file")
    parser.add_argument("--rnumber", help="Add random number to folder")
    parser.add_argument("--md5", help="File to append to md5sum")
    parser.add_argument("--dry-run", dest='dry_run', action='store_true', default=False,
                        help="Don't actually tar anything, just output the logs")
    args = parser.parse_args()
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)
    return args


def get_output_name(args):
    # Return the appropriate name for the tar folder
    folder_basename = os.path.basename(os.path.normpath(args.folder_path))
    output_name = [folder_basename.zfill(4)]
    if args.flowcellID:
        output_name.append(args.flowcellID)
    if args.rnumber:
        output_name.append(args.rnumber)
    output_name = '_'.join(output_name)
    logging.info("Folder name is %s" % output_name)
    return output_name


def get_output_path(folder_path, output_name):
    # Get the .tar.gz. output to the folder
    folder_dirpath = os.path.dirname(os.path.normpath(folder_path))
    output_path = os.path.join(folder_dirpath, output_name + ".tar.gz")
    logging.info("Tar will be written to %s" % output_path)
    return output_path


def tar_up_folder(folder_path, output_path, dry_run=False):
    # Tar up the folder provided
    # Get the output path
    logging.info("Output path is %s" % output_path)

    # Normalise the output path, it's come straight from args
    folder_path = os.path.normpath(folder_path)

    # Get number of files in the path
    fast5_files = [fast5_file
                   for fast5_file in os.listdir(folder_path)
                   if fast5_file.endswith(".fast5")]

    if len(fast5_files) == 0:
        logger.error("No fast5 files in %s" % folder_path)
        sys.exit(1)

    if not dry_run:
        # Get time
        start_time = datetime.now()
        logging.info("Starting tarring %s into %s" % (folder_path, output_path))
        logging.info("%d files to tar" % len(fast5_files))
        # Open up the output_path file
        archive = tarfile.open(output_path, "w|gz")
        # Add each of the fast5 files to the archive
        for fast5_file in fast5_files:
            input_file = os.path.join(folder_path, fast5_file)
            output_file = os.path.join(os.path.basename(folder_path), fast5_file)
            archive.add(input_file, arcname=output_file)
        archive.close()
        end_time = datetime.now()
        diff_time = end_time - start_time
        logging.info("Finished tarring %s in %s" % (folder_path, output_path))
        logging.info("Added %d files to the tar archive" % len(fast5_files))
        logging.info("Process completed in %s" % round(diff_time.total_seconds(), 2))
    else:
        logging.info("Would have tarred %s into %s" % (folder_path, output_path))


def get_md5sum(output_path):
    logging.info("Obtaining the md5sum for %s" % output_path)
    # Grab the md5sum of the file
    md5_command = ['md5sum', output_path]
    md5_proc = subprocess.run(md5_command, stdout=subprocess.PIPE)
    md5_output = md5_proc.stdout.decode().splitlines()[0]
    logging.info("Obtained %s as md5 for %s" % (md5_output, output_path))
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
    output_path = get_output_path(args.folder_path, output_name)
    # Tar up folder
    tar_up_folder(args.folder_path, output_path, args.dry_run)
    # Get md5
    if not args.dry_run:
        md5sum = get_md5sum(output_path)
    # Write md5
    if not args.dry_run:
        write_md5sum(md5sum, args.md5)


if __name__ == "__main__":
    main()
