#!/usr/bin/env python3

import os
import sys
import argparse
import poreduck.version


# Return version number
def get_version():
    return poreduck.__version__


# Run the function we select through argparse
def run_function(args):
    # Which function of the three did we choose?
    if args.command == "tarMyFast5":
        import poreduck.minion_starter as command_to_run
    if args.command == "albacoreHpc":
        import poreduck.albacore_server_scaled as command_to_run
    if args.command == "compare_runs":
        import poreduck.compare_runs as command_to_run
    # Now run it!
    command_to_run.main(args)


# Globals for albacore server
FLOWCELLS = ["FLO-MIN107", "FLO-MIN106", "FLO-PRO001"]

KITS = ["SQK-LWP001", "SQK-NSK007", "VSK-VBK001", "SQK-RAS201", "SQK-RBK001", "SQK-LWB001",
        "SQK-RNA001", "SQK-RLI001", "SQK-RAD002", "SQK-RLB001", "SQK-RAB201", "SQK-LSK208",
        "SQK-LSK108", "SQK-RAD003", "SQK-DCS108", "SQK-PCS108", "SQK-LSK308"]  # More kits to come

QSUB_TYPES = ["SLURM"]

# Define main script:


def main():
    # Create poreduck parser
    parser = argparse.ArgumentParser(prog='poreduck', description="poreduck package")
    parser.add_argument("--version", help="Get version of poreduck",
                        action="version",
                        version=poreduck.__version__)

    subparsers = parser.add_subparsers(help="Callable poreduck functions", dest="command")

    # Transfer to server arguments
    tar_parser = subparsers.add_parser('tarMyFast5',
                                       help="The transfer_fast5_to_server transfers MiNION data from a " +
                                                 "laptop in realtime. The process will finish when the script " +
                                                 "believes MinKNOW is no longer running.")
    tar_parser.add_argument("--samplesheet", type=str, required=True,
                        help="Path to tab delimited samplesheet. "
                             "Columns are SampleName, UTCMuxStartDate, UTCMuxStartTime, "
                             "UTCSeqStartDate, UTCSeqStartTime SlurmID")
    tar_parser.add_argument("--reads_path", type=str, required=True,
                        help="path/to/reads. "
                             "Ubuntu: /var/lib/MinKNOW/data/reads"
                             "Mac: /Library/MinKNOW/data"
                             "Windows: C:\\data\\reads")
    tar_parser.set_defaults(func=run_function)

    # Albacore server arguments:
    albacore_parser = subparsers.add_parser('albacoreHPC',
                                            help="The albacore-server scaled command incorporates qsub to spread " +
                                                 "the server load and rapidly " + "generate data off the MinION.")
    albacore_parser.add_argument("--fast5_dir", type=str, required=True,
                                 help="/path/to/reads/fast5, " +
                                      "should have a bunch of tar zipped files in it.")
    albacore_parser.add_argument("--flowcell", type=str, choices=FLOWCELLS, required=True,
                                 help="Pick a flowcell version")
    albacore_parser.add_argument("--kit", type=str, choices=KITS, required=True,
                                 help="Pick a kit type")
    albacore_parser.add_argument("--output_dir", type=str, required=False, default=None,
                                 help="Will be called 'albacore'" +
                                      " and sit adjacent to the reads folder if left blank.")
    albacore_parser.add_argument("--num_threads", type=int, required=False, default=5,
                                 help="How many threads did you wish to use " +
                                      "per parallel output default is 5.")
    albacore_parser.add_argument("--fastq_dir", type=str, required=False, default=None,
                                 help="Where should the fastq data be placed." +
                                      "Will be called 'fastq'" +
                                      "and sit adjacent to reads folder if left blank.")
    albacore_parser.add_argument("--hpc_dir", type=str, required=False, default=None,
                                 help="Where would you like to place the hpc batch files?")
    albacore_parser.add_argument("--max_processes", type=int, required=False, default=10,
                                 help="Limit the number of jobs that can be processed at any given moment." +
                                      "This command will prevent extraction/albacore jobs from being submitted while" +
                                      "there exists 'max_processes' jobs running / in the queue.")
    albacore_parser.add_argument("--log_dir", type=str, required=False, default=None,
                                 help="Where do you wish to store the poreduck logs? " +
                                      "Will be called 'poreduck_logs' and sit adjacent to reads folder if left blank")
    albacore_parser.add_argument("--albacore_version", type=str, required=True,
                                 help="Albacore 2 comes with the pass and fail folders." +
                                      "Expressed as <majore_version>.<minor_version>.<patch>")
    albacore_parser.add_argument("--qsub_albacore_template", type=str, required=True,
                                 help="The qsub albacore template for your qsub command. " +
                                      "Check the poreduck examples for more information.")
    albacore_parser.set_defaults(func=run_function)

    # Compare arguments
    compare_parser = subparsers.add_parser('compare_runs',
                                           help="This command takes in a samplesheet and plots each of the samples")
    compare_parser.add_argument("--samplesheet", type=str, required=True,
                                help="Similar to the samplesheet found in tarMyfast5 or albacoreHPC" +
                                     "But with an additional column 'Path'")

    compare_parser.set_defaults(func=run_function)
    args = parser.parse_args()

    # Print help if just 'poreduck' is typed in.
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
