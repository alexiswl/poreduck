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
    if args.command == "transfer_to_server":
        import poreduck.transfer_fast5_to_server as command_to_run
    if args.command == "albacore_on_server":
        import poreduck.albacore_server_scaled as command_to_run
    if args.command == "qc_plots":
        import poreduck.plot_yields as command_to_run
    if args.command == "compare_plots":
        import poreduck.compare_runs as command_to_run
    # Now run it!
    command_to_run.main(args)


# Globals for albacore server
CONFIGS = {"FC106_RAD001": "r94_250bps_linear.cfg",  # Rapid sequencing
           "FC106_LSK208_2d": "r94_250bps_2d.cfg",   # 2D unsure which one
           "FC106_LSK108": "r94_450bps_linear.cfg",  # For 1D ligation sequencing.
           "FC106_RAD002": "r94_450bps_linear.cfg"}  # Second Rapid sequencing kit.
FLOWCELLS = ["FLO-MIN107", "FLO-MIN106"]
KITS = ["SQK-LWP001", "SQK-NSK007", "VSK-VBK001", "SQK-RAS201", "SQK-RBK001", "SQK-LWB001",
        "SQK-RNA001", "SQK-RLI001", "SQK-RAD002", "SQK-RLB001", "SQK-RAB201", "SQK-LSK208",
        "SQK-LSK108", "SQK-RAD003", "SQK-DCS108", "SQK-PCS108", "SQK-LSK308"]  # More kits to come
QSUB_TYPES = ["SGE", "TORQUE", "SLURM"]

# Define main script:


def main():
    # Create poreduck parser
    parser = argparse.ArgumentParser(prog='poreduck', description="poreduck package")
    parser.add_argument("--version", help="Get version of poreduck",
                        action="version",
                        version=poreduck.__version__)

    subparsers = parser.add_subparsers(help="Callable poreduck functions", dest="command")

    # Transfer to server arguments
    transfer_parser = subparsers.add_parser('transfer_to_server',
                                            help="The transfer_fast5_to_server transfers MiNION data from a " +
                                                 "laptop in realtime. The process will finish when the script " +
                                                 "believes MinKNOW is no longer running.")

    transfer_parser.add_argument("--reads_dir", type=str, required=True,
                                 help="/path/to/reads, " +
                                      "should have a bunch of runs labelled <YYYYMMDD_HHMM_SAMPLE_NAME>")
    transfer_parser.add_argument("--server_name", type=str, required=True,
                                 help="If you were to ssh username@server, " +
                                      "please type in the server bit.")
    transfer_parser.add_argument("--user_name", type=str, required=True,
                                 help="If you were to ssh username@server, " +
                                      "please type in the username bit")
    transfer_parser.add_argument("--dest_directory", type=str, required=True,
                                 help="Where abouts on the server do you wish to place these files?")
    transfer_parser.add_argument("--sample_name", type=str, required=True,
                                 help="Sample name that you typed into MinKNOW.")
    transfer_parser.add_argument("--suffix", type=str, required=False, default=None,
                                 help="Would you like a suffix at the end of each of your csv and tar files?")
    transfer_parser.add_argument("--sshpass", default=False, dest='sshpass', action='store_true',
                                 help="ssh-pass is rather poor practise and quite a security risk." +
                                      "Tick this option and set up an id_rsa key if you'd prefer." +
                                      "You will still be required to enter your password for set-up purposes.")
    transfer_parser.add_argument("--local", default=False, dest='local', action='store_true',
                                 help="Has local basecalling been used, changes directory from fast5 to fast5/pass")
    transfer_parser.set_defaults(func=run_function)

    # Albacore server arguments:
    albacore_parser = subparsers.add_parser('albacore_on_server',
                                            help="The albacore-server scaled command incorporates qsub to spread " +
                                                 "the server load and rapidly " + "generate data off the MinION.")
    albacore_parser.add_argument("--reads_dir", type=str, required=True,
                                 help="/path/to/reads, " +
                                      "should have a bunch of tar zipped files in it.")
    albacore_parser.add_argument("--flowcell", type=str, choices=FLOWCELLS,
                                 help="Pick a flowcell version")
    albacore_parser.add_argument("--kit", type=str, choices=KITS,
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
    albacore_parser.add_argument("--resume", type=str, required=False, default=None,
                                 help="Resume the albacore run, need a csv file")
    albacore_parser.add_argument("--qsub_directory", type=str, required=False, default=None,
                                 help="Where would you like to place the qsub files?")
    albacore_parser.add_argument("--qsub_host", type=str, required=False, default=None,
                                 help="Where would you like the qsub jobs to be run?")
    albacore_parser.add_argument("--max_processes", type=int, required=False, default=None,
                                 help="Limit the number of jobs that can be processed at any given moment." +
                                      "This command will prevent extraction/albacore jobs from being submitted while" +
                                      "there exists 'max_processes' jobs running / in the queue.")
    albacore_parser.add_argument("--log_directory", type=str, required=False, default=None,
                                 help="Where do you wish to store the poreduck logs? " +
                                      "Will be called 'poreduck_logs' and sit adjacent to reads folder if left blank")
    albacore_parser.add_argument("--barcoding", default=False, dest='barcoding', action='store_true',
                                 help="Use this option to demultiplex library?")
    albacore_parser.add_argument("--qsub_type", choices=QSUB_TYPES, default="SGE",
                                 help="What qsub system are you using?")
    albacore_parser.add_argument("--qsub_extraction_template", type=str, required=True,
                                 help="The qsub extraction template for your qsub command. " +
                                      "Check the poreduck examples for more information.")
    albacore_parser.add_argument("--qsub_albacore_template", type=str, required=True,
                                 help="The qsub albacore template for your qsub command. " +
                                      "Check the poreduck examples for more information.")
    albacore_parser.set_defaults(func=run_function)

    # Plotting arguments
    plot_parser = subparsers.add_parser('qc_plots',
                                        help="This plot takes in the csv files along with the fastq files to produce some " +
                                             "qc plots of the data")
    csv_args = plot_parser.add_mutually_exclusive_group(required=True)
    csv_args.add_argument("--csv_dir", type=str,
                             help="/path/to/csv_dir" +
                                  "Should have a bunch of csv files in it.")
    csv_args.add_argument("--no_csv", default=False, action='store_true', dest="no_csv",
                          help="I didn't use poreduck to upload my data to my server." +
                               "I promise to do so next time to reap all of the benefits.")
    plot_parser.add_argument("--fastq_dir", type=str, required=True,
                             help="/path/to/fastq/files")
    plot_parser.add_argument("--output_dir", type=str, required=True,
                             help="/path/to/plots_dir" +
                                  "By default will be created in current working directory")
    plot_parser.add_argument("--clip", default=False, action='store_true', dest="clip",
                             help="Remove outliers from the historgram ( >3 std.dev)")
    plot_parser.add_argument("--sample_name", type=str, required=False,
                             help="Name to add onto each of the plots")
    plot_parser.set_defaults(func=run_function)

    # Compare arguments
    compare_parser = subparsers.add_parser('compare_plots',
                                           help="This command takes in two fastq directories, and creates comparison plots")
    compare_parser.add_argument("--fastq_dirs", type=str, required=True,
                                help="/path/to/fastq_1,/path/to/fastq_2, etc")
    compare_parser.add_argument("--run_names", type=str, required=True,
                                help="\"name of run 1\",\"name of run 2\", etc")
    compare_parser.add_argument("--plots_dir", type=str, required=True,
                                help="where do you wish these plots to go?")
    compare_parser.add_argument("--clip", default=False, action='store_true', dest="clip",
                                help="Remove outliers from the historgram (0.999th percentile)")

    compare_parser.set_defaults(func=run_function)
    args = parser.parse_args()

    # Print help if just 'poreduck' is typed in.
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
