#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

"""
This script will create a yield plot of the data that has been created by the
albacore_server_scaled.py script and placed in the folder transfer_fast5_to_server.py
"""

# Check python version


# Arguments required
"""
csv directory
fastq directory
output directory
"""

# Get csv's
# Required to generate the yield plot.
# Each csv is written to a pandas dataframe.

# Use SeqIO to map reads to length
# Add this to end of dataframe.

# Concatenate csvs

# Sort by date column.

# Group by minute, aggregate read-lengths

# Aggregate yield by minute

# Plot yield using matploy lib
