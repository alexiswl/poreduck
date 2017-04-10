#!/usr/bin/env python

import pandas as pd
import os
import time

dir = 'Z:\\bioinfo-proj-alexis/PoreCampAU_data/PoreCampAU\data/fast5/e_coli_R9'

fast5_files = [fast5_file for fast5_file in os.listdir(dir)
               if fast5_file.endswith(".fast5")]
os.chdir(dir)
fast5_pd = pd.DataFrame(columns=['filename', 'ctime', 'rnumber', 'mux', 'channel', 'read_no'])
fast5_pd['filename'] = fast5_files
fast5_pd['ctime'] = [time.ctime(os.path.getmtime(fast5_file)) for fast5_file in fast5_files]
fast5_pd['rnumber'] = [fast5_file.split('_')[-4] for fast5_file in fast5_files]
fast5_pd['mux'] = ["YES" if "mux_scan" in fast5_file else "NO" for fast5_file in fast5_files]
fast5_pd['channel'] = [fast5_file.split('_')[-3] for fast5_file in fast5_files]
fast5_pd['read_no'] = [fast5_file.split('_')[-2] for fast5_file in fast5_files]

runs = fast5_pd['rnumber'].unique().tolist()
print(runs)
for run in runs:
    print(fast5_pd.loc[(fast5_pd.rnumber == run) & (fast5_pd.mux == "NO")]['filename'])

fast5_pd.to_csv('C:\\Users\\lucatta\\PycharmProjects\\poreduck\\fast5.csv', header=True, index=False)

sub_int = "0_7276"

try:
    print(int(sub_int))
except ValueError:
    print(sub_int)

number = 0

print("%05d" % number)