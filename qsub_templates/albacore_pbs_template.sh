#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#PBS -N albacore
#PBS -o %STDOUT%
#PBS -e %STDERR%
#PBS -S /bin/bash
#PBS -l h_vmem=%MEM%g
#PBS -wd %WORKING_DIRECTORY%
#PBS -v OMP_NUM_THREADS=1
# Insert your source files here around our 'command'. An example would be:
source activate albacore_%VER%
%COMMAND%
source deactivate