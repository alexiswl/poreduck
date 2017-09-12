#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#$ -N albacore
#$ -terse
#$ -o %STDOUT%
#$ -e %STDERR%
#$ -S /bin/bash
#$ -l h_vmem=%MEM%g
#$ -wd %WORKING_DIRECTORY%
#$ -v OMP_NUM_THREADS=1
# Insert your source files here around our 'command'. An example would be:
source activate albacore_env
%COMMAND%
source deactivate