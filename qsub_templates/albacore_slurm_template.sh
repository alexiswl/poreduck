#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#SBATCH -J albacore
#SBATCH --parsable
#SBATCH -o %STDOUT%
#SBATCH --time=200
#SBATCH -e %STDERR%
#SBATCH --mem=%MEM%g
#SBATCH -D %WORKING_DIRECTORY%
# Insert your source files here around our 'command'. An example would be:
export OMP_NUM_THREADS=1
source activate albacore_%ALBACORE_VER%
%COMMAND%
source deactivate
