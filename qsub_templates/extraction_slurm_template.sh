#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#SBATCH -J extraction
#SBATCH --parsable
#SBATCH -o %STDOUT%
#SBATCH -e %STDERR%
#SBATCH -D %WORKING_DIRECTORY%
# Insert your source files here around our 'command'.
%COMMAND%
