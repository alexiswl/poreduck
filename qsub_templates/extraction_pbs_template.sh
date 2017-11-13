#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#PBS -N extraction
#PBS -o %STDOUT%
#PBS -e %STDERR%
#PBS -S /bin/bash
#PBS -wd %WORKING_DIRECTORY%
# Insert your source files here around our 'command'.
%COMMAND%