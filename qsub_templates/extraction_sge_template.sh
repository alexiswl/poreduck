#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#$ -N extraction
#$ -o %STDOUT%
#$ -e %STDERR%
#$ -S /bin/bash
#$ -wd %WORKING_DIRECTORY%
# Insert your source files here around our 'command'.
%COMMAND%
