#!/bin/bash
# We have a list of variables in capital letters that will be replaced using sed.
# When we run albacore. One needs to manipulate the source files to suit their particular server.
# This is the parameter list.
#SBATCH --mem
#SBATCH --parsable
#SBATCH --time=200
#SBATCH --output
#SBATCH --error
#SBATCH --workdir
#SBATCH --job-name
#SBATCH --exclude=balder-wn06
# Other variables for analysis
# ALBACORE_DIR
# ALBACORE_VER
# FASTQ_DIR
# NUM_THREADS
# KIT
# FLOWCELL
# SUBFOLDER_NAME

# Create tmp directories to untar and save the folders
TMP_EXT=`mktemp -d /tmp/fast5.XXXXXXX`
TMP_SAVE=`mktemp -d /tmp/albacore.XXXXXXX`

# Return error if tar file is corrupted
tar_cmd="tar xf ${SUBFOLDER_NAME}.fast5.tar.gz -C ${TMP_EXT}"
eval $tar_cmd
ret_code=$?
if [ ${ret_code} != 0 ]; then
      printf "Error exit code [%d] when extracting tar file: ''${tar_cmd}'" ${ret_code}
      printf "Moving subfolder to .corrupted"
      touch ${SUBFOLDER_NAME}.fast5.tar.gz.corrupted
      exit ${ret_code}
fi

# Insert your source files here An example would be:
export OMP_NUM_THREADS=1

# Run albacore through CONDA
unset PYTHONPATH
CONDA_ROOTDIR=~/anaconda2

source ${CONDA_ROOTDIR}/bin/activate albacore_${ALBACORE_VER}
read_fast5_basecaller.py --input ${TMP_EXT}/${SUBFOLDER_NAME} \
                         --worker_threads ${NUM_THREADS} \
                         --save_path ${TMP_SAVE} \
                         --flowcell ${FLOWCELL} \
                         --kit ${KIT} \
			 --output_format fastq
source deactivate

# Move albacore log file to log directory.
mv ${TMP_SAVE}/pipeline.log ${ALBACORE_DIR}/${SUBFOLDER_NAME}.pipeline.log
mv ${TMP_SAVE}/sequencing_summary.txt ${ALBACORE_DIR}/${SUBFOLDER_NAME}.sequencing_summary.txt

# Move pass fastq file to fastq directory
mkdir -p ${FASTQ_DIR}/pass
cat ${TMP_SAVE}/workspace/pass/*.fastq | gzip > ${FASTQ_DIR}/pass/${SUBFOLDER_NAME}.fastq.gz

# Move failed fastq file to failed directory
mkdir -p ${FASTQ_DIR}/fail
cat ${TMP_SAVE}/workspace/fail/*.fastq | gzip > ${FASTQ_DIR}/fail/${SUBFOLDER_NAME}.fastq.gz

# Delete tmp directory
rm -rf ${TMP_SAVE} ${TMP_EXT}
