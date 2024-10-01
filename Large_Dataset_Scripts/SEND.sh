#!/bin/bash
#SBATCH --job-name=SEND_PY              # Job name
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=32G                       # Total memory allocation
#SBATCH --time=7-00:00:00               # Time limit for the job

# Load dependencies (adjust as needed for your cluster)
module purge
module load csds

# Retrieve arguments passed from the SUBMIT.sh script
# Example: sbatch $SEND_SCRIPT $TEMP_BATCH_DIR $BATCH
INPUT_DIR=$1
BATCH_NUM=$2

PYTHON_FILE=/path/to/python/code/MOF_Dimensions_CIF_Parallel.py
OUTPUT_FILE=/path/to/output/directory/output_${BATCH_NUM}.csv
NUM_CORES=8

# Executes the Python script with specified input and output files, using 8 cpu cores
python $PYTHON_FILE -i $INPUT_DIR -o $OUTPUT_FILE -n $NUM_CORES
