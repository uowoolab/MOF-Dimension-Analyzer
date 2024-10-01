#!/bin/bash
#SBATCH --job-name=SPLITTER      # Job name
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=1        # Number of CPUs cores per task
#SBATCH --mem=2G                 # Total memory allocation
#SBATCH --time=0-24:00:00        # Time limit for the job

# Modify this for your use case
INPUT_DIR=/path/to/directory/containing/cif/files
BATCH_SIZE=50000                                     # Number of files per batch
SEND_SCRIPT=/path/to/SEND.sh/script/SEND.sh

# Calculate total number of files in the input directory
TOTAL_FILES=$(find $INPUT_DIR -type f | wc -l)

# Determine how many batches are needed
NUM_BATCHES=$((TOTAL_FILES / BATCH_SIZE))       # Calculate the base number of batches (integer division)
if [ $((TOTAL_FILES % BATCH_SIZE)) -ne 0 ]; then
    NUM_BATCHES=$((NUM_BATCHES + 1))            # If there are leftover files, add one more batch
fi

# Loop through the number of batches to process files in chunks
for (( BATCH=1; BATCH<=NUM_BATCHES; BATCH++ ))
do
    # Calculate the range of file indices for the current batch
    START=$(( (BATCH - 1) * BATCH_SIZE + 1 ))
    END=$(( BATCH * BATCH_SIZE ))
    if [ $END -gt $TOTAL_FILES ]; then
        END=$TOTAL_FILES
    fi
    
    # Create a temporary directory for the current batch
    TEMP_BATCH_DIR="$INPUT_DIR/batch_$BATCH"
    mkdir -p $TEMP_BATCH_DIR
    
    # Link files to the temporary batch directory
    # Assuming files in INPUT_DIR are named in a sortable order
    CURRENT_FILE=0
    find $INPUT_DIR -type f | sort | while read FILENAME; do
        CURRENT_FILE=$((CURRENT_FILE + 1))
        if [ $CURRENT_FILE -ge $START ] && [ $CURRENT_FILE -le $END ]; then
            ln -s "$FILENAME" "$TEMP_BATCH_DIR/"    # Create a symlink in the batch directory for the current file
        fi
        if [ $CURRENT_FILE -gt $END ]; then
            break
        fi
    done

    # Submit the batch for processing
    sbatch $SEND_SCRIPT $TEMP_BATCH_DIR $BATCH
    
done

echo "All batches submitted."
