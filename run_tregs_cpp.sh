#!/bin/bash
#SBATCH --job-name=treg_array_cpp
#SBATCH --partition=epyc2
#SBATCH --array=0-999
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=logs/treg_cpp_%A_%a.out
#SBATCH --error=logs/treg_cpp_%A_%a.err

module load Anaconda3
source activate r-env

N_CHUNKS=1000

# SLURM_ARRAY_TASK_ID gives 0..999
CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

echo "Running chunk $CHUNK_ID on node $(hostname)"

Rscript /storage/homefs/bt25p365/tregs/UBX_datagen_cpp.R \
    $N_CHUNKS \
    $CHUNK_ID

echo "Chunk $CHUNK_ID completed"
