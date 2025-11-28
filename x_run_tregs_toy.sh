#!/bin/bash
#SBATCH --job-name=treg_array
#SBATCH --partition=epyc2
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=logs/treg_%A_%a.out
#SBATCH --error=logs/treg_%A_%a.err

module load Anaconda3
source activate r-env

N_CHUNKS=10

# SLURM_ARRAY_TASK_ID gives 0..9
CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

echo "Running chunk $CHUNK_ID on node $(hostname)"

Rscript /storage/homefs/bt25p365/tregs/UBX_datagen_toy.R \
    $N_CHUNKS \
    $CHUNK_ID

echo "Chunk $CHUNK_ID completed"
