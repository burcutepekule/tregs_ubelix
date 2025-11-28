#!/bin/bash
#SBATCH --job-name=treg0_array
#SBATCH --partition=epyc2
#SBATCH --array=0-499
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=logs/treg0_%A_%a.out
#SBATCH --error=logs/treg0_%A_%a.err

module load Anaconda3
source activate r-env

# Fixed parameters for this run
STERILE=1
ALLOW_TREGS=0
RANDOMIZE_TREGS=0
N_CHUNKS=500

# SLURM_ARRAY_TASK_ID gives 0..499
CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

echo "Running chunk $CHUNK_ID on node $(hostname)"

Rscript /storage/homefs/bt25p365/tregs/UBX_datagen.R \
    $STERILE \
    $ALLOW_TREGS \
    $RANDOMIZE_TREGS \
    $N_CHUNKS \
    $CHUNK_ID

echo "Chunk $CHUNK_ID completed"
