#!/bin/bash
# Number of MPI tasks
#SBATCH --job-name=polcal-RA
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=480GB
#SBATCH --output=polcal-RA_%A_%a.log
#SBATCH --nodelist=c005
#SBATCH --array=1-60


python polcal.G182.RA.py $SLURM_ARRAY_TASK_ID
