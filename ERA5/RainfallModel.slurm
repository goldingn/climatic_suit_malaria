#!/bin/bash

#SBATCH --job-name=punim0361TempModel
#SBATCH --time=15-00 # days-hours format
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=51190

#SBATCH --array=1-25

#SBATCH --mail-type=ALL
#SBATCH --mail-user=owen.brown@unimelb.edu.au

module load R

Rscript -e "source('RainfallModel_Run.R');project(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')))"
