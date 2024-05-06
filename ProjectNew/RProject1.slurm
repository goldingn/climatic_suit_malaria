#!/bin/bash

#SBATCH --job-name=RProject1New
#SBATCH -p cloud
#SBATCH --time=15-00 # days-hours format
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=51190

#SBATCH --array=1-25

#SBATCH --mail-type=ALL
#SBATCH --mail-user=owenb696@gmail.com

module load R

##Rscript -e "install.packages('tempsuitcalc_0.1.tar.gz', type = 'source', repos = NULL)"

Rscript -e "source('ProjectInputNew.R');project(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')))"

