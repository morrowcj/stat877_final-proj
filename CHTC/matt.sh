#!/bin/bash
#SBATCH -o matt.out
#SBATCH -e matt_error.out
#SBATCH -J matt
#SBATCH -D /workspace/ting
#SBATCH -p long
#SBATCH --mem 2000
module load R/R-3.5.3
R CMD BATCH --no-save matt.R matt.out