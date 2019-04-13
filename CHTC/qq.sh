#!/bin/bash
#SBATCH -o qq.out
#SBATCH -e qq_error.out
#SBATCH -J qq
#SBATCH -D /workspace/ting
#SBATCH -p short
#SBATCH --mem 4000
module load R/R-3.5.3
R CMD BATCH --no-save qq.R qq.out