#!/bin/bash
#SBATCH --job-name=run_c_theta
#SBATCH --output=run_c_theta.out
#SBATCH --error=../logs/run_c_theta.err
#SBATCH --time=03:00:00             # Adjust based on expected runtime
#SBATCH --mem=4G                    # Adjust memory based on needs
#SBATCH --cpus-per-task=1          # Change if your script is multi-threaded


module load R/4.3.2b
module load gcc/9.2.0

Rscript different_c_theta_v2.R