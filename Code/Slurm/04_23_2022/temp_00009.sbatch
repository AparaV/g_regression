#!/bin/bash
#SBATCH --job-name prelims_V50_A1_n100 
#SBATCH --partition medium                         # Slurm partition to use
#SBATCH --ntasks 5                                 # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 6-23:59                             # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=4000                         # Memory limit for each tasks (in MB)
#SBATCH -o /homes/aparav/slurm_logs/%x_%j.out      # File to which STDOUT will be written
#SBATCH -e /homes/aparav/slurm_logs/%x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL                            # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=aparav@uw.edu                  # Email to which notifications will be sent

module load R/4.1.0-foss-2021a

cd ~/prelims/lsq_causal_effects/Code

Rscript run_simulations.R \
    --num_nodes 50 \
    --treatment_size 1 \
    --n 100 \
    --num_reps 1000 \
    --output ../Results/V50_A1_n100.rds 