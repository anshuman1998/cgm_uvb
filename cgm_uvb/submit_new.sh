#!/bin/bash -login
# Propogate environment variables to compute node
#SBATCH --export=ALL

# set the number of nodes and processes per node
#SBATCH --nodes=6

# set the number of tasks (processes) per node.
#SBATCH --ntasks-per-node=28

# set name of job
#SBATCH --job-name=vikram_NHI14_18

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=vikramkhaire@iist.ac.in


echo $SLURM_JOB_ID
echo $SLURM_NPROCS


export PATH="/home/vikram/anaconda3/bin:$PATH"

python vikram_cloudy_run_metal_grid_14.py

