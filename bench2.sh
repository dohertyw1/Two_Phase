#!/bin/bash --login
###
#SBATCH --job-name=MPLS1
#SBATCH -A scw1239
#SBATCH --output=bench2.out.%J
#SBATCH --error=bench2.err.%J
#SBATCH --mem-per-cpu=4000
#SBATCH -n 60
#SBATCH --ntasks-per-node=40
#SBATCH -t 20:00:00
###

module purge

module load anaconda/2020.02

eval "$(/apps/languages/anaconda/2020.02/bin/conda shell.bash hook)"

module list

conda activate fenicsproject

which mpirun

which python

python --version

export OMP_NUM_THREADS=1

NTASKS=$SLURM_NTASKS

mpirun -n $NTASKS python3 droplet_ve_dim_con.py
