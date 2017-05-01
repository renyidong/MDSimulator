#!/bin/sh
#SBATCH --time=00:59:59
#SBATCH --overcommit
#SBATCH --ntasks-per-node 64
#SBATCH --workdir=/gpfs/u/home/PCP6/PCP6rndn/scratch/
#

export BGLOCKLESSMPIO_F_TYPE=0x47504653
srun /gpfs/u/home/PCP6/PCP6rndn/barn/assignment4/a.out "$@"
