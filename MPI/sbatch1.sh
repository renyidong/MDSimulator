#!/bin/sh
#SBATCH --time=00:59:59
#SBATCH --overcommit
#SBATCH --workdir=/gpfs/u/home/PCP6/PCP6lnrh/scratch/
#

srun --ntasks-per-node 16 /gpfs/u/home/PCP6/PCP6lnrh/scratch/project/a.out 100 16 1
srun --ntasks-per-node 32 /gpfs/u/home/PCP6/PCP6lnrh/scratch/project/a.out 100 16 1
srun --ntasks-per-node 64 /gpfs/u/home/PCP6/PCP6lnrh/scratch/project/a.out 100 16 1

