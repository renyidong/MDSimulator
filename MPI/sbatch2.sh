#!/bin/sh
#SBATCH --time=00:59:59
#SBATCH --overcommit
#SBATCH --workdir=/gpfs/u/home/PCP6/PCP6rndn/scratch/
#

srun --ntasks-per-node 16 /gpfs/u/home/PCP6/PCP6rndn/barn/project/a.out 100 16 1
srun --ntasks-per-node 16 /gpfs/u/home/PCP6/PCP6rndn/barn/project/a.out 100 16 2
srun --ntasks-per-node 16 /gpfs/u/home/PCP6/PCP6rndn/barn/project/a.out 100 16 4
