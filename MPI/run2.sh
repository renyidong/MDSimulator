#!/bin/bash
set -ex

function get_partition_size {
	if [ "$1" -le 64 ]; then
		PARTITION=small
	elif [ "$1" -le 512 ]; then
		PARTITION=medium
	elif [ "$1" -le 2048 ]; then
		PARTITION=large
	elif [ "$1" -le 4096 ]; then
		PARTITION=verylarge
	else
		echo "No such partition"
		exit
	fi
}

for nodes in 16 32 64 ; do
	get_partition_size $nodes
	OUTFILE="/gpfs/u/home/PCP6/PCP6rndn/scratch/proj2_${nodes}.out"
	sbatch --partition "$PARTITION" \
	       --nodes "$nodes" \
	       --output "$OUTFILE" /gpfs/u/home/PCP6/PCP6rndn/barn/project/sbatch2.sh 
done
