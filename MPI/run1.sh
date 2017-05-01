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

for nodes in 1 2 4 8 16 ; do
	get_partition_size $nodes
	OUTFILE="/gpfs/u/home/PCP6/PCP6lnrh/scratch/proj1_${nodes}.out"
	sbatch --partition "$PARTITION" \
	       --nodes "$nodes" \
	       --output "$OUTFILE" /gpfs/u/home/PCP6/PCP6lnrh/scratch/sbatch1.sh 
done

