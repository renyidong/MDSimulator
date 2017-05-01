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

function wait_for_queue {
	while true; do
		[ $QLEN -gt 8 ] || return
		sleep 300
	done
}

QLEN=`squeue | wc -l`
for nfiles in 1 16 64 128 ; do
for block_size in 1 2 4 8 16 ; do
for nodes in 8 16 32 64 128 ; do
	while [ $QLEN -gt 8 ]; do
		sleep 300
		QLEN=`squeue | wc -l`
	done
	get_partition_size $nodes
	OUTFILE="/gpfs/u/home/PCP6/PCP6rndn/scratch/hw4_${nodes}_${nfiles}_${block_size}.out"
	[ -e $OUTFILE ] && continue
	touch $OUTFILE
	sbatch --partition "$PARTITION" \
	       --nodes "$nodes" \
	       --output "$OUTFILE" ./sbatch.sh "$nfiles" "$block_size" "tmp"
	QLEN=$((QLEN+1))
	[ $QLEN -le 8 ] && continue
done
done
done
