#!/bin/bash

bamall="/home/robin/Documents/Project/Samples/bam/all"
bamsingle="/home/robin/Documents/Project/Samples/bam/single"
bamlist="/home/robin/Documents/Project/TFM_code/Files/Other/singlebams.txt"

for bam in `cat "$bamlist"`
do
	bamfile=`echo $bam | sed 's/$/.bam/g'`
	baifile=`echo $bam | sed 's/$/.bam.bai/g'`

	cp "$bamall"/"$bamfile" "$bamsingle"
	cp "$bamall"/"$baifile" "$bamsingle"
done
	
