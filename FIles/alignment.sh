#!/bin/bash

# Obtain filelist.txt: ls "fastq/" | cat | sed 's/_.*//g' | sort | uniq > filelist.txt

for file in `cat filelist.txt`
do
	R1=`echo $file | sed 's/$/_R1.fastq.gz/g'`
	R2=`echo $file | sed 's/$/_R2.fastq.gz/g'`
	bwa mem -P -t 8 hg38ref/hg38index example/"$R1" example/"$R2" | samtools sort -@ 8 -o example/"$file".bam -
	samtools index example/"$file".bam
done

