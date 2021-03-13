#!/bin/bash

# Obtain filelist.txt: 

ls "/home/robin/Documents/Project/Samples/example/fastq_ex/" | cat | sed 's/_.*//g' | sort | uniq > filelist.txt

#Obtain BAM

## Define routes

hg38="/home/robin/Documents/Project/Samples/hg38new"
fastq="/home/robin/Documents/Project/Samples/example/fastq_ex"
bam="/home/robin/Documents/Project/Samples/example/bam_ex"

for file in `cat filelist.txt`
do
	R1=`echo $file | sed 's/$/_R1.fastq.gz/g'`
	R2=`echo $file | sed 's/$/_R2.fastq.gz/g'`
	bwa mem -P -t 8 "$hg38"/hg38.fa "$fastq"/"$R1" "$fastq"/"$R2" | samtools sort -@ 8 -o "$bam"/"$file".bam -
	samtools index "$bam"/"$file".bam
done

