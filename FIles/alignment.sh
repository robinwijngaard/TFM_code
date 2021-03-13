#!/bin/bash

# Obtain filelist.txt: ls "/home/robin/Documents/Project/Samples/example/fastq_ex/" | cat | sed 's/_.*//g' | sort | uniq > filelist.txt

for file in `cat filelist.txt`
do
	R1=`echo $file | sed 's/$/_R1.fastq.gz/g'`
	R2=`echo $file | sed 's/$/_R2.fastq.gz/g'`
	bwa mem -P -t 8 /home/robin/Documents/Project/Samples/hg38/hg38.fa /home/robin/Documents/Project/Samples/example/fastq_ex/"$R1" /home/robin/Documents/Project/Samples/example/fastq_ex/example/fastq_ex/"$R2" | samtools sort -@ 8 -o /home/robin/Documents/Project/Samples/example/bam_ex/"$file".bam -
	samtools index /home/robin/Documents/Project/Samples/example/bam_ex/"$file".bam
done

