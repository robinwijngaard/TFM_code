#!/bin/bash


## Define routes

hg38="/home/robin/Documents/Project/Samples/hg38"
fastq="/home/robin/Documents/Project/Samples/example/fastq_ex"
bam="/home/robin/Documents/Project/Samples/example/bam_ex"

# Obtain filelist.txt: 

ls "$fastq" | cat | sed 's/_.*//g' | sort | uniq > "$fastq"/fastqlist.txt

#Obtain BAM



for file in `cat "$fastq"/fastqlist.txt`
do
	R1=`echo $file | sed 's/$/_R1.fastq.gz/g'`
	R2=`echo $file | sed 's/$/_R2.fastq.gz/g'`
	bwa mem -P -t 8 "$hg38"/hg38.fa "$fastq"/"$R1" "$fastq"/"$R2" | samtools sort -@ 8 -o "$bam"/"$file".bam -
	samtools index "$bam"/"$file".bam
done

