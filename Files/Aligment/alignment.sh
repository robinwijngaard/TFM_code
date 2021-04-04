#!/bin/bash


## Define routes

hg38="/home/robin/Documents/Project/Samples/hg38"
fastq="/home/robin/Documents/Project/Samples/fastq"
bam="/home/robin/Documents/Project/Samples/bam"
alignmdir="/home/robin/Documents/Project/TFM_code/Files/Aligment"


# Obtain filelist.txt: 

ls "$fastq" | cat | sed 's/_.*//g' | sort | uniq > "$alignmdir"/fastqlist.txt

#Obtain BAM



for file in `cat "$alignmdir"/fastqlist.txt`
do
	R1=`echo $file | sed 's/$/_R1.fastq.gz/g'`
	R2=`echo $file | sed 's/$/_R2.fastq.gz/g'`
	bwa mem -t 8 "$hg38"/hg38.fa "$fastq"/"$R1" "$fastq"/"$R2" | samtools sort -@ 8 -o "$bam"/"$file".bam -
	samtools index "$bam"/"$file".bam
done

