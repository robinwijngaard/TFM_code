#!/bin/bash

find ~/Documents/Project/Samples/example/bam_ex -name *.bam > bamfiles.txt

#paths

bams="/home/robin/Documents/Project/TFM_code/FIles"
bed="/home/robin/Documents/Project/Samples/bedfiles"
ref="/home/robin/Documents/Project/Samples/hg38new"


##/home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

cd /home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

/home/robin/Downloads/R-3.1.2/bin/Rscript ReadInBams.R --bams "$bams"/bamfiles.txt --bed "$bed"/ICR96_hg38.bed --fasta "$ref"/hg38.fa --out DECoNtest


#/home/robin/Downloads/R-3.1.2/bin/Rscript ReadInBams.R --bams /home/robin/Documents/Project/TFM_code/FIles/bamfile.txt --bed /home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38.bed --fasta /home/robin/Documents/Project/Samples/hg38ncbi/hg38.fna --out DEConTest

