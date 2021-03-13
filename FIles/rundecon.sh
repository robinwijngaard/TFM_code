#!/bin/bash

find ~/Documents/Project/Samples/example/bam_ex -name *.bam > bamfiles.txt

#paths

decon="/home/robin/Documents/Project/TFM_code/DECoN-master/Linux"
bams="/home/robin/Documents/Project/TFM_code/FIles/"
bed="/home/robin/Documents/Project/Samples/bedfiles"
ref="/home/robin/Documents/Project/Samples/hg38ncbi"

/home/robin/Downloads/R-3.1.2/bin/Rscript "$decon"/ReadInBams.R --bams "$bams"/bamfiles.txt --bed "$bed"/ICR96_hg38.bed --fasta "$ref"/hg38.fna --out DECoNtest