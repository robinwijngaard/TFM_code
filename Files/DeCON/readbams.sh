#!/bin/bash

#paths

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
bed="/home/robin/Documents/Project/Samples/bedfiles"
ref="/home/robin/Documents/Project/Samples/hg38"
output="/home/robin/Documents/Project/Results/DeCON"


##find "$bams" -name *.bam > "$bams"/bamfiles.txt



##/home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

cd /home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

Rscript ReadInBams.R --bams "$bams" --bed "$bed"/ICR96_hg38_noSNP.bed --fasta "$ref"/hg38.fa --out "$output"/test


#/home/robin/Downloads/R-3.1.2/bin/Rscript ReadInBams.R --bams /home/robin/Documents/Project/TFM_code/FIles/bamfile.txt --bed /home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38.bed --fasta /home/robin/Documents/Project/Samples/hg38ncbi/hg38.fna --out DEConTest

