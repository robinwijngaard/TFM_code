#!/bin/bash

#paths

bams="/home/robin/Documents/Project/Samples/example/bam_ex"                                                                   #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
bed="/home/robin/Documents/Project/Samples/bedfiles"                                                                       #name of bed file
fasta="/home/robin/Documents/Project/Samples/hg38"                                                                     #name of fasta file
exonfile="/home/robin/Documents/Project/Samples/bedfiles"
output="/home/robin/Documents/Project/Results/Panelcnmops"

##find "$bams" -name *.bam > "$bams"/bamfiles.txt



##/home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

cd /home/robin/Documents/Project/TFM_code/Panelcnmops

Rscript runPanelcnmops.R --bams "$bams" --bed "$bed"/ICR96_hg38.bed --fasta "$ref"/hg38.fa --exon "$exonfile"/exons.hg38.bed --out "$output"


#/home/robin/Downloads/R-3.1.2/bin/Rscript ReadInBams.R --bams /home/robin/Documents/Project/TFM_code/FIles/bamfile.txt --bed /home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38.bed --fasta /home/robin/Documents/Project/Samples/hg38ncbi/hg38.fna --out DEConTest

