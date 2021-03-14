#!/bin/bash

#paths

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
bed="/home/robin/Documents/Project/Samples/bedfiles"
ref="/home/robin/Documents/Project/Samples/hg38new"
output="//home/robin/Documents/Project/Results/DeCON"



find "$bams" -name *.bam > "$bams"/bamfiles.txt


##/home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

cd /home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

/home/robin/Downloads/R-3.1.2/bin/Rscript IdentifyFailures.R --Rdata "$output"/DECoNtest.RData --mincorr .98 --mincov 100 --custom FALSE --out QCtest
