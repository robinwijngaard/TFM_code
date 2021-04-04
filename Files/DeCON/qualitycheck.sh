#!/bin/bash

#paths

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
bed="/home/robin/Documents/Project/Samples/bedfiles"
output="/home/robin/Documents/Project/Results/DeCON"


##/home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

cd /home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

Rscript IdentifyFailures.R --Rdata "$output"/test.RData --mincorr .98 --mincov 100 --custom FALSE --out "$output"/QCtest
