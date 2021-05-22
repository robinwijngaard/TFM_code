#!/bin/bash

#paths

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
bed="/home/robin/Documents/Project/Samples/bedfiles"
output="/home/robin/Documents/Project/Results/DeCON"


##/home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

cd /home/robin/Documents/Project/TFM_code/DECoN-master/Linux/

Rscript makeCNVcalls.R --Rdata "$output"/test.RData --custom FALSE --out "$output"/test --plot All --plotFolder "$output"/TestPlots