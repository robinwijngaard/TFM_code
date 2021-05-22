#!/bin/bash

rundir="/home/robin/Documents/Project/TFM_code/CoNVaDING-master"

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
targets="/home/robin/Documents/Project/Samples/bedfiles"
counts="/home/robin/Documents/Project/Results/Convading/CreateCounts"
ref="/home/robin/Documents/Project/Samples/hg38"


cd "$rundir"

perl ./CoNVaDING.pl \
-mode StartWithBam \
-inputDir "$bams" \
-useSampleAsControl \
-controlsDir "$counts" \
-outputDir "$counts" \
-bed "$targets"/ICR96_hg38.bed