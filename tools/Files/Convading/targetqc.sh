#!/bin/bash

rundir="/home/robin/Documents/Project/TFM_code/CoNVaDING-master"

counts="/home/robin/Documents/Project/Results/Convading/CreateCounts"
qcsamples="/home/robin/Documents/Project/Results/Convading/qcsamples"
cnv="/home/robin/Documents/Project/Results/Convading/cnvdetection"
targetqc="/home/robin/Documents/Project/Results/Convading/targetqc"
counts="/home/robin/Documents/Project/Results/Convading/CreateCounts"

cd "$rundir"

perl ./CoNVaDING.pl \
-mode GenerateTargetQcList \
-outputDir "$targetqc" \
-controlsDir "$counts" \
-inputDir "$counts"