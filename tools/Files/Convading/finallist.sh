#!/bin/bash

rundir="/home/robin/Documents/Project/TFM_code/CoNVaDING-master"

counts="/home/robin/Documents/Project/Results/Convading/CreateCounts"
qcsamples="/home/robin/Documents/Project/Results/Convading/qcsamples"
cnv="/home/robin/Documents/Project/Results/Convading/cnvdetection"
targetqc="/home/robin/Documents/Project/Results/Convading/targetqc"
finallist="/home/robin/Documents/Project/Results/Convading/finallist"

cd "$rundir"

perl ./CoNVaDING.pl \
-mode CreateFinalList \
-inputDir "$cnv" \
-outputDir "$finallist" \
-targetQcList "$targetqc"