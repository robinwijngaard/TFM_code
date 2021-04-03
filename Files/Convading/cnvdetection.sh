#!/bin/bash

rundir="/home/robin/Documents/Project/TFM_code/CoNVaDING-master"

counts="/home/robin/Documents/Project/Results/Convading/CreateCounts"
qcsamples="/home/robin/Documents/Project/Results/Convading/qcsamples"
cnv="/home/robin/Documents/Project/Results/Convading/cnvdetection"

cd "$rundir"

perl ./CoNVaDING.pl \
-mode StartWithBestScore \
-inputDir "$qcsamples" \
-outputDir "$cnv" \