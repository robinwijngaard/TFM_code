#!/bin/bash

rundir="/home/robin/Documents/Project/TFM_code/CoNVaDING-master"

counts="/home/robin/Documents/Project/Results/Convading/CreateCounts"
qcsamples="/home/robin/Documents/Project/Results/Convading/qcsamples"

cd "$rundir"

perl ./CoNVaDING.pl \
-mode StartWithMatchScore \
-inputDir "$counts" \
-outputDir "$qcsamples" \
-controlsDir "$counts" \
-controlSamples 6