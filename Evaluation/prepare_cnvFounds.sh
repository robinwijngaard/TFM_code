#!/bin/bash

## For in-house dataset
cnvFounds=/home/robin/Documents/Project/Files/Clinic/cnvFounds  # dir where results are exported
outputbenchmarker=/home/robin/Documents/Project/Files/Clinic/outputbechmarker # dir of the output of the CNVbechmarkeR

Rscript prepare_cnvFounds.R --cnvFounds $cnvFounds --outputbenchmarker $outputbenchmarker


