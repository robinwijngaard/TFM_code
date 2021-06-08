#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  #outputDir
failedRois=/home/robin/Documents/Project/Files/ICR96/cnvFounds #failed ROIs file (output of CNVbenchmarkeR)
roisFile=/home/robin/Documents/Project/Files/ICR96/validated_rois.bed #validated ROIs file

Rscript failedrois.r --outputDir $outputDir --failedRois $failedRois --roisFile $roisFile