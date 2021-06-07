#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  #outputDir
failedRois=/home/robin/Documents/Project/Files/ICR96/cnvFounds
roisFile=/home/robin/Documents/Project/Files/ICR96/validated_rois.bed

Rscript failedrois.r --outputDir $outputDir --failedRois $failedRois --roisFile $roisFile