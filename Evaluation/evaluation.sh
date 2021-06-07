#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  #outputDir
bedFile=/home/robin/Documents/Project/Files/ICR96/ICR96_hg38_noSNP.bed #bedFile 
roisFile=/home/robin/Documents/Project/Files/ICR96/validated_rois.bed
cnvFounds=/home/robin/Documents/Project/Files/ICR96/cnvFounds

Rscript evaluation.r --outputDir $outputDir --bedFile $bedFile --roisFile $roisFile --cnvFounds $cnvFounds


## For clinic dataset
outputDir=/home/robin/Documents/Project/Files/Clinic  #outputDir
bedFile=/home/robin/Documents/Project/Files/Clinic/clinic.bed #bedFile 
roisFile=/home/robin/Documents/Project/Files/Clinic/validated_rois.bed
cnvFounds=/home/robin/Documents/Project/Files/Clinic/cnvFounds

Rscript evaluation.r --outputDir $outputDir --bedFile $bedFile --roisFile $roisFile --cnvFounds $cnvFounds

