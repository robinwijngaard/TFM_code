#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  #outputDir
bedFile=/home/robin/Documents/Project/Files/ICR96/ICR96_hg38_noSNP.bed #bedFile 
roisFile=/home/robin/Documents/Project/Files/ICR96/validated_rois.bed #validated ROI file (output obtain_validated_ROIS)
cnvFounds=/home/robin/Documents/Project/Files/ICR96/cnvFounds #cnvFounds file (output CNVbenchmarkeR, one file per algorithm)

Rscript evaluation.r --outputDir $outputDir --bedFile $bedFile --roisFile $roisFile --cnvFounds $cnvFounds


## For in-house dataset
outputDir=/home/robin/Documents/Project/Files/Clinic  #outputDir
bedFile=/home/robin/Documents/Project/Files/Clinic/clinic.bed #bedFile 
roisFile=/home/robin/Documents/Project/Files/Clinic/validated_rois.bed #validated ROI file (output obtain_validated_ROIS)
cnvFounds=/home/robin/Documents/Project/Files/Clinic/cnvFounds #cnvFounds file (output CNVbenchmarkeR, one file per algorithm)

Rscript evaluation.r --outputDir $outputDir --bedFile $bedFile --roisFile $roisFile --cnvFounds $cnvFounds

