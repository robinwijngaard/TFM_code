#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  # outputDir
bedFile=/home/robin/Documents/Project/Files/ICR96/ICR96_hg38_noSNP.bed # bedFile 
validatedFile=/home/robin/Documents/Project/Files/ICR96/ICR96_validated_regions38.bed # Validated file ICR96 format

Rscript obtain_validated_ROIS.R --outputDir $outputDir --bedFile $bedFile --validatedFile $validatedFile --EPCAM TRUE


## For in-house dataset
outputDir=/home/robin/Documents/Project/Files/Clinic  # outputDir
bedFile=/home/robin/Documents/Project/Files/Clinic/clinic.bed # bedFile 
validatedFile=/home/robin/Documents/Project/Files/Clinic/validated_ICR96format.bed # Validated file ICR96 format

Rscript obtain_validated_ROIS.R --outputDir $outputDir --bedFile $bedFile --validatedFile $validatedFile --EPCAM FALSE

