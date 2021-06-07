#!/bin/bash

## For clinic dataset
outputDir=/home/robin/Documents/Project/Files/Clinic  #outputDir
bedFile=/home/robin/Documents/Project/Files/Clinic/clinic.bed #bedFile 
mlpaFile=/home/robin/Documents/Project/Files/Clinic/MLPA.xlsx #Validated file ICR96 format

Rscript transform_ICR96_validation_format.R --outputDir $outputDir --bedFile $bedFile --mlpaFile $mlpaFile

