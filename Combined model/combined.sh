#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  #outputDir
combinedFile=/home/robin/Documents/Project/Files/ICR96/results/combinedData.txt #combinedData file (output evaluation) 
inhouseFile=/home/robin/Documents/Project/Files/Clinic/results/combinedData.txt

Rscript evaluation.r --outputDir $outputDir --combinedFile $combinedFile --inhouseFile $inhouseFile

