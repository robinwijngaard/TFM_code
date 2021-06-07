#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96/exomedepth  #outputDir
resultFile=/home/robin/Documents/Project/Files/ICR96/results/resultData.txt #bedFile 
cnvFounds=/home/robin/Documents/Project/Files/ICR96/cnvFounds/cnvFounds_exomedepth.txt #

Rscript exomedepth.r --outputDir $outputDir --resultFile $resultFile --cnvFounds $cnvFounds

## For Clinic dataset
outputDir=/home/robin/Documents/Project/Files/Clinic/exomedepth  #outputDir
resultFile=/home/robin/Documents/Project/Files/Clinic/results/resultData.txt #bedFile 
cnvFounds=/home/robin/Documents/Project/Files/Clinic/cnvFounds/cnvFounds_exomedepth.txt #

Rscript exomedepth.r --outputDir $outputDir --resultFile $resultFile --cnvFounds $cnvFounds