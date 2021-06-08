#!/bin/bash

## For ICR96 dataset
outputDir=/home/robin/Documents/Project/Files/ICR96  #outputDir
resultFile=/home/robin/Documents/Project/Files/ICR96/results/resultData.txt #results file (obtained from evaluation output)

Rscript graphs.r --outputDir $outputDir --resultFile $resultFile