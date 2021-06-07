#!/bin/bash

outputDir=/home/robin/Documents/Project/Files/Template #outputDir
annotatedFile=/home/robin/Documents/Project/Files/Template/annotatedFile.bed 
panellsFile=/home/robin/Documents/Project/Files/Template/gens_nm_panell.xlsx
samplesFile=/home/robin/Documents/Project/Files/Template/samples.txt
RDSfolder=/home/robin/Documents/Project/Files/Template/RDSfiles
markdownfile=/home/robin/Documents/Project/TFM_code/Template/report.Rmd

Rscript report.R --outputDir $outputDir --annotatedFile $annotatedFile --panellsFile $panellsFile --samplesFile $samplesFile --RDSfolder $RDSfolder --markdownfile $markdownfile