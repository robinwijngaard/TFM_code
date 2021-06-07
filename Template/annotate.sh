#!/bin/bash

outputDir=/home/robin/Documents/Project/Files/Template #outputDir
bedFile=/home/robin/Documents/Project/Files/ICR96/ICR96_hg38_noSNP.bed #bedFile 
panellsFile=/home/robin/Documents/Project/Files/Template/gens_nm_panell.xlsx

Rscript annotate.R --outputDir $outputDir --bedFile $bedFile --panellsFile $panellsFile
