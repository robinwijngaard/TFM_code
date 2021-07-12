#!/bin/bash

outputDir=/home/robin/Documents/Project/TFM_code/Implementation/annotate #outputDir
bedFile=/home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38_noSNP.bed #bedFile 
panellsFile=/home/robin/Documents/Project/TFM_code/Implementation/files/gens_nm_panell_ICR96.xlsx #genes per gene panell

Rscript annotate.R --outputDir $outputDir --bedFile $bedFile --panellsFile $panellsFile
