#!/bin/bash

# Input files
bamsDir=/home/robin/Documents/Project/Samples/example/bam_ex
bedFile=/home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38_noSNP.bed # clinic.bed for clinic dataset
fastaFile=/home/robin/Documents/Project/Samples/hg38/hg38.fa

# Output dirs
resultDir=/home/robin/Documents/Project/TFM_code/Implementation/results

# Parameters
run=R1234 #run name
readLength=150 # sequence length

# Optimitzable parametres
phi_bins=1 
transition_probability=0.0001  
expected_CNV_length=50000

# For sample report
annotatedFile=/home/robin/Documents/Project/TFM_code/Implementation/annotate/annotatedFile_ICR96.bed
panellsFile=/home/robin/Documents/Project/TFM_code/Implementation/files/gens_nm_panell_ICR96.xlsx
markdownfile=/home/robin/Documents/Project/TFM_code/Implementation/scripts/report.Rmd 
samplesFile=/home/robin/Documents/Project/TFM_code/Implementation/files/samples.txt 


cd /home/robin/Documents/Project/TFM_code/Implementation/scripts/

Rscript runExomedepth.r --bamsDir $bamsDir --bedFile $bedFile --fastaFile $fastaFile --resultDir $resultDir --run $run \
--readLength $readLength --phi_bins $phi_bins --transition_probability $transition_probability --expected_CNV_length $expected_CNV_length \

Rscript runReport.r --resultDir $resultDir --run $run --annotatedFile $annotatedFile --panellsFile $panellsFile \
--markdownfile $markdownfile --samplesFile $samplesFile