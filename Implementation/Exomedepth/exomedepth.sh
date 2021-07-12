#!/bin/bash

# Input files
utils=/home/robin/Documents/Project/TFM_code/Implementation/Exomedepth/utils.r  #utils R script
bamsDir=/home/robin/Documents/Project/Samples/example/bam_ex
bedFile=/home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38_noSNP.bed
fastaFile=/home/robin/Documents/Project/Samples/hg38/hg38.fa

# Output dirs
resultDir=/home/robin/Documents/Project/TFM_code/Implementation/Exomedepth/results


# Parameters
run=R1234 #run name
readLength=150 # sequence length

# Optimitzable parametres
phi_bins=1 
transition_probability=0.0001  
expected_CNV_length=50000

cd /home/robin/Documents/Project/TFM_code/Implementation/ExomeDepth/

Rscript runExomedepth.r --utils $utils --bamsDir $bamsDir --bedFile $bedFile --fastaFile $fastaFile --resultDir $resultDir --run $run \
--readLength $readLength --phi_bins $phi_bins --transition_probability $transition_probability --expected_CNV_length $expected_CNV_length