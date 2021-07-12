#!/bin/bash

# Files
utils=/route/to/utils.r  #utils R script
bamsDir=/route/to/bamdirs
bedFile=/route/to/bedFile
fastaFile=/route/to/hg38

# Parameters
run=R1234 #run name
readLength=150 # sequence length

# Optimitzable parametres
phi_bins=1 
transition_probability=0.0001  
expected_CNV_length=50000

Rscript evaluation.r --outputDir $outputDir --bedFile $bedFile --roisFile $roisFile --cnvFounds $cnvFounds