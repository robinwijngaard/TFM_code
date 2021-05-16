#!/bin/bash

bedfile="/home/robin/Documents/Project/Samples/bedfiles"
validatedfile="/home/robin/Documents/Project/Samples/validatedresults"
filter="/home/robin/Documents/Project/TFM_code/Files/Other"

# Filter validated results

#head -1 "$validatedfile"/ICR96_validated_regions38.txt > "$validatedfile"/ICR96_validated_regions38_clinic.txt
#grep -f "$filter"/geneList.txt "$validatedfile"/ICR96_validated_regions38.txt >> "$validatedfile"/ICR96_validated_regions38_clinic.txt

head -1 "$validatedfile"/ICR96_validated_regions38.txt > "$validatedfile"/ICR96_validated_regions38_single.txt
grep -f "$filter"/singlebams.txt "$validatedfile"/ICR96_validated_regions38.txt >> "$validatedfile"/ICR96_validated_regions38_single.txt

# Filter genes on bedfile

#grep -f "$filter"/geneList.txt "$bedfile"/ICR96_hg38_noSNP.bed > "$bedfile"/ICR96_hg38_noSNP_clinic.bed
