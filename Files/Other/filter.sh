#!/bin/bash

bedfile="/home/robin/Documents/Project/Samples/bedfiles"
validatedfile="/home/robin/Documents/Project/Samples/validatedresults"
filter="/home/robin/Documents/Project/TFM_code/Files/Other"

#awk -F '\t' 'NR==FNR {A[$1];next} $3 in A {print $0}' "$filter"/geneList.txt "$validatedfile"/ICR96_validated_regions38.txt >"$validatedfile"/ICR96_validated_regions38_filter.txt
#awk -F '\t' 'NR==FNR {A[$1];next} $1 in A {print $0}' "$filter"/singlebams.txt "$validatedfile"/ICR96_validated_regions38.txt >"$validatedfile"/ICR96_validated_regions38_single.txt
#awk -F '\t' 'NR==FNR {A[$1];next} $1 in A {print $0}' "$filter"/singlebams.txt "$validatedfile"/ICR96_validated_regions38_filter.txt >"$validatedfile"/ICR96_validated_regions38_filter_single.txt
#awk -F '\t' 'FNR==NR {A[$1];next}$4 in A'  "$filter"/geneList.txt "$bedfile"/ICR96_hg38_noSNP.txt >"$bedfile"/ICR96_hg38_noSNP_filter.txt

head -1 "$validatedfile"/ICR96_validated_regions38.txt > "$validatedfile"/ICR96_validated_regions38_filter.txt
grep -f "$filter"/geneList.txt "$validatedfile"/ICR96_validated_regions38.txt >> "$validatedfile"/ICR96_validated_regions38_filter.txt

head -1 "$validatedfile"/ICR96_validated_regions38.txt > "$validatedfile"/ICR96_validated_regions38_single.txt
grep -f "$filter"/singlebams.txt "$validatedfile"/ICR96_validated_regions38.txt >> "$validatedfile"/ICR96_validated_regions38_single.txt

head -1 "$validatedfile"/ICR96_validated_regions38_filter.txt > "$validatedfile"/ICR96_validated_regions38_filter_single.txt
grep -f "$filter"/singlebams.txt "$validatedfile"/ICR96_validated_regions38_filter.txt >> "$validatedfile"/ICR96_validated_regions38_filter_single.txt


#grep -f "$filter"/geneList.txt "$validatedfile"/ICR96_validated_regions38.txt > "$validatedfile"/ICR96_validated_regions38_filter.txt

grep -f "$filter"/geneList.txt "$bedfile"/ICR96_hg38_noSNP.bed > "$bedfile"/ICR96_hg38_noSNP_filter.bed
