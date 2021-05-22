#!/bin/bash


# Enter in environment

eval "$(conda shell.bash hook)"
conda activate ~/anaconda3/envs/cnvkit

# Control: if v = 3.8.8, OK (normal is 3.8.5)
python --version

rundir="/home/robin/Documents/Project/TFM_code/cnvkit-master/"

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
targets="/home/robin/Documents/Project/Samples/bedfiles"
out="/home/robin/Documents/Project/Results/cnvkit"
ref="/home/robin/Documents/Project/Samples/hg38"

cd "$rundir"

cnvkit.py target "$targets"/ICR96_hg38.bed --split -o "$out"/targets.split.bed

# Create a blank file to substitute for antitargets
touch "$out"/MT

# For each sample

cd "$bams"

for bam in `ls *.bam`
do
	name=`echo "$bam" | sed 's/\..*//g'`
	echo "$name"

	cd "$rundir"
	cnvkit.py coverage "$bams"/"$bam" "$out"/targets.split.bed -o "$out"/"$name".targetcoverage.cnn
done


cd "$rundir"
cnvkit.py reference "$out"/*.targetcoverage.cnn -f "$ref"/hg38.fa --no-edge -o "$out"/ref-tas.cnn

cd "$out"

for sample in `ls *.targetcoverage.cnn`
do
	name=`echo "$sample" | sed 's/.targetcoverage.cnn//g'`
	echo "$name"

	cd "$rundir"
	cnvkit.py fix "$out"/"$sample" "$out"/MT "$out"/ref-tas.cnn --no-edge -o "$out"/"$name".cnr
	cnvkit.py segment "$out"/"$name".cnr -o "$out"/"$name".cns
done