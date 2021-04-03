#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate ~/anaconda3/envs/mantaenv

#paths

# Germline configuration examples

bams="/home/robin/Documents/Project/Samples/example/bam_ex"
bed="/home/robin/Documents/Project/Samples/bedfiles"
ref="/home/robin/Documents/Project/Samples/hg38"

manta_conf="/home/robin/Documents/Project/manta/Install"
manta_exe="/home/robin/Documents/Project/Results/Manta"


# Bamfiles

cd "$manta_exe"
ls /home/robin/Documents/Project/Samples/example/bam_ex/*.bam > bamfiles.txt

declare -a lines
readarray -t line <bamfiles.txt
echo "${line[0]}"

python "$manta_conf"/bin/configManta.py \
--bam `echo "${line[0]}"` \
--bam `echo "${line[1]}"` \
--bam `echo "${line[2]}"` \
--bam `echo "${line[3]}"` \
--bam `echo "${line[4]}"` \
--bam `echo "${line[5]}"` \
--bam `echo "${line[6]}"` \
--referenceFasta "$ref"/hg38.fa \
--runDir "$manta_exe"/"$1" \
--exome \
--callRegions "$bed"/ICR96_hg38.bed