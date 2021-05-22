#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate ~/anaconda3/envs/mantaenv

#paths
bams="/home/robin/Documents/Project/Samples/example/bam_ex"
bed="/home/robin/Documents/Project/Samples/bedfiles"
ref="/home/robin/Documents/Project/Samples/hg38"
manta_conf="/home/robin/Documents/Project/manta/Install"
manta_exe="/home/robin/Documents/Project/Results/Manta"

# number of open files
ulimit -n 75000

cd "$manta_exe"
ls /home/robin/Documents/Project/Samples/example/bam_ex/*.bam > bamfiles.txt

# array of samples
declare -a lines
readarray -t line <bamfiles.txt
echo "${line[0]}"

# manta configuration
python "$manta_conf"/bin/configManta.py \
--bam `echo "${line[0]}"` \
--bam `echo "${line[1]}"` \
--bam `echo "${line[2]}"` \
--bam `echo "${line[3]}"` \
--bam `echo "${line[4]}"` \
--bam `echo "${line[5]}"` \
--bam `echo "${line[6]}"` \
--referenceFasta "$ref"/hg38.fa \
--runDir "$manta_exe" \
--exome \
--callRegions "$bed"/ICR96_hg38_noSNP.bed.gz

# manta execution
python "$manta_exe"/runWorkflow.py 








#--bam `echo "${line[7]}"` \
#--bam `echo "${line[8]}"` \
#--bam `echo "${line[9]}"` \
#--bam `echo "${line[10]}"` \
#--bam `echo "${line[11]}"` \