#!/bin/bash


eval "$(conda shell.bash hook)"
conda activate ~/anaconda3/envs/mantaenv

#paths

# Germline configuration examples

manta_conf="/home/robin/Documents/Project/manta/Install"
manta_exe="/home/robin/Documents/Project/Results/Manta"

python "$manta_exe"/"$1"/runWorkflow.py 