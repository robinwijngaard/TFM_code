#!/bin/bash
# Description: Run the benchmark analysis for those algorithms and dataset defined in the algorithms.yaml and datasets.yaml file

# Read params
source ./utils/parse_yaml.sh
eval $(parse_yaml algorithms.yaml "pars_")

# create logs/output folder if not exists
mkdir -p logs
mkdir -p output

#Execute algorithms over selected datasets
#if [ "$pars_algorithms_panelcn" == "true" ]; then
#    echo "[$(date)] Executing panelcn.MOPS"
#	Rscript ./algorithms/panelcnmops/runPanelcnmops.r ./algorithms/panelcnmops/panelcnmopsParams.yaml datasets.yaml  > logs/panelcnmops.log 2>&1
#fi

#if [ "$pars_algorithms_decon" == "true" ]; then
#    echo "[$(date)] Executing DECoN"
#	Rscript ./algorithms/decon/runDecon.r ./algorithms/decon/deconParams.yaml datasets.yaml  > logs/decon.log 2>&1
#fi

#if [ "$pars_algorithms_exomedepth" == "true" ]; then
#    echo "[$(date)] Executing ExomeDepth"
#	Rscript ./algorithms/exomedepth/runExomedepth.r ./algorithms/exomedepth/exomedepthParams.yaml datasets.yaml  > logs/exomedepth.log 2>&1
#fi

#if [ "$pars_algorithms_codex2" == "true" ]; then
#    echo "[$(date)] Executing CODEX2"
#	Rscript ./algorithms/codex2/runCodex2.r ./algorithms/codex2/codex2Params.yaml datasets.yaml  > logs/codex2.log 2>&1
#fi

#if [ "$pars_algorithms_convading" == "true" ]; then
#    echo "[$(date)] Executing CoNVaDING"
#    export LC_ALL=en_US.UTF-8
#	Rscript ./algorithms/convading/runConvading.r ./algorithms/convading/convadingParams.yaml datasets.yaml  > logs/convading.log 2>&1
#fi

#if [ "$pars_algorithms_manta" == "true" ]; then
#    echo "[$(date)] Executing Manta"
#    Rscript ./algorithms/manta/runManta.r ./algorithms/manta/mantaParams.yaml datasets.yaml  > logs/manta.log 2>&1
#fi

if [ "$pars_algorithms_manta2" == "true" ]; then
    echo "[$(date)] Executing Manta2"
    Rscript ./algorithms/manta2/runManta2.r ./algorithms/manta2/mantaParams2.yaml datasets.yaml  > logs/manta2.log 2>&1
fi

if [ "$pars_algorithms_manta3" == "true" ]; then
    echo "[$(date)] Executing Manta3"
    Rscript ./algorithms/manta3/runManta3.r ./algorithms/manta3/mantaParams3.yaml datasets.yaml  > logs/manta3.log 2>&1
fi

#if [ "$pars_algorithms_cnvkit" == "true" ]; then
#    echo "[$(date)] Executing CNVkit1"
#    Rscript ./algorithms/cnvkit/runCNVkit.r ./algorithms/cnvkit/cnvkitParams.yaml datasets.yaml  > logs/cnvkit.log 2>&1
#fi

#if [ "$pars_algorithms_cnvkit2" == "true" ]; then
#    echo "[$(date)] Executing CNVkit2"
#    Rscript ./algorithms/cnvkit2/runCNVkit2.r ./algorithms/cnvkit2/cnvkitParams2.yaml datasets.yaml  > logs/cnvkit2.log 2>&1
#fi

#if [ "$pars_algorithms_cnvkit3" == "true" ]; then
#    echo "[$(date)] Executing CNVkit3"
#    Rscript ./algorithms/cnvkit3/runCNVkit3.r ./algorithms/cnvkit3/cnvkitParams3.yaml datasets.yaml  > logs/cnvkit3.log 2>&1
#fi

#if [ "$pars_algorithms_cnvkit4" == "true" ]; then
#    echo "[$(date)] Executing CNVkit4"
#    Rscript ./algorithms/cnvkit4/runCNVkit4.r ./algorithms/cnvkit4/cnvkitParams4.yaml datasets.yaml  > logs/cnvkit4.log 2>&1
#fi

#if [ "$pars_algorithms_cnvkit5" == "true" ]; then
#    echo "[$(date)] Executing CNVkit5"
#    Rscript ./algorithms/cnvkit5/runCNVkit5.r ./algorithms/cnvkit5/cnvkitParams5.yaml datasets.yaml  > logs/cnvkit5.log 2>&1
#fi

#Generate summary file
echo "[$(date)] Generating summary file"
Rscript ./utils/summary.r algorithms.yaml datasets.yaml  > logs/summary.log 2>&1