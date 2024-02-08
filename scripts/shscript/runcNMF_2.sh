#!/bin/bash
# Submits cNMF jobs in parallel on farm.
# runtime: ~20m
#
# conda activate cnmf_env
# fork install of cNMF, fixing 0 count boolean vector error for sparse matrices:
# pip install -e git+https://github.com/skoplev/cNMF#egg=cnmf

# number of workers, split for each k to evaluate
n_workers=36

# Raw counts h5ad file
h5ad_file="/lustre/scratch126/cellgen/team205/sk29/megagut/cNMF_rebuttal/h5ad/pooled_healthy_disease.plus_additional_epi.18485genes.with_fineannot.nodoublets.20230322.small_largeintestine.h5ad"

name="megagutintestine_cNMF_v2"

mkdir -p logs  # folder for LSF logs, if not already

# prepare cNMF
# -K wait bsub, bash script continues only when finished
#
# Normalize gene counts and prepare run parameters
# -c path to h5ad file
# -k space separated list of K values to be tested
# --numgenes high variance genes to be used

# 80GB for full dataset
# numgenes 2000 -> 5000

bsub -q long -K \
    -n 2 \
    -M 40GB -R "select[mem>40GB] rusage[mem=40GB]" \
    -o logs/prepare_megagutintestine_output.log -e logs/prepare_megagutintestine_error.log \
    cnmf prepare \
        --output-dir ./results \
        --name $name \
        -c "$h5ad_file" \
        --n-iter 100 \
        -k {10..80..2} \
        --seed 14 \
        --total-workers "$n_workers" \
        --numgenes 2000


# Loop over jobs, running cnmf factorize
for (( i = 0; i < $n_workers; i++ )); do
    bsub -q basement \
        -n 2 \
        -M 6GB -R "select[mem>6GB] rusage[mem=6GB]" \
        -o logs/megagutintestine_output.log -e logs/megagutintestine_error.log \
        cnmf factorize \
            --output-dir ./results \
            --name $name \
            --worker-index "$i" \
            --total-workers "$n_workers"
done