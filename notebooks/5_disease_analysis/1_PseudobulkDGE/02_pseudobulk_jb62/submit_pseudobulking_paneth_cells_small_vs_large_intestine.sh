#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_pseudobulking_paneth_cells_small_vs_large_intestine.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_pseudobulking_paneth_cells_small_vs_large_intestine.e
# BSUB -J submit_pseudobulking_paneth_cells_small_vs_large_intestine
# BSUB -R "select[mem>40GB] rusage[mem=40GB] span[hosts=1]"
# BSUB -M 40GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

python /lustre/scratch126/cellgen/team205/jb62/scripts/pseudobulk_generation/generate_pseudobulk_decoupler.py \
    --input /nfs/team205/ao15/Megagut/Annotations_v3/h5ad/DGE_analysis_objects/Paneth_cells.small_vs_large_intestine.4211cells_18485genes.20231114.h5ad \
    --output-prefix /nfs/team205/jb62/megagut/results/pseudobulk/PSBULK_Paneth_cells.small_vs_large_intestine \
    --min-number-cells 2 \
    --min-number-counts 0 \
    --donor-column-id donorID_unified \
    --cell-type-column-id level_3_annot \
    --keep-columns study,disease,organ_groups,organ_disease_simple \
    --no-raw-slot