#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_pseudobulking_stomach_small_intestine_organ_disease_comparison.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_pseudobulking_stomach_small_intestine_organ_disease_comparison.e
# BSUB -J submit_pseudobulking_stomach_small_intestine_organ_disease_comparison
# BSUB -R "select[mem>40GB] rusage[mem=40GB] span[hosts=1]"
# BSUB -M 40GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

python /lustre/scratch126/cellgen/team205/jb62/scripts/pseudobulk_generation/generate_pseudobulk_decoupler.py \
    --input /nfs/team205/ao15/Megagut/Annotations_v3/h5ad/DGE_analysis_objects/MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison.22574cells_18485genes.20231114.h5ad \
    --output-prefix /nfs/team205/jb62/megagut/results/pseudobulk/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison \
    --min-number-cells 2 \
    --min-number-counts 0 \
    --donor-column-id donorID_unified \
    --cell-type-column-id level_3_annot \
    --keep-columns study,disease,organ_groups,disease_organ_simple \
    --no-raw-slot