#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_pseudobulking_oral_mucosa.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_pseudobulking_oral_mucosa.e
# BSUB -J submit_pseudobulking_oral_mucosa
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

python /lustre/scratch126/cellgen/team205/jb62/scripts/pseudobulk_generation/generate_pseudobulk_decoupler.py \
    --input /nfs/team205/ao15/Megagut/Annotations_v3/h5ad/DGE_analysis_objects/Oral_mucosa_fibro.disease_comparison.8620cells_21932genes.20231206.h5ad \
    --output-prefix /nfs/team205/jb62/megagut/results/pseudobulk/PSBULK_Oral_mucosa_fibro.disease_comparison_ \
    --min-number-cells 2 \
    --min-number-counts 0 \
    --donor-column-id donorID_unified \
    --cell-type-column-id level_3_annot \
    --keep-columns study,disease,control_vs_disease_simple \
    --no-raw-slot