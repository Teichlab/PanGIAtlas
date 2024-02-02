#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_deseq2_muc6_muc5ac_stomach_smallintestine_organ_disease_comparison.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_deseq2_muc6_muc5ac_stomach_smallintestine_organ_disease_comparison.e
# BSUB -J submit_deseq2_muc6_muc5ac_stomach_smallintestine_organ_disease_comparison
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DATADIR="/nfs/team205/jb62/megagut/results/pseudobulk"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/SmallInstestineIBD_vs_StomachControl_MUC6_MUC5AC"

if [ ! -d ${OUTDIR1} ]; then
    mkdir ${OUTDIR1}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_metadata.csv \
    --output-dir ${OUTDIR1} \
    --output-prefix PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design disease_organ_simple \
    --column-to-compare disease_organ_simple \
    --condition-of-interest Small_intestine_IBD \
    --control-condition Stomach_control


OUTDIR2="/nfs/team205/jb62/megagut/results/dge/SmallIntestineControl_vs_StomachControl_MUC6_MUC5AC"

if [ ! -d ${OUTDIR2} ]; then
    mkdir ${OUTDIR2}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_metadata.csv \
    --output-dir ${OUTDIR2} \
    --output-prefix PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design disease_organ_simple \
    --column-to-compare disease_organ_simple \
    --condition-of-interest Small_intestine_control \
    --control-condition Stomach_control