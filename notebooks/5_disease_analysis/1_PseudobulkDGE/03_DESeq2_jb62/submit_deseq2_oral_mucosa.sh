#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_deseq2_oral_mucosa.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_deseq2_oral_mucosa.e
# BSUB -J submit_deseq2_oral_mucosa
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DATADIR="/nfs/team205/jb62/megagut/results/pseudobulk"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_oral_mucosa"

if [ ! -d ${OUTDIR1} ]; then
    mkdir ${OUTDIR1}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_Oral_mucosa_fibro.disease_comparison__count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_Oral_mucosa_fibro.disease_comparison__metadata.csv \
    --output-dir ${OUTDIR1} \
    --output-prefix PSBULK_Oral_mucosa_fibro.disease_comparison_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design control_vs_disease_simple \
    --column-to-compare control_vs_disease_simple \
    --condition-of-interest IBD \
    --control-condition control

OUTDIR2="/nfs/team205/jb62/megagut/results/dge/periodontitis_vs_control_in_oral_mucosa"

if [ ! -d ${OUTDIR2} ]; then
    mkdir ${OUTDIR2}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_Oral_mucosa_fibro.disease_comparison__count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_Oral_mucosa_fibro.disease_comparison__metadata.csv \
    --output-dir ${OUTDIR2} \
    --output-prefix PSBULK_Oral_mucosa_fibro.disease_comparison_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design control_vs_disease_simple \
    --column-to-compare control_vs_disease_simple \
    --condition-of-interest periodontitis \
    --control-condition control

OUTDIR3="/nfs/team205/jb62/megagut/results/dge/IBD_vs_periodontitis_in_oral_mucosa"

if [ ! -d ${OUTDIR3} ]; then
    mkdir ${OUTDIR3}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_Oral_mucosa_fibro.disease_comparison__count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_Oral_mucosa_fibro.disease_comparison__metadata.csv \
    --output-dir ${OUTDIR3} \
    --output-prefix PSBULK_Oral_mucosa_fibro.disease_comparison_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design control_vs_disease_simple \
    --column-to-compare control_vs_disease_simple \
    --condition-of-interest IBD \
    --control-condition periodontitis