#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_deseq2_paneth_cells_small_vs_large_intestine.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_deseq2_paneth_cells_small_vs_large_intestine.e
# BSUB -J submit_deseq2_paneth_cells_small_vs_large_intestine
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DATADIR="/nfs/team205/jb62/megagut/results/pseudobulk"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/LargeInstestineMetaplasia_vs_SmallIntestineControl_in_paneth_cells"

if [ ! -d ${OUTDIR1} ]; then
    mkdir ${OUTDIR1}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_Paneth_cells.small_vs_large_intestine_count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_Paneth_cells.small_vs_large_intestine_metadata.csv \
    --output-dir ${OUTDIR1} \
    --output-prefix PSBULK_Paneth_cells.small_vs_large_intestine_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design organ_disease_simple \
    --column-to-compare organ_disease_simple \
    --condition-of-interest Large_intestine_metaplasia \
    --control-condition Small_intestine_control


OUTDIR2="/nfs/team205/jb62/megagut/results/dge/LargeInstestineMetaplasia_vs_SmallIntestineNeighbouringInflamed_in_paneth_cells"

if [ ! -d ${OUTDIR2} ]; then
    mkdir ${OUTDIR2}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_Paneth_cells.small_vs_large_intestine_count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_Paneth_cells.small_vs_large_intestine_metadata.csv \
    --output-dir ${OUTDIR2} \
    --output-prefix PSBULK_Paneth_cells.small_vs_large_intestine_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design organ_disease_simple \
    --column-to-compare organ_disease_simple \
    --condition-of-interest Large_intestine_metaplasia \
    --control-condition Small_intestine_neighbouring_inflamed