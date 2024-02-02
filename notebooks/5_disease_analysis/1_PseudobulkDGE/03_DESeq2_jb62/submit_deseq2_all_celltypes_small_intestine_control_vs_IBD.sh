#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_deseq2_all_celltypes_small_intestine_control_vs_IBD.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_deseq2_all_celltypes_small_intestine_control_vs_IBD.e
# BSUB -J submit_deseq2_all_celltypes_small_intestine_control_vs_IBD
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DATADIR="/nfs/team205/jb62/megagut/results/pseudobulk"
OUTDIR="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_all_celltypes_of_small_intestine"

if [ ! -d ${OUTDIR} ]; then
    mkdir ${OUTDIR}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_All_celltypes.small_intestine.control_vs_IBD_count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_All_celltypes.small_intestine.control_vs_IBD_metadata.csv \
    --output-dir ${OUTDIR} \
    --output-prefix PSBULK_All_celltypes.small_intestine.control_vs_IBD_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design study,control_vs_disease_simple \
    --column-to-compare control_vs_disease_simple \
    --condition-of-interest IBD \
    --control-condition control