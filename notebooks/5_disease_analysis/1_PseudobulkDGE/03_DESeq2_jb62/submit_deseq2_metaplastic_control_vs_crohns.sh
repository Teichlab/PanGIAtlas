#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_deseq2_metaplastic_control_vs_crohns.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_deseq2_metaplastic_control_vs_crohns.e
# BSUB -J submit_deseq2_metaplastic_control_vs_crohns
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DATADIR="/nfs/team205/jb62/megagut/results/pseudobulk"
OUTDIR="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_metaplastic_cells"

if [ ! -d ${OUTDIR} ]; then
    mkdir ${OUTDIR}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/PSBULK_Metaplastic_cells.control_vs_crohns_count_matrix.csv \
    --metadata ${DATADIR}/PSBULK_Metaplastic_cells.control_vs_crohns_metadata.csv \
    --output-dir ${OUTDIR} \
    --output-prefix PSBULK_Metaplastic_cells.control_vs_crohns_ \
    --cell-type-column-id level_3_annot \
    --columns-for-design control_vs_disease_simple \
    --column-to-compare control_vs_disease_simple \
    --condition-of-interest IBD \
    --control-condition control