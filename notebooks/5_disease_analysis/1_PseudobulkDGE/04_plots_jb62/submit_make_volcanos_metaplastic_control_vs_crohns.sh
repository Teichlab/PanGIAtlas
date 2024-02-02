#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_make_volcanos_metaplastic_control_vs_crohns.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_make_volcanos_metaplastic_control_vs_crohns.e
# BSUB -J submit_make_volcanos_metaplastic_control_vs_crohns
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DDS_DIR="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_metaplastic_cells/degs"
METADATA="/nfs/team205/jb62/megagut/results/pseudobulk/PSBULK_Metaplastic_cells.control_vs_crohns_metadata.csv"
OUTDIR="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_metaplastic_cells"


Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results ${DDS_DIR}/control_vs_disease_simple/PSBULK_Metaplastic_cells.control_vs_crohns_DESeq2_control_vs_disease_simple_IBD_vs_control_level_3_annot_design_control_vs_disease_simple_DEGs_noFilter.csv \
    --deseq2-design control_vs_disease_simple \
    --deseq2-comparison control_vs_disease_simple_IBD_vs_control \
    --metadata ${METADATA} \
    --cell-type-metadata level_3_annot \
    --donor-metadat donorID_unified \
    --lfc-column-id log2FoldChange \
    --lower-lfc-threshold -0.5 \
    --upper-lfc-threshold 0.5 \
    --padj-column-id padj \
    --cell-type-dds cell_type \
    --volcano-width 10.0 \
    --volcano-height 6.0 \
    --gene-column-id gene \
    --output-dir ${OUTDIR}