#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_make_volcanos_muc6_muc5ac_stomach_smallintestine_organ_disease_comparison.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_make_volcanos_muc6_muc5ac_stomach_smallintestine_organ_disease_comparison.e
# BSUB -J submit_make_volcanos_muc6_muc5ac_stomach_smallintestine_organ_disease_comparison
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DDS_DIR1="/nfs/team205/jb62/megagut/results/dge/SmallInstestineIBD_vs_StomachControl_MUC6_MUC5AC/degs"
METADATA="/nfs/team205/jb62/megagut/results/pseudobulk/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_metadata.csv"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/SmallInstestineIBD_vs_StomachControl_MUC6_MUC5AC"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results ${DDS_DIR1}/disease_organ_simple/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_DESeq2_disease_organ_simple_Small_intestine_IBD_vs_Stomach_control_level_3_annot_design_disease_organ_simple_DEGs_noFilter.csv \
    --deseq2-design disease_organ_simple \
    --deseq2-comparison disease_organ_simple_Small_intestine_IBD_vs_Stomach_control \
    --metadata ${METADATA} \
    --cell-type-metadata level_3_annot \
    --donor-metadata donorID_unified \
    --lfc-column-id log2FoldChange \
    --lower-lfc-threshold -0.5 \
    --upper-lfc-threshold 0.5 \
    --padj-column-id padj \
    --cell-type-dds cell_type \
    --volcano-width 10.0 \
    --volcano-height 6.0 \
    --gene-column-id gene \
    --output-dir ${OUTDIR1}

DDS_DIR2="/nfs/team205/jb62/megagut/results/dge/SmallIntestineControl_vs_StomachControl_MUC6_MUC5AC/degs"
OUTDIR2="/nfs/team205/jb62/megagut/results/dge/SmallIntestineControl_vs_StomachControl_MUC6_MUC5AC"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results ${DDS_DIR2}/disease_organ_simple/PSBULK_MUC6_MUC5AC.stomach_smallintestine.organ_disease_comparison_DESeq2_disease_organ_simple_Small_intestine_control_vs_Stomach_control_level_3_annot_design_disease_organ_simple_DEGs_noFilter.csv \
    --deseq2-design disease_organ_simple \
    --deseq2-comparison disease_organ_simple_Small_intestine_control_vs_Stomach_control \
    --metadata ${METADATA} \
    --cell-type-metadata level_3_annot \
    --donor-metadata donorID_unified \
    --lfc-column-id log2FoldChange \
    --lower-lfc-threshold -0.5 \
    --upper-lfc-threshold 0.5 \
    --padj-column-id padj \
    --cell-type-dds cell_type \
    --volcano-width 10.0 \
    --volcano-height 6.0 \
    --gene-column-id gene \
    --output-dir ${OUTDIR2}