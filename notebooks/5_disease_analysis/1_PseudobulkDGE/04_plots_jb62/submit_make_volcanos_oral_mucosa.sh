#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_make_volcanos_paneth_cells_small_vs_large_intestine.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_make_volcanos_paneth_cells_small_vs_large_intestine.e
# BSUB -J submit_make_volcanos_paneth_cells_small_vs_large_intestine
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DDS_DIR1="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_oral_mucosa/degs"
METADATA="/nfs/team205/jb62/megagut/results/pseudobulk/PSBULK_Oral_mucosa_fibro.disease_comparison__metadata.csv"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/IBD_vs_control_in_oral_mucosa"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results ${DDS_DIR1}/control_vs_disease_simple/PSBULK_Oral_mucosa_fibro.disease_comparison_DESeq2_control_vs_disease_simple_IBD_vs_control_level_3_annot_design_control_vs_disease_simple_DEGs_noFilter.csv \
    --deseq2-design control_vs_disease_simple \
    --deseq2-comparison control_vs_disease_simple_IBD_vs_control \
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

DDS_DIR2="/nfs/team205/jb62/megagut/results/dge/IBD_vs_periodontitis_in_oral_mucosa/degs"
OUTDIR2="/nfs/team205/jb62/megagut/results/dge/IBD_vs_periodontitis_in_oral_mucosa"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results ${DDS_DIR2}/control_vs_disease_simple/PSBULK_Oral_mucosa_fibro.disease_comparison_DESeq2_control_vs_disease_simple_IBD_vs_periodontitis_level_3_annot_design_control_vs_disease_simple_DEGs_noFilter.csv \
    --deseq2-design control_vs_disease_simple \
    --deseq2-comparison control_vs_disease_simple_IBD_vs_periodontitis \
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

DDS_DIR3="/nfs/team205/jb62/megagut/results/dge/periodontitis_vs_control_in_oral_mucosa/degs"
OUTDIR3="/nfs/team205/jb62/megagut/results/dge/periodontitis_vs_control_in_oral_mucosa"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results ${DDS_DIR3}/control_vs_disease_simple/PSBULK_Oral_mucosa_fibro.disease_comparison_DESeq2_control_vs_disease_simple_periodontitis_vs_control_level_3_annot_design_control_vs_disease_simple_DEGs_noFilter.csv \
    --deseq2-design control_vs_disease_simple \
    --deseq2-comparison control_vs_disease_simple_periodontitis_vs_control \
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
    --output-dir ${OUTDIR3}