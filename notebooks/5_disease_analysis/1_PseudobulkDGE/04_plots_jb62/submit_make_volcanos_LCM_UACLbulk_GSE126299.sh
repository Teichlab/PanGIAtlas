#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_make_volcanos_LCM_UACLbulk_GSE126299.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_make_volcanos_LCM_UACLbulk_GSE126299.e
# BSUB -J submit_make_volcanos_LCM_UACLbulk_GSE126299
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DDS_DIR1="/nfs/team205/jb62/megagut/results/dge/LCM_UACLbulk_GSE126299_UACL_vs_healthy_control/degs"
METADATA="/nfs/team205/ao15/Megagut/LCM_UACLbulk_GSE126299/GSE126299_meta.csv"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/LCM_UACLbulk_GSE126299_UACL_vs_healthy_control"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results "${DDS_DIR1}/Sample_info/LCM_UACLbulk_GSE126299_UACL_vs_healthy_control_DESeq2_Sample_info_CD_UACL_Intestinal_epithelium_monlayer_vs_healthy_control_non_UACL_Intestinal_epithelium_monlayer_Enriched cell type_design_Sample_info_DEGs_noFilter.csv" \
    --deseq2-design Sample_info \
    --deseq2-comparison Sample_info_CD_UACL_Intestinal_epithelium_monlayer_vs_healthy_control_non_UACL_Intestinal_epithelium_monlayer \
    --metadata ${METADATA} \
    --metadata-row-names 1 \
    --cell-type-metadata "Enriched cell type" \
    --donor-metadata Sample_ID \
    --lfc-column-id log2FoldChange \
    --lower-lfc-threshold -0.5 \
    --upper-lfc-threshold 0.5 \
    --padj-column-id padj \
    --cell-type-dds cell_type \
    --volcano-width 6.0 \
    --volcano-height 6.0 \
    --gene-column-id gene \
    --output-dir ${OUTDIR1}

DDS_DIR2="/nfs/team205/jb62/megagut/results/dge/LCM_UACLbulk_GSE126299_UACL_vs_nonUACL/degs"
OUTDIR2="/nfs/team205/jb62/megagut/results/dge/LCM_UACLbulk_GSE126299_UACL_vs_nonUACL"

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/make_deseq2_plots.R \
    --deseq2-results "${DDS_DIR2}/Sample_info/LCM_UACLbulk_GSE126299_GSE126299_UACL_vs_nonUACL_DESeq2_Sample_info_CD_UACL_Intestinal_epithelium_monlayer_vs_CD_non_UACL_Intestinal_epithelium_monlayer_Enriched cell type_design_Sample_info_DEGs_noFilter.csv" \
    --deseq2-design Sample_info \
    --deseq2-comparison Sample_info_CD_UACL_Intestinal_epithelium_monlayer_vs_CD_non_UACL_Intestinal_epithelium_monlayer \
    --metadata ${METADATA} \
    --metadata-row-names 1 \
    --cell-type-metadata "Enriched cell type" \
    --donor-metadata Sample_ID \
    --lfc-column-id log2FoldChange \
    --lower-lfc-threshold -0.5 \
    --upper-lfc-threshold 0.5 \
    --padj-column-id padj \
    --cell-type-dds cell_type \
    --volcano-width 6.0 \
    --volcano-height 6.0 \
    --gene-column-id gene \
    --output-dir ${OUTDIR2}