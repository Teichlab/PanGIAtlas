#!/bin/bash

# BSUB -q normal
# BSUB -o /nfs/team205/jb62/megagut/logs/submit_deseq2_LCM_UACLbulk_GSE126299.o
# BSUB -e /nfs/team205/jb62/megagut/logs/submit_deseq2_LCM_UACLbulk_GSE126299.e
# BSUB -J submit_deseq2_LCM_UACLbulk_GSE126299
# BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
# BSUB -M 10GB

source /lustre/scratch126/cellgen/team205/jb62/miniconda3/etc/profile.d/conda.sh 
conda activate jb62-single-cell-env

DATADIR="/nfs/team205/ao15/Megagut/LCM_UACLbulk_GSE126299"
OUTDIR1="/nfs/team205/jb62/megagut/results/dge/LCM_UACLbulk_GSE126299_UACL_vs_nonUACL"

if [ ! -d ${OUTDIR1} ]; then
    mkdir ${OUTDIR1}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/GSE126299_featurecounts_upload.gene.quant.txt \
    --sample-location cols \
    --metadata ${DATADIR}/GSE126299_meta.csv \
    --metadata-row-names 2 \
    --output-dir ${OUTDIR1} \
    --output-prefix LCM_UACLbulk_GSE126299_GSE126299_UACL_vs_nonUACL_ \
    --cell-type-column-id "Enriched cell type" \
    --columns-for-design Sample_info \
    --column-to-compare Sample_info \
    --condition-of-interest "CD_UACL_Intestinal epithelium monlayer" \
    --control-condition "CD_non-UACL_Intestinal epithelium monlayer"


OUTDIR2="/nfs/team205/jb62/megagut/results/dge/LCM_UACLbulk_GSE126299_UACL_vs_healthy_control"

if [ ! -d ${OUTDIR2} ]; then
    mkdir ${OUTDIR2}
fi

Rscript /lustre/scratch126/cellgen/team205/jb62/scripts/dge_analysis/run_dge_deseq2_per_celltype.R \
    --count-matrix ${DATADIR}/GSE126299_featurecounts_upload.gene.quant.txt \
    --sample-location cols \
    --metadata ${DATADIR}/GSE126299_meta.csv \
    --metadata-row-names 2 \
    --output-dir ${OUTDIR2} \
    --output-prefix LCM_UACLbulk_GSE126299_UACL_vs_healthy_control_ \
    --cell-type-column-id "Enriched cell type" \
    --columns-for-design Sample_info \
    --column-to-compare Sample_info \
    --condition-of-interest "CD_UACL_Intestinal epithelium monlayer" \
    --control-condition "healthy control_non-UACL_Intestinal epithelium monlayer"
