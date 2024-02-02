#!/bin/bash
set -e pipefail

#navigate to this folder, activate environment with docopt for jsub
#/software/team205/bin/jsub lsf -n cpdb -q normal -c 20 -m 150g "bash cpdb_run_pooleddisease18485genes_20230526.sh" | bsub -G teichlab

PATH="/software/singularity-v3.6.4/bin:$PATH"

META=/nfs/team205/ao15/Megagut/Annotations_v3/disease_analysis/interactions_cellphoneDB/meta/pooled_disease.remapped.18485genes.AP_SI.meta.csv
COUNTS=/nfs/team205/ao15/Megagut/Annotations_v3/disease_analysis/interactions_cellphoneDB/counts/pooled_disease.remapped.18485genes.AP_SI.counts_normalisedonly.csv
# must exist
OUTPUT_PATH=/nfs/team205/ao15/Megagut/Annotations_v3/disease_analysis/interactions_cellphoneDB/cpdb_output/pooled_disease_18485genes

# create output folder if it does not exist
if [[ ! -d "${OUTPUT_PATH}" ]]; then
    mkdir -p "${OUTPUT_PATH}"
fi

singularity run --bind /nfs,/lustre  \
  /nfs/cellgeni/singularity/images/cellphonedb-v3.0.2.sif \
  cellphonedb method statistical_analysis \
  $META $COUNTS \
  --output-path=$OUTPUT_PATH \
  --threads=20 \
  --counts-data=hgnc_symbol \
  --threshold 0.1 \
  --database "v4.0.0" \
  --verbose
