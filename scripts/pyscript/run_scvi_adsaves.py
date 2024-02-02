#!/usr/bin/env python
"""Run scvi from command line

Usage: run_scvi [options] <input> <out_prefix>

Options:
  --n-hvg <int>                number of highly variable genes to use [default: 5000]
  --remove-from-hvg <str>      gene sets to remove from HVG, must be boolean columns in `adata.var`
  --batch <str>                batch variable in `adata.obs`
  --min-batch-size <int>       batches smaller than this will be removed [default: 0]
  --categorical <str>          categorical variables in `adata.obs`
  --continuous <str>           continuous variables in `adata.obs`
  --n-latent <int>             number of latent variables to return [default: 20]
  --batch-size <int>           batch size for training [default: 512]
  --limit-train-batches <int>  number of batches checked each epoch
  --train-reference            use scArches' style of training
  --make-umap                  compute UMAP use scvi latent variables
  --n-neighbors <int>          number of nearest neighbors [default: 15]
  --make-plot                  save UMAP plot
  --color-by <str>             column in `adata.obs` that cells are colored by
  --debug                      print debug information
  --profile                    print profile information
  <input>                      input h5ad
  <out_prefix>                 output prefix

Used commands (reruns without categorical covariate [study]) 16th Jan 2023:

### Full object ###
./run_scvi_adsaves.py --n-hvg 7500 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pooled_healthy.gene_cellbender.good_qc_cluster_mito80.raw.subsetted.manual_doublet_removed.fine_annot.20230110.h5ad \
            ~/Annotations_v3/scvi_output20230116/pooled_healthy.hvg7500_noCC.nodoublets.scvi_output.20230113
            
./run_scvi_adsaves.py --n-hvg 10000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pooled_healthy.gene_cellbender.good_qc_cluster_mito80.raw.subsetted.manual_doublet_removed.fine_annot.20230110.h5ad \
            ~/Annotations_v3/scvi_output20230116/pooled_healthy.hvg10000_noCC.nodoublets.scvi_output.20230120
            
./run_scvi_adsaves.py --n-hvg 15000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pooled_healthy.gene_cellbender.good_qc_cluster_mito80.raw.subsetted.manual_doublet_removed.fine_annot.20230110.h5ad \
            ~/Annotations_v3/scvi_output20230116/pooled_healthy.hvg15000_noCC.nodoublets.scvi_output.20230120

### Compartments ###
ENDO
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Endo.h5ad \
            ~/Annotations_v3/scvi_output20230116/Endo.hvg5000_noCC.scvi_output.20230116

NEURAL 
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Neural.h5ad \
            ~/Annotations_v3/scvi_output20230116/Neural.hvg5000_noCC.scvi_output.20230116
            
            
T CELLS
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.T.h5ad \
            ~/Annotations_v3/scvi_output20230116/T.hvg5000_noCC.scvi_output.20230116
            
B CELLS           
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv,feature_list/Ig_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.B.h5ad \
            ~/Annotations_v3/scvi_output20230116/B.hvg5000_noCC_noIgVar.scvi_output.20230116

MYELOID
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Myeloid.h5ad \
            ~/Annotations_v3/scvi_output20230116/Myeloid.hvg5000_noCC.scvi_output.20230116
            
MESENCHMYE SUBSETTED BY AGE
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Mes_AP.h5ad \
            ~/Annotations_v3/scvi_output20230116/Mes_AP.hvg5000_noCC.scvi_output.20230116
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Mes_FT.h5ad \
            ~/Annotations_v3/scvi_output20230116/Mes_FT.hvg5000_noCC.scvi_output.20230116
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Mes_ST.h5ad \
            ~/Annotations_v3/scvi_output20230116/Mes_ST.hvg5000_noCC.scvi_output.20230116

EPITHELIAL   
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221206.Epi_all.donorA33_removed.h5ad \
            ~/Annotations_v3/scvi_output20230116/Epi_all.hvg5000_noRPL.scvi_output.20230116


EPITHELIAL SUBSETTED BY AGE/REGION
ORAL MUCOSA
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Oral_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/Oral_epi.hvg5000_noRPL.scvi_output.20230116
            
SALIVARY GLANDS
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Salivary_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/Salivary_epi.hvg5000_noRPL.scvi_output.20230116
            
OESOPHAGUS
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Oesophagus_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/Oesophagus_epi.hvg5000_noRPL.scvi_output.20230116

STOMACH
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Stomach_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/Stomach_epi.hvg5000_noRPL.scvi_output.20230116
            
SMALL INTESTINE
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_AP_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/SI_AP_epi.hvg5000_noRPL.scvi_output.20230116
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_FT_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/SI_FT_epi.hvg5000_noRPL.scvi_output.20230116
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_ST_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/SI_ST_epi.hvg5000_noRPL.scvi_output.20230116
            
LARGE INTESTINE
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.LI_AP_epi.donorA33_removed.h5ad \
            ~/Annotations_v3/scvi_output20230116/LI_AP_epi.hvg5000_noRPL.scvi_output.20230116
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.LI_FT_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/LI_FT_epi.hvg5000_noRPL.scvi_output.20230116
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.LI_ST_epi.h5ad \
            ~/Annotations_v3/scvi_output20230116/LI_ST_epi.hvg5000_noRPL.scvi_output.20230116
            

### Compartments WITH DONOR CORRECTED 26/01/2023 ###
ENDO
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Endo.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Endo.hvg5000_noCC.scvi_output.20230126

NEURAL 
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Neural.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Neural.hvg5000_noCC.scvi_output.20230126
            
            
T CELLS
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.T.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/T.hvg5000_noCC.scvi_output.20230126
            
B CELLS           
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv,feature_list/Ig_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.B.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/B.hvg5000_noCC_noIgVar.scvi_output.20230126

MYELOID
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Myeloid.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Myeloid.hvg5000_noCC.scvi_output.20230126
            
MESENCHMYE SUBSETTED BY AGE
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Mes_AP.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Mes_AP.hvg5000_noCC.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Mes_FT.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Mes_FT.hvg5000_noCC.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Mes_ST.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Mes_ST.hvg5000_noCC.scvi_output.20230126

EPITHELIAL   
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221206.Epi_all.donorA33_removed.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Epi_all.hvg5000_noRPL.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221206.Epi_adult.donorA33_removed.h5ad \
            ~/Annotations_v3/scvi_output20230126/Epi_adult.hvg5000_noRPL.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221206.Epi_fetal.donorA33_removed.h5ad \
            ~/Annotations_v3/scvi_output20230126/Epi_dev.hvg5000_noRPL.scvi_output.20230126


EPITHELIAL SUBSETTED BY AGE/REGION
ORAL MUCOSA
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Oral_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Oral_epi.hvg5000_noRPL.scvi_output.20230126
            
SALIVARY GLANDS
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Salivary_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Salivary_epi.hvg5000_noRPL.scvi_output.20230126
            
OESOPHAGUS
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Oesophagus_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Oesophagus_epi.hvg5000_noRPL.scvi_output.20230126

STOMACH
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.Stomach_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Stomach_epi.hvg5000_noRPL.scvi_output.20230126
            
SMALL INTESTINE
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_AP_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/SI_AP_epi.hvg5000_noRPL.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_FT_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/SI_FT_epi.hvg5000_noRPL.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_ST_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/SI_ST_epi.hvg5000_noRPL.scvi_output.20230126
            
LARGE INTESTINE
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.LI_AP_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/LI_AP_epi.hvg5000_noRPL.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.LI_FT_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/LI_FT_epi.hvg5000_noRPL.scvi_output.20230126
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.LI_ST_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/LI_ST_epi.hvg5000_noRPL.scvi_output.20230126 
            
./run_scvi_adsaves.py --n-hvg 7500 \
            --remove-from-hvg feature_list/cc_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/pooled_healthy.gene_cellbender.good_qc_cluster_mito80.raw.subsetted.manual_doublet_removed.fine_annot.20230110.donorIDfixed.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/pooled_healthy.hvg7500_noCC.nodoublets.scvi_output.20230126
            
            
GROUPED EPITHELIAL OBJECTS:

./run_scvi_adsaves.py --n-hvg 7500 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_LI_combined_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Epi_SI_LI_adult.hvg7500_noRPL.scvi_output.20230126
            
            
./run_scvi_adsaves.py --n-hvg 5000 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_LI_combined_epi.20230126.18485genes.h5ad \
            ~/Annotations_v3/scvi_output20230126/Epi_SI_LI_adult.hvg5000_noRPL.scvi_output.20230218
            

./run_scvi_adsaves.py --n-hvg 7500 \
            --remove-from-hvg feature_list/RPL_genes.csv \
            --batch donorID_unified \
            --min-batch-size 0 \
            --continuous log1p_n_counts,percent_mito \
            --train-reference \
            ~/Annotations_v3/pre_scvi20221124_updated_donorID_20230126/with_doublets/pooled_healthy.lv20_batch256.with_broad_annotation.20220917.with_countlayers.20221124.SI_LI_stomach_combined_epi.20230126.h5ad \
            ~/Annotations_v3/scvi_output20230126/Epi_SI_LI_stomach_adult.hvg7500_noRPL.scvi_output.20230126



Note:
    It shouldn't take long to run, so submitting to gpu-normal is sufficient.

    --make-umap is slow for large datasets and doesn't make use of GPU, so
    better not to enable it and instead do your own UMAP calculation on a cpu
    notebook.

"""


import logging
import signal
import sys
import os

import numpy as np
import pandas as pd
import scanpy as sc
import sctk as sk
import scvi

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


ARCHES_PARAMS = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
)


def main(args):
    logging.debug(args)

    # Get input options
    input_h5ad = args["input"]
    output_prefix = args["out_prefix"]
    n_hvg = int(args["n-hvg"])
    remove_genes = (
        args["remove-from-hvg"].split(",")
        if args["remove-from-hvg"] is not None
        else None
    )
    batch_key = args["batch"]
    min_batch_size = int(args["min-batch-size"])
    categorical_keys = (
        args["categorical"].split(",") if args["categorical"] is not None else None
    )
    continuous_keys = (
        args["continuous"].split(",") if args["continuous"] is not None else None
    )
    n_latent = int(args["n-latent"])
    batch_size = int(args["batch-size"]) if args["n-latent"] is not None else None
    limit_train_batches = (
        float(args["limit-train-batches"])
        if args["limit-train-batches"] is not None
        else None
    )
    train_reference = args["train-reference"]
    make_umap = args["make-umap"]
    n_neighbors = int(args["n-neighbors"])
    make_plot = args["make-plot"]
    color_by = args["color-by"]

    # Read input
    ad = sc.read(input_h5ad)

    if min_batch_size > 1:
        batch_sizes = ad.obs[batch_key].value_counts()
        batch_to_keep = batch_sizes.index[batch_sizes >= min_batch_size].astype(str).values
        ad = ad[ad.obs[batch_key].isin(batch_to_keep)].copy()

    # Subset to HVG
    sc.pp.highly_variable_genes(
        ad,
        flavor="seurat_v3",
        n_top_genes=n_hvg,
        subset=False,
    )
    ad1 = ad[:, ad.var.highly_variable].copy()

    # Remove certain gene sets if required
    if remove_genes:
        k_remove = np.zeros(ad1.n_vars).astype(bool)
        for g in remove_genes:
            if g in ad1.var.columns:
                k_remove = k_remove | ad1.var[g].values
            elif os.path.exists(g):
                genes = sk.read_list(g)
                k_remove = k_remove | ad1.var_names.isin(genes)
        ad1 = ad1[:, ~k_remove].copy()

    # Run scvi
    scvi.model.SCVI.setup_anndata(
        ad1,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_keys,
        continuous_covariate_keys=continuous_keys,
    )
    if limit_train_batches is None:
        if ad1.n_obs / batch_size <= 20:
            limit_train_batches = 1.0
        elif ad1.n_obs / batch_size >= 100:
            limit_train_batches = 0.2
        else:
            limit_train_batches = 20
    extra_params = ARCHES_PARAMS if train_reference else {}
    vae = scvi.model.SCVI(
        ad1, n_layers=2, dropout_rate=0.2, n_latent=n_latent, **extra_params
    )
    vae.train(
        train_size=0.9,
        early_stopping_patience=30,
        max_epochs=400,
        batch_size=batch_size,
        limit_train_batches=limit_train_batches,
        use_gpu=True,
    )
    ad.obsm["X_scvi"] = vae.get_latent_representation()

    vae.save(output_prefix + ".vae", overwrite=True)
    #ad.write(output_prefix + "-ad.h5ad", compression="gzip")
    ad1.write(output_prefix + "-ad1.h5ad",compression="gzip")
    if make_umap:
        sc.pp.neighbors(
            ad,
            use_rep="X_scvi",
            n_pcs=ad.obsm["X_scvi"].shape[1],
            n_neighbors=n_neighbors,
            key_added="neighbors_scvi",
        )
        sc.tl.umap(ad, neighbors_key="neighbors_scvi", min_dist=0.1)
        ad.write(output_prefix + ".h5ad", compression="gzip")
        if make_plot:
            sk.set_figsize((6, 6))
            ax = sc.pl.umap(ad, color=color_by, legend_loc="on data", show=False)
            ax.get_figure().savefig(output_prefix + ".umap.png", bbox_inches="tight")
            sk.set_figsize((4, 3))
            plt.plot(vae.history["elbo_train"])
            plt.savefig(output_prefix + ".elbo.png", bbox_inches="tight")
    else:
        X_scvi = pd.DataFrame(
            ad.obsm["X_scvi"],
            index=ad.obs_names,
            columns=[f"LV{(i+1):02d}" for i in range(n_latent)],
        )
        X_scvi.to_csv(output_prefix + ".csv.gz")

    return 0


if __name__ == "__main__":
    from docopt import docopt

    args = docopt(__doc__)
    args = {k.lstrip("-<").rstrip(">"): args[k] for k in args}
    try:
        if args.get("debug"):
            logLevel = logging.DEBUG
        else:
            logLevel = logging.WARN
        logging.basicConfig(
            level=logLevel,
            format="%(asctime)s; %(levelname)s; %(funcName)s; %(message)s",
            datefmt="%y-%m-%d %H:%M:%S",
        )
        if args.get("profile"):
            import cProfile

            cProfile.run("main(args)")
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning("Interrupted")
        sys.exit(1)
