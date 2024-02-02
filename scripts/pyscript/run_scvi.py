#!/usr/bin/env python
"""Run scVI/scANVI from command line

Usage:
  run_scvi train [option] <input> <out_prefix>
  run_scvi query [option] <input> <out_prefix>
  run_scvi integrate [option] <input> <out_prefix>

Options:
  --use-model <path>           trained scVI/scANVI model
  --n-hvg <int>                number of highly variable genes to use [default: 5000]
  --remove-from-hvg <str>      gene sets to remove from HVG, must be boolean columns in `adata.var` or one gene per row files
  --batch <str>                batch variable in `adata.obs` to correct for
  --min-batch-size <int>       batches smaller than this will be removed [default: 0]
  --train-reference            use scArches' style of training
  --use-annotation <str>       use this column in `adata.obs` as annotation and trains scANVI model
  --categorical <str>          categorical variables in `adata.obs`, ignored if --train-reference is specified
  --continuous <str>           continuous variables in `adata.obs`
  --n-latent <int>             number of latent variables to return [default: 20]
  --batch-size <int>           batch size for training [default: 256]
  --limit-train-batches <int>  number of batches checked each epoch
  --make-embed                 compute MDE use scVI/scANVI latent variables
  --n-neighbors <int>          number of nearest neighbors [default: 15]
  --debug                      print debug information
  --profile                    print profile information
  <input>                      path to input h5ad or a list of h5ad
  <out_prefix>                 output prefix

# Example:
#     run_scvi.py train \\
#             --n-hvg 5000 \\
#             --remove-from-hvg feature_list/JP_cycle_genes.list,feature_list/ig_genes.list,feature_list/tcr_genes.list \\
#             --batch donor \\
#             --min-batch-size 30 \\
#             --categorical study \\
#             --continuous log1p_n_counts,percent_mito \\
#             --train-reference \\
#             h5ad/compartment/pooled_healthy.gene_cellbender.qc_filtered.post_scrublet.B_Plasma.processed.h5ad \\
#             h5ad/compartment/pooled_healthy.gene_cellbender.qc_filtered.post_scrublet.B_Plasma.scvi_output

# Note:
#     It shouldn't take long to run, so submitting to gpu-normal is sufficient.
# 
#     --make-embed only makes MDE. To make UMAP, do it on a cpu notebook.

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
    print(args)

    # Get input options
    input_h5ad = args["input"]
    output_prefix = args["out_prefix"]
    use_model = args["use-model"]
    n_hvg = int(args["n-hvg"])
    remove_genes = (
        args["remove-from-hvg"].split(",")
        if args["remove-from-hvg"] is not None
        else None
    )
    batch_key = args["batch"]
    min_batch_size = int(args["min-batch-size"])
    train_reference = args["train-reference"]
    use_annotation = args["use-annotation"]
    categorical_keys = (
        args["categorical"].split(",") if args["categorical"] is not None and not train_reference else None
    )
    continuous_keys = (
        args["continuous"].split(",") if args["continuous"] is not None else None
    )
    n_latent = int(args["n-latent"])
    batch_size = int(args["batch-size"]) if args["n-latent"] is not None else None
    limit_train_batches = (
        "auto"
        if args["limit-train-batches"] == "auto"
        else float(args["limit-train-batches"])
        if args["limit-train-batches"] is not None
        else None
    )
    make_embed = args["make-embed"]
    n_neighbors = int(args["n-neighbors"])

    ad = read_input(input_h5ad)

    if args["train"]:
        train(ad, output_prefix, )
    elif args["query"]:
        query()
    elif args["integrate"]:
        integrate()

def read_input():
    pass


def train(
        ad,
        output_prefix,
        batch_key,
        min_batch_size=0,
        n_hvg=5000,
        remove_genes=None,
        train_reference=False,
        use_annotation=None,
        categorical_keys=None,
        continuous_keys=None,
        make_embed=False,
        ):
    if min_batch_size > 1:
        batch_sizes = ad.obs[batch_key].value_counts()
        batch_to_keep = batch_sizes.index[batch_sizes >= min_batch_size].astype(str).values
        ad = ad[ad.obs[batch_key].isin(batch_to_keep)].copy()

    # Remove certain gene sets if required
    k_remove = np.zeros(ad.n_vars).astype(bool)
    if remove_genes:
        for g in remove_genes:
            if g in ad.var.columns:
                k_remove = k_remove | ad.var[g].values
            elif os.path.exists(g):
                genes = sk.read_list(g)
                k_remove = k_remove | ad.var_names.isin(genes)

    # Subset to HVG
    sc.pp.highly_variable_genes(
        ad,
        flavor="seurat_v3",
        n_top_genes=n_hvg,
        subset=False,
    )
    min_variance = np.sort(-ad.var.variances_norm[~k_remove])[n_hvg]
    k_keep = (ad.var.variances_norm >= min_variance) & ~k_remove
    logging.info("hvg done")

    ad1 = ad[:, ~k_keep].copy()

    # Run scvi
    scvi.model.SCVI.setup_anndata(
        ad1,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_keys,
        continuous_covariate_keys=continuous_keys,
    )
    if limit_train_batches == "auto":
        if ad1.n_obs / batch_size <= 50:
            limit_train_batches = 1.0
        elif ad1.n_obs / batch_size >= 100:
            limit_train_batches = 0.5
        else:
            limit_train_batches = 50
    extra_params = ARCHES_PARAMS if train_reference else {}
    scvi_vae = scvi.model.SCVI(
        ad1, n_layers=2, dropout_rate=0.2, n_latent=n_latent, **extra_params
    )
    scvi_vae.train(
        train_size=0.9,
        early_stopping_patience=30,
        max_epochs=400,
        batch_size=batch_size,
        limit_train_batches=limit_train_batches,
        use_gpu=True,
    )
    scvi_vae.save(output_prefix + ".scvi_vae", overwrite=True)
    ad.obsm["X_scvi"] = scvi_vae.get_latent_representation()
    logging.info("scvi done")

    del ad1

    if use_annotation and use_annotation in ad.obs.keys():
        scan_vae = scvi.model.SCANVI.from_scvi_model(
            scvi_vae,
            unlabeled_category="Unknown",
            labels_key=use_annotation,
        )
        scan_vae.train(
            max_epochs=20,
            n_samples_per_label=100,
            batch_size=batch_size,
            use_gpu=True,
        )
        scan_vae.save(output_prefix + ".scanvi_vae", overwrite=True)
        ad.obsm["X_scanvi"] = scan_vae.get_latent_representation()
        logging.info("scvi done")

    if make_embed:
        ad.obsm["X_mde_scvi"] = scvi.model.utils.mde(ad.obsm["X_scvi"])
        if "X_scanvi" in ad.obsm.keys():
            ad.obsm["X_mde_scanvi"] = scvi.model.utils.mde(ad.obsm["X_scanvi"])
        logging.info("mde done")
        ad.write(output_prefix + ".h5ad", compression="gzip")
    else:
        X_scvi = pd.DataFrame(
            ad.obsm["X_scvi"],
            index=ad.obs_names,
            columns=[f"LV{(i+1):02d}" for i in range(n_latent)],
        )
        X_scvi.to_csv(output_prefix + ".X_scvi.csv.gz")
        if "X_scanvi" in ad.obsm.keys():
            X_scanvi = pd.DataFrame(
                ad.obsm["X_scanvi"],
                index=ad.obs_names,
                columns=[f"LV{(i+1):02d}" for i in range(n_latent)],
            )
            X_scanvi.to_csv(output_prefix + ".X_scanvi.csv.gz")

    print("all done")
    return 0


if __name__ == "__main__":
    from docopt import docopt

    args = docopt(__doc__)
    args = {k.lstrip("-<").rstrip(">"): args[k] for k in args}
    try:
        if args.get("debug"):
            logLevel = logging.DEBUG
        else:
            logLevel = logging.INFO
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
