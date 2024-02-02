#!/usr/bin/env python
"""Program description.

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
  --make-plot                  plot UMAP
  --color-by <str>             column in `adata.obs` that cells are colored by
  --n-neighbors <int>          number of nearest neighbors [default: 15]
  --debug                      print debug information
  --profile                    print profile information
  <input>                      input h5ad
  <out_prefix>                 output prefix
"""


import logging
import signal
import sys
import os

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import scanpy as sc
import sctk as sk
import scvi

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


ARCHES_PARAMS = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
)


def run_mde(
    data,
    n_neighbors=15,
    random_state=0,
    **kwargs,
):
    """
    Util to run :func:`pymde.preserve_neighbors` for visualization of scvi-tools embeddings.

    Parameters
    ----------
    data
        The data of shape (n_obs, k), where k is typically defined by one of the models
        in scvi-tools that produces an embedding (e.g., :class:`~scvi.model.SCVI`.)
    n_neighbors
        Number of nearest neighbors, if set to None, a sensible number will be chosen according
        to the number of observations in `data`.
    random_state
        Random seed passed to mde for reproducibility
    kwargs
        Keyword args to :func:`pymde.preserve_neighbors`
    Returns
    -------
    The pymde embedding, defaults to two dimensions.

    Notes
    -----
    This function is a modification of scvi.model.utils.mde().

    If you use this function in your research please cite:

    Agrawal, Akshay, Alnur Ali, and Stephen Boyd. "Minimum-distortion embedding." arXiv preprint arXiv:2103.02559 (2021).
    """
    try:
        import pymde
        import torch
    except ImportError:
        raise ImportError("Please install pymde package via `pip install pymde`")

    if isinstance(data, pd.DataFrame):
        data = data.values

    device = "cuda"

    _kwargs = dict(
        embedding_dim=2,
        constraint=pymde.Standardized(),
        repulsive_fraction=0.5,
        verbose=False,
        device=device,
        n_neighbors=n_neighbors,
    )
    _kwargs.update(kwargs)

    pymde.seed(random_state)
    mde = pymde.preserve_neighbors(data, **_kwargs)

    return mde


def mde_umap(mde):
    emb = mde.embed(verbose=False)

    emb = emb.cpu().numpy()

    return emb


def mde_neighbors(mde, n):
    edges = mde.edges.cpu().numpy()
    weights = mde.distortion_function.weights

    connectivities = csr_matrix((np.ones(weights.size), (edges[:, 0], edges[:, 1])), shape=(n, n))
    distances = csr_matrix((weights, (edges[:, 0], edges[:, 1])), shape=(n, n))

    return connectivities, distances


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
