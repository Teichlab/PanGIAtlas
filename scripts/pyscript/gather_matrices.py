#!/usr/bin/env python
"""Gather starsolo and cellbender outputs into a single h5ad

Usage: gather_matrices.py [options] <metadata> <sample_id>

Options:
  --debug            print debug information
  --profile          print profile information
  --dry              print starsolo and cellbender outputs to gather
  --force            force overwrite if output exists
  --rerun-cb <tsv>   rerun cellbender summary table
  <metadata>         metadata csv
  <sample_id>        sample_id
"""


import logging
import signal
import sys
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os.path
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import sctk as sk


MAPPED_DATA_ROOT = "/lustre/scratch117/cellgen/cellgeni/TIC-starsolo"
CB_RERUN_ROOT = "/lustre/scratch117/cellgen/team205/nh3/20220125_digestive_tract_integration/v2/h5ad/cellbender"


def gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5):
    cr_gene_filtered_ad = sc.read_10x_mtx(cr_gene_filtered_mtx)
    logging.info("cr_gene_filtered_mtx done")
    cr_velo_filtered_ad = sk.read_velocyto(cr_velo_filtered_mtx)
    logging.info("cr_velo_filtered_mtx done")
    cb_gene_filtered_ad = sk.read_cellbender(cb_filtered_h5)
    logging.info("cb_filtered_h5 done")

    common_cells = list(
        set(cr_gene_filtered_ad.obs_names.tolist())
        & set(cb_gene_filtered_ad.obs_names.tolist())
    )
    k_cr = cr_gene_filtered_ad.obs_names.isin(common_cells)
    k_cb = cb_gene_filtered_ad.obs_names.isin(common_cells)
    ad = anndata.AnnData(
        X=cb_gene_filtered_ad.X[np.where(k_cb)[0], :],
        obs=cb_gene_filtered_ad.obs[k_cb].copy(),
        var=cb_gene_filtered_ad.var.copy(),
        layers={
            "raw": cr_gene_filtered_ad.X[np.where(k_cr)[0], :],
            "spliced": cr_velo_filtered_ad.X[np.where(k_cr)[0]],
            "unspliced": cr_velo_filtered_ad.layers["unspliced"][np.where(k_cr)[0]],
            "ambiguous": cr_velo_filtered_ad.layers["ambiguous"][np.where(k_cr)[0]],
        }
    )
    for layer in ("spliced", "unspliced", "ambiguous"):
        ad.layers[layer].eliminate_zeros()
    return ad


def main(args):
    logging.debug(args)

    metadata_csv = args["metadata"]
    sid = args["sample_id"]
    rerun_cb = args["rerun-cb"]

    metadata_df = pd.read_csv(metadata_csv, index_col=0)

    if rerun_cb:
        rerun_cb_df = pd.read_csv(rerun_cb, sep="\t", index_col=0, names=["choice", "quality", "comments", "action"])

    logging.info("metadata read")

    k = np.where(metadata_df.sampleID.values == sid)[0]

    if k.size > 0:
        logging.info(f"{sid} found")
        k = k[0]
    else:
        raise ValueError(f"{sid} not found in {metadata_csv}")

    starsolo_out = metadata_df.starsolo_out.values[k]
    if rerun_cb and sid in rerun_cb_df.index:
        choice = int(rerun_cb_df.loc[sid, "choice"])
        if choice == 2:
            cb_filtered_h5 = f"{CB_RERUN_ROOT}/{sid}.cellbender.out/cellbender_out_filtered.h5"
        elif choice == 3:
            cb_filtered_h5 = f"{CB_RERUN_ROOT}/{sid}.cellbender.out3/cellbender_out_filtered.h5"
        else:
            cb_filtered_h5 = f"{MAPPED_DATA_ROOT}/{metadata_df.cellbender_out.values[k]}"
    else:
        cb_filtered_h5 = f"{MAPPED_DATA_ROOT}/{metadata_df.cellbender_out.values[k]}"

    starsolo_out_path = f"{MAPPED_DATA_ROOT}/{starsolo_out}"
    cr_gene_filtered_mtx = f"{starsolo_out_path}/Gene/filtered"
    cr_velo_filtered_mtx = f"{starsolo_out_path}/Velocyto/filtered"

    output_h5ad = f"h5ad/input/{sid}.gene_velo_cellbender.filtered.h5ad"

    if not os.path.exists(cr_gene_filtered_mtx):
        raise FileNotFoundError(cr_gene_filtered_mtx)
    if not os.path.exists(cr_velo_filtered_mtx):
        raise FileNotFoundError(cr_velo_filtered_mtx)
    if not os.path.exists(cb_filtered_h5):
        raise FileNotFoundError(cb_filtered_h5)

    if args["dry"]:
        print(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5)
    else:
        if os.path.exists(output_h5ad) and not args["force"]:
            logging.info(f"{output_h5ad} exists, skip gathering without --force")
        else:
            ad = gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5)
            ad.write(output_h5ad, compression="gzip")
            logging.info("done")

    return 0


if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    args = {k.lstrip('-<').rstrip('>'):args[k] for k in args}
    try:
        if args.get('debug'):
            logLevel = logging.DEBUG
        else:
            logLevel = logging.INFO
        logging.basicConfig(
            level=logLevel,
            format='%(asctime)s; %(levelname)s; %(funcName)s; %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')
        if args.get('profile'):
            import cProfile
            cProfile.run('main(args)')
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning('Interrupted')
        sys.exit(1)
