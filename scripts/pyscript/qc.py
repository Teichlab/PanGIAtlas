#!/usr/bin/env python
"""Perform automatic QC

Usage: qc.py [options] <sample_id>

Options:
  --debug             print debug information
  --profile           print profile information
  --plot_only         plot only
  --qc_metrics <str>  comma-separated list of QC metrics [default: log1p_n_counts,log1p_n_genes,percent_mito,percent_ribo,percent_hb,percent_top50,percent_soup,percent_spliced]
  --models <str>      comma-separated <name>:<model.pkl> pairs giving celltypist models to use [default: ctp_pred:Immune_All_High.pkl]
  --clst_res <float>  resolution for QC clustering [default: 0.2]
  --min_frac <float>  min frac of pass_auto_filter for a cluster to be called good [default: 0.5]
  <sample_id>         sample id
"""


import logging
import signal
import sys
import os

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

import celltypist
import anndata
import scanpy as sc
import sctk as sk


INPUT_H5AD_DIR = "h5ad/input"
POSTQC_H5AD_DIR = "h5ad/post_qc"
QC_PLOT_DIR = "qc_plot"


celltypist_models = {
    "cecilia22_predH": "Immune_All_High.pkl",
    "cecilia22_predL": "Immune_All_Low.pkl",
    "elmentaite21_pred": "Cells_Intestinal_Tract.pkl",
    "suo22_pred": "Pan_Fetal_Human.pkl",
}


def calculate_qc(ad, run_scrublet=True):
    """Calculate QC metrics assuming:
    1) .X contains post-soup-removal counts
    2) .raw contains pre-soup-removal counts
    3) .spliced contains spliced counts
    4) .unspliced contains unspliced counts
    """
    sk.calculate_qc(ad, log1p=False)
    sk.calculate_qc(ad, suffix="_raw", layer="raw", log1p=False)
    sk.calculate_qc(ad, suffix="_spliced", layer="spliced", log1p=False)
    sk.calculate_qc(ad, suffix="_unspliced", layer="unspliced", log1p=False)
    ad.obs["percent_soup"] = (1 - ad.obs["n_counts"] / ad.obs["n_counts_raw"]) * 100
    ad.obs["percent_spliced"] = (
        ad.obs["n_counts_spliced"]
        / (ad.obs["n_counts_spliced"] + ad.obs["n_counts_unspliced"])
        * 100
    )
    if run_scrublet:
        sk.run_scrublet(ad)


def run_QC(
    ad,
    qc_metrics=None,
    models=None,
    res=0.2,
    threshold=0.5,
    relabel_only=False,
    plot_only=False,
):
    if qc_metrics is None:
        qc_metrics = [
            "log1p_n_counts",
            "log1p_n_genes",
            "percent_mito",
            "percent_ribo",
            "percent_hb",
            "percent_top50",
            "percent_soup",
            "percent_spliced",
        ]

    metric_pairs = [("log1p_n_counts", qm) for qm in qc_metrics[1:]]
    if "percent_mito" in qc_metrics and "percent_ribo" in qc_metrics:
        metric_pairs.append(("log(percent_mito)", "percent_ribo"))
    if "percent_soup" in qc_metrics and "percent_spliced" in qc_metrics:
        metric_pairs.append(("log(percent_soup)", "percent_spliced"))
    if "percent_spliced" in qc_metrics and "scrublet_score" in qc_metrics:
        metric_pairs.append(("percent_spliced", "scrublet_score"))

    if models is None:
        models = celltypist_models
    ctp_name = list(models.keys())[0]

    if not plot_only and not relabel_only:
        for name, mod in models.items():
            run_celltypist(ad, mod, min_prob=0.3, key_added=name)

        calculate_qc(ad, run_scrublet=("scrublet_score" in qc_metrics))
        sk._pipeline.generate_qc_clusters(ad, metrics=qc_metrics, res=res)

    mito_thresholds = [20, 50, 80]

    if not plot_only:
        for max_mito in mito_thresholds:
            sk._pipeline.find_good_qc_cluster(
                ad,
                metrics={
                    "n_counts": (1000, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                    "percent_spliced": (50, 97.5, "log", "both", 0.1),
                },
                key_added=f"good_qc_cluster_mito{max_mito}",
            )
            ad.obs.rename(
                columns={
                    f"good_qc_cluster_mito{max_mito}_fqo2": f"pass_auto_filter_mito{max_mito}"
                },
                inplace=True,
            )
        ad.obs["pass_auto_filter"] = 100
        for max_mito in mito_thresholds[::-1]:
            ad.obs.loc[
                ad.obs[f"pass_auto_filter_mito{max_mito}"], "pass_auto_filter"
            ] = max_mito
        # ad.obs["pass_auto_filter"] = ad.obs["pass_auto_filter"].astype("category")
        ad.obs["good_qc_cluster"] = 100
        for max_mito in mito_thresholds[::-1]:
            ad.obs.loc[
                ad.obs[f"good_qc_cluster_mito{max_mito}"], "good_qc_cluster"
            ] = max_mito
        # ad.obs["good_qc_cluster"] = ad.obs["good_qc_cluster"].astype("category")

    ad.uns["qc_cluster_colors"] = sk._plot.make_palette(
        ad.obs["qc_cluster"].cat.categories.size
    )
    sk.set_figsize((3, 3))
    metric_ufig = sc.pl.embedding(
        ad,
        basis="umap_qc",
        color=qc_metrics + ["good_qc_cluster", "qc_cluster"],
        ncols=int(np.ceil((len(qc_metrics) + 2) / 2)),
        show=False,
        vmax=[
            80
            if qm in ("percent_mito", "percent_hb", "percent_soup")
            else 0.5
            if qm == "scrublet_score"
            else None
            for qm in qc_metrics
        ]
        + [None, None],
        size=10,
        return_fig=True,
    )

    metric_vfig = sk.plot_qc_violin(
        ad, groupby="qc_cluster", metrics=qc_metrics, figsize=(3, 3), return_fig=True
    )

    good_cluster_sfig = sk.plot_qc_scatter(
        ad,
        color_by="good_qc_cluster",
        metric_pairs=metric_pairs,
        use_hexbin=False,
        figsize=(3.5, 3),
        wspace=0.3,
        vmin=0,
        vmax=100,
        s=10 ** (2 - np.round(np.log10(ad.n_obs))),
        return_fig=True,
    )
    good_cluster_sfig.suptitle("good qc cluster")

    ad.obs["pass_default"] = (
        (ad.obs.n_counts >= 1e3) & (ad.obs.n_genes >= 200) & (ad.obs.percent_mito < 20)
    )
    ad.obs["pass_default"] = ad.obs["pass_default"].astype("category")

    cell_size = 10 ** (2 - np.round(np.log10(ad.n_obs)))

    pass_default_sfig = sk.plot_qc_scatter(
        ad,
        metric_pairs=metric_pairs,
        color_by="pass_default",
        s=cell_size,
        figsize=(3.5, 3),
        return_fig=True,
    )
    pass_default_sfig.suptitle("pass default filter")

    pass_auto_sfig = sk.plot_qc_scatter(
        ad,
        metric_pairs=metric_pairs,
        color_by="pass_auto_filter",
        s=cell_size,
        figsize=(3.5, 3),
        vmin=0,
        vmax=100,
        return_fig=True,
    )
    pass_auto_sfig.suptitle("pass auto filter")

    ctp_prob_sfig = sk.plot_qc_scatter(
        ad,
        metric_pairs=metric_pairs,
        color_by=f"{ctp_name}_prob",
        s=cell_size,
        figsize=(3.5, 3),
        return_fig=True,
    )
    ctp_prob_sfig.suptitle(f"{ctp_name} prob")

    sk.set_figsize((3, 3))
    sc.pl.embedding(
        ad,
        basis="umap_qc",
        color=[f"{ctp_name}_prob", f"{ctp_name}_uncertain", ctp_name],
        size=30,
        legend_fontsize=12,
        wspace=0.4,
        show=False,
    )
    ctp_pred_ufig = plt.gcf()

    sc.pl.embedding(
        ad,
        basis="umap_qc",
        color=["good_qc_cluster", "pass_auto_filter", "pass_default", "qc_cluster"],
        size=30,
        legend_fontsize=12,
        show=False,
    )
    qc_cluster_ufig = plt.gcf()

    return (
        metric_vfig,
        metric_ufig,
        good_cluster_sfig,
        pass_default_sfig,
        pass_auto_sfig,
        ctp_prob_sfig,
        ctp_pred_ufig,
        qc_cluster_ufig,
    )


def run_celltypist(ad, model, min_prob=0.3, key_added="ctp_pred"):
    aux_ad = ad.copy()
    sc.pp.normalize_total(aux_ad, target_sum=1e4)
    sc.pp.log1p(aux_ad)
    try:
        pred = celltypist.annotate(aux_ad, model=model, majority_voting=True)
        ad.obs[key_added] = pred.predicted_labels.majority_voting
    except ValueError:
        pred = celltypist.annotate(aux_ad, model=model, majority_voting=False)
        ad.obs[key_added] = pred.predicted_labels.predicted_labels
    ad.obs[f"{key_added}_prob"] = pred.probability_matrix.max(axis=1)
    ad.obs[f"{key_added}_uncertain"] = ad.obs[key_added].astype(str)
    ad.obs.loc[
        ad.obs[f"{key_added}_prob"] < min_prob, f"{key_added}_uncertain"
    ] = "Uncertain"
    ad.obs[f"{key_added}_uncertain"] = ad.obs[f"{key_added}_uncertain"].astype(
        "category"
    )
    del aux_ad


def parse_model_option(model_str):
    models = {}
    for model in model_str.split(","):
        if ":" in model:
            name, mod = model.split(":")
        else:
            name = "ctp_pred"
            mod = model
        models[name] = mod

    return models


def process_sample(ad, qc_metrics, models, clst_res, min_frac, plot_only=False):

    qc_figs = run_QC(
        ad,
        qc_metrics=qc_metrics,
        models=models,
        res=clst_res,
        threshold=min_frac,
        plot_only=plot_only,
    )
    return qc_figs


def main(args):
    logging.debug(args)

    sid = args["sample_id"]
    qc_metrics = args["qc_metrics"].split(",")
    models = parse_model_option(args["models"])
    clst_res = float(args["clst_res"])
    min_frac = float(args["min_frac"])
    plot_only = args["plot_only"]

    input_h5ad = os.path.join(
        INPUT_H5AD_DIR, f"{sid}.gene_velo_cellbender.filtered.h5ad"
    )
    postqc_h5ad = os.path.join(
        POSTQC_H5AD_DIR, f"{sid}.gene_velo_cellbender.post_qc.h5ad"
    )
    if plot_only:
        input_h5ad = postqc_h5ad

    ad = sc.read(input_h5ad)

    (
        metric_vfig,
        metric_ufig,
        good_cluster_sfig,
        pass_default_sfig,
        pass_auto_sfig,
        ctp_prob_sfig,
        ctp_pred_ufig,
        qc_cluster_ufig,
    ) = process_sample(ad, qc_metrics, models, clst_res, min_frac, plot_only)

    sid_qc_plot_dir = os.path.join(QC_PLOT_DIR, sid)
    os.makedirs(sid_qc_plot_dir, exist_ok=True)
    metric_vfig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.metric_vfig.png", bbox_inches="tight"
    )
    metric_ufig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.metric_ufig.png", bbox_inches="tight"
    )
    good_cluster_sfig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.good_sfig.png", bbox_inches="tight"
    )
    pass_default_sfig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.default_sfig.png", bbox_inches="tight"
    )
    pass_auto_sfig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.auto_sfig.png", bbox_inches="tight"
    )
    ctp_prob_sfig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.ctp_sfig.png", bbox_inches="tight"
    )
    ctp_pred_ufig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.ctp_ufig.png", bbox_inches="tight"
    )
    qc_cluster_ufig.savefig(
        f"{sid_qc_plot_dir}/{sid}.qc_plot.cluster_ufig.png", bbox_inches="tight"
    )

    if not plot_only:
        ad.write(postqc_h5ad, compression="gzip")

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
