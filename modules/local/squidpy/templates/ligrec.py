#!/usr/bin/env python3
import anndata as ad
import squidpy as sq
import logging
import scanpy as sc
import os
import pickle
from matplotlib import pyplot as plt


def default_ligrec(adata, **kwargs):
    if "gene_symbols" in kwargs:
        gene_symbols = kwargs.pop("gene_symbols")
    else:
        gene_symbols = None
    ligrec=sq.gr.ligrec(
        adata,
        n_perms=par["sq_gr_ligrec_nperms"],
        n_jobs=cpus,
        cluster_key="cell_type",
        copy=True,
        use_raw=False,
        corr_method="fdr_bh",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        alpha=par["sq_gr_ligrec_alpha"],
        gene_symbols=gene_symbols,
        threshold=par["sq_gr_ligrec_threshold"],
        numba_parallel=False)
    
    return ligrec

def default_ligrec_pl(ligrec, **kwargs): 
    source_groups = kwargs.pop("source_groups")
    target_groups = kwargs.pop("target_groups")
    save = kwargs.pop("save")

    sq.pl.ligrec(ligrec,
            source_groups=source_groups,
            target_groups=target_groups,
            save=save,
            title="", # cell names and title overlap sometimes
            alpha=0.01,
            pvalue_threshold=par["sq_pl_ligrec_pvalue"],
            remove_nonsig_interactions=True,
            swap_axes=False)
    plt.close('all')
    return None
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    log = logging.getLogger()


    adata_path = "${adata}"
    sample = "${prefix}"
    process = "${task.process}"
    cpus = ${task.cpus}

    par = {"sq_gr_ligrec_threshold": ${params.sq_gr_ligrec_threshold},
        "sq_gr_ligrec_alpha": ${params.sq_gr_ligrec_alpha},
        "sq_gr_ligrec_nperms": ${params.sq_gr_ligrec_nperms},
        "sq_pl_ligrec_pvalue": ${params.sq_pl_ligrec_pvalue}}

    log.info("received params:{}".format(par))

    log.info("loading {}".format(adata_path))
    adata = ad.read_h5ad(adata_path)
    log.info("adata is {}".format(adata))

    #squidpy insists on figure dir naming, not creating plot outdir as usually
    base_path = os.getcwd()
    save_path = os.path.join(sample, "ligrec")
    os.makedirs(save_path, exist_ok=True)
    os.chdir(save_path)

    #filter non-na cell_types
    adata = adata[~adata.obs["cell_type"].isna()].copy()
    clusters = adata.obs["cell_type"].unique()

    #normalize and log1p
    log.info("normalizing and log1p transforming data")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # try running ligrec
    res = None
    try:
        res = default_ligrec(adata)
    except Exception as e:
        log.error("default ligrec failed: {}".format(e))
        try:
            res = default_ligrec(adata, gene_symbols="index")
        except Exception as e2:
            log.error("ligrec without gene_symbols failed: {}".format(e2))

    if res is None:
        log.error("ligrec did not return results.")
    else:
        log.info("ligres completed successfully, saving to pickle")
        # dictionary of pandas frames: means, pvalues, metadata
        pickle.dump(res, open("ligrec_interactions.pickle", "wb"))

        log.info("saving ligrec metadata")
        # just metadata, as it's easy to comprehend, one row per interaction
        res["metadata"].to_csv(
            "ligrec_metadata.csv"
        )

        log.info("saving ligrec interaction plot")
        for s in clusters:
            try:
                log.info("plotting ligrec for source cluster {}".format(s))
                default_ligrec_pl(ligrec=res, source_groups=s, target_groups=clusters, save="source_{}.png".format(s))
            except Exception as e:
                log.error("ligrec plot for source {} failed: {}".format(s, e))
                continue
        
        for t in clusters:
            try:
                log.info("plotting ligrec for target cluster {}".format(t))
                default_ligrec_pl(ligrec=res, source_groups=clusters, target_groups=t, save="target_{}.png".format(t))
            except Exception as e:
                log.error("ligrec plot for target {} failed: {}".format(t, e))
                continue

    os.chdir(base_path)
    with open ("versions.yml", "w") as f:
        f.write("{}:\\n".format(process))
        f.write("    squidpy: {}\\n".format(sq.__version__))
        f.write("    anndata: {}\\n".format(ad.__version__))
        f.write("    scanpy: {}\\n".format(sc.__version__))