#!/usr/bin/env python3
import anndata as ad
import squidpy as sq
import logging
import scanpy as sc
import os
import pickle

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


adata_path = "${adata}"
sample = "${prefix}"

sq_gr_ligrec_threshold = ${params.sq_gr_ligrec_threshold}
sq_gr_ligrec_alpha = ${params.sq_gr_ligrec_alpha}
sq_gr_ligrec_nperms = ${params.sq_gr_ligrec_nperms}
sq_pl_ligrec_pvalue = ${params.sq_pl_ligrec_pvalue}


log.info("loading {}".format(adata_path))
adata = ad.read_h5ad(adata_path)
log.info("adata is {}".format(adata))

#squidpy insists on figure dir naming, not creating plot outdir as usually
os.makedirs(sample, exist_ok=True)
os.chdir(sample)

#filter non-na cell_types
adata = adata[~adata.obs["cell_type"].isna()]

#normalize and log1p
log.info("normalizing and log1p transforming data")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

def default_ligrec(**kwargs):
    if "gene_symbols" in kwargs:
        gene_symbols = kwargs.pop("gene_symbols")
    else:
        gene_symbols = None
    ligrec=sq.gr.ligrec(
        adata,
        n_perms=sq_gr_ligrec_nperms,
        cluster_key="cell_type",
        copy=True,
        use_raw=False,
        corr_method="fdr_bh",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        alpha=sq_gr_ligrec_alpha,
        gene_symbols=gene_symbols,
        threshold=sq_gr_ligrec_threshold
    )
    return ligrec
    

# try running ligrec
res = None
try:
    res = default_ligrec()
except Exception as e:
    log.error("ligrec failed: {}".format(e))
    try:
        log.error("trying ligrec without gene_symbols")
        res = default_ligrec(gene_symbols="index")
    except Exception as e2:
        raise e2

if res is None:
    log.error("ligrec did not return results, exiting.")
    exit(0)
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
    sq.pl.ligrec(res,
                dendrogram = 'interacting_clusters',
                save="ligrec_interactions_by_clusters_{}.png".format(sample),
                title="{} Ligand-Receptor Interaction".format(sample),
                alpha=0.01,
                pvalue_threshold=sq_pl_ligrec_pvalue,
                remove_nonsig_interactions=True,
                swap_axes=True)
    sq.pl.ligrec(res,
                dendrogram = 'interacting_molecules',
                save="ligrec_interactions_by_molecules_{}.png".format(sample),
                title="{} Ligand-Receptor Interaction".format(sample),
                alpha=0.01,
                pvalue_threshold=sq_pl_ligrec_pvalue,
                remove_nonsig_interactions=True,
                swap_axes=False)

os.chdir("..")
with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    scanpy: {}\\n".format(sc.__version__))

