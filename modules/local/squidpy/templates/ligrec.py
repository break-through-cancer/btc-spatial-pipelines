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

# try running ligrec
res = None
try:
    res = sq.gr.ligrec(
        adata,
        n_perms=1000,
        cluster_key="cell_type",
        copy=True,
        use_raw=False,
        corr_method="fdr_bh",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        gene_symbols='index'
    )
except Exception as e:
    log.error("ligrec failed: {}".format(e))

if res is None:
    log.error("ligrec did not return results, exiting.")
    exit(0)
else:
    log.info("ligres completed successfully, saving to pickle")
    pickle.dump(res, open("ligrec_interactions.pickle", "wb"))
    
    log.info("saving ligrec mean interation values")
    res["means"].to_json(
        "ligrec_means.json",
        index=False
    )

    log.info("saving ligrec pvalues")
    res["pvalues"].to_json(
        "ligrec_pvalues.json",
        index=False
    )

    log.info("saving ligrec metadata")
    res["metadata"].to_json(
        "ligrec_metadata.json",
        index=False
    )

    log.info("saving ligrec interaction plot")
    sq.pl.ligrec(res,
                #source_groups=
                dendrogram = 'interacting_clusters',
                save="ligrec_interactions_{}.png".format(sample),
                title="{} Ligand-Receptor Interaction".format(sample),
                alpha=0.005,
                swap_axes=True)

os.chdir("..")
with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    scanpy: {}\\n".format(sc.__version__))