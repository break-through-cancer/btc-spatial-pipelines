#!/usr/bin/env python3
import anndata as ad
import squidpy as sq
import logging
import scanpy as sc
import os

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


adata_path = "${adata}"
sample = "${prefix}"

os.makedirs(sample, exist_ok=True)

#squidpy insists on dir naming, not creating outdir as usually
log.info("loading {}".format(adata_path))
adata = ad.read_h5ad(adata_path)
log.info("adata is {}".format(adata))

#filter non-na cell_types
adata = adata[~adata.obs["cell_type"].isna()]

#normalize and log1p
log.info("normalizing and log1p transforming data")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

res = sq.gr.ligrec(
    adata,
    n_perms=1000,
    cluster_key="cell_type",
    copy=True,
    use_raw=False,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},
    gene_symbols='index'
)

log.info("saving ligrec mean interation values")
res["means"].to_csv(
    "{sample}/ligrec_means.csv".format(sample),
    index=False
)

log.info("saving ligrec pvalues")
res["pvals"].to_csv(
    "{sample}/ligrec_pvalues.csv".format(sample),
    index=False
)

log.info("saving ligrec metadata")
res["metadata"].to_csv(
    "{sample}/ligrec_metadata.csv".format(sample),
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

with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    scanpy: {}\\n".format(sc.__version__))