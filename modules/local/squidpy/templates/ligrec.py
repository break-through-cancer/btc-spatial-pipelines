#!/usr/bin/env python3
import anndata as ad
import squidpy as sq
import logging
import os

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


adata_path = "${adata}"
sample = "${prefix}"

#squidpy insists on dir naming, not creating outdir as usually
log.info("loading {}".format(adata_path))
adata = ad.read_h5ad(adata_path)
log.info("adata is {}".format(adata))

#filter non-na cell_types
adata = adata[~adata.obs["cell_type"].isna()]

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
             alpha=0.005)
