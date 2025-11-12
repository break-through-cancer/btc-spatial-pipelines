#!/usr/bin/env python3
import os
import anndata as ad
import squidpy as sq

os.makedirs("${prefix}", exist_ok=True)
adata = ad.read_h5ad("$adata")
process = "${task.process}"
seed = ${params.seed}
nperms = ${params.sq_gr_spatial_autocorr_nperms}
n_jobs = ${task.cpus}

if nperms <= 0:
    nperms = None
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode="moran", seed=seed, n_perms=nperms, n_jobs=n_jobs)

svgs = adata.uns["moranI"]
svgs = svgs[svgs["pval_norm"] < 0.05]
svgs = svgs[svgs["I"].notnull()]
svgs.to_csv("${prefix}/spatially_variable_genes.csv")

#versions
with open ("versions.yml", "w") as f:
    f.write("{}:\\n".format(process))
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
