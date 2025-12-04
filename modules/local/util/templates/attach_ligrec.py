#!/usr/bin/env python
import pickle
import os
import anndata as ad
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

adata_path = "${adata}"
ligrec_path = "${ligrec}"
sample = "${meta.id}"

#output directory
os.makedirs(sample, exist_ok=True)

log.info(f"attaching ligand-receptor results to adata.uns")
adata = ad.read_h5ad(adata_path)
with open(ligrec_path, "rb") as f:
    ligrec = pickle.load(f)

log.info(f"got ligrec with keys: {list(ligrec.keys())}")

# flatten indices to store in adata
# https://github.com/scverse/squidpy/issues/533
for k, v in list(ligrec.items()):
    v.index = v.index.to_flat_index().str.join("-")
    if k != "metadata":
        v.columns = v.columns.to_flat_index().str.join("-")
adata.uns['ligrec_means'] = ligrec['means'].sparse.to_dense()
adata.uns['ligrec_pvalues'] = ligrec['pvalues'].sparse.to_dense()

log.info(f"saving adata with attached ligand-receptor results to {sample}/staple.h5ad")
adata.write_h5ad(f"{sample}/staple.h5ad", compression='gzip')

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(ad.__version__))
