import os
import anndata as ad
import sys
import scanpy as sc

adata = ad.read_h5ad("../../out/staple/sample1/staple.h5ad")

# 1. Check main object size (overall memory footprint)
sys.getsizeof(adata)/1e6  # in MB

for key, val in adata.uns.items():
    print(f"{key}:{sys.getsizeof(val)/1e6}")  # in MB

sys.getsizeof(adata.var)/1e6  # in MB
sys.getsizeof(adata.obs)/1e6  # in MB

for key, val in adata.obsp.items():
    print(f"{key}:{sys.getsizeof(val)/1e6}")  # in MB


#get hightly variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)

adata.write_h5ad("../../out/staple/sample1/staple_small.h5ad", compression="gzip", compression_opts=9)