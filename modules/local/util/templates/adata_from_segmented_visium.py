#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc


sample = "${prefix}"
data = "${data}"
table = "square_002um"
visium_hd = "${params.visium_hd}"

os.makedirs(sample, exist_ok=True)

#read visium_hd dataset
ds = sd.visium_hd(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
adata = to_legacy_anndata(ds, coordinate_system=sample, 
                          table_name = table, include_images=True)
adata.var_names_make_unique()
adata.uns['layout'] = 'IRREGULAR'


#read bin to cell mappings
barcode_mappings = pd.read_parquet(os.path.join(data, "barcode_mappings.parquet"))

#aggregate sparse adata X by summing over cell_id, keep sparse
adata.obs['cell_id'] = barcode_mappings.set_index(table).loc[adata.obs_names]['cell_id'].values
adata = adata[~adata.obs['cell_id'].isnull(), :]

cell_adata = sc.get.aggregate(
    adata,
    by='cell_id',
    func='sum',
    axis='obs',
    layer=None
)


#for each cell, compute cell centers from adata.uns['spatial']
cell_centers_array = []
for cell_id in cell_adata.obs_names:
    if pd.isnull(cell_id):
        continue
    cell_barcodes = barcode_mappings[barcode_mappings['cell_id'] == cell_id][table].tolist()
    cell_indexes = np.where(adata.obs_names.isin(cell_barcodes))[0]
    cell_centers = adata.obsm['spatial'][cell_indexes,:].mean(axis=0).tolist()
    cell_centers_array.append(cell_centers)

#compose a new adata with cell-level gene data
adata = ad.AnnData(
    X=cell_adata.layers['sum'],
    obs=cell_adata.obs.copy(),
    var=adata.var.copy(),
    uns=adata.uns.copy(),
    obsm={'spatial': np.array(cell_centers_array)}
)

#save
outname = os.path.join(sample, f"{mod_name}.h5ad")
adata.write_h5ad(filename=outname)

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))