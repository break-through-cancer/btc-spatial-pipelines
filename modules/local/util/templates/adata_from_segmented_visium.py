#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

sample = "${prefix}"
data = "${data}"
table = "square_002um"
visium_hd = "${params.visium_hd}"

os.makedirs(sample, exist_ok=True)

#read visium_hd dataset
log.info(f"loading Visium HD table {table} for sample {sample}")
ds = sd.visium_hd(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
log.info(f"converting to anndata")
adata = to_legacy_anndata(ds, coordinate_system=sample, 
                          table_name = table, include_images=True)
adata.var_names_make_unique()
adata.uns['layout'] = 'IRREGULAR'
log.info(f"adata is {adata}")

#read bin to cell mappings
log.info(f"loading segmentation mappings")
barcode_mappings = pd.read_parquet(os.path.join(data, "barcode_mappings.parquet"))

#aggregate sparse adata X by summing over cell_id, keep sparse
log.info(f"aggregating binned data to cell level")
adata.obs['cell_id'] = barcode_mappings.set_index(table).loc[adata.obs_names]['cell_id'].values
adata = adata[~adata.obs['cell_id'].isnull(), :]

cell_adata = sc.get.aggregate(
    adata,
    by='cell_id',
    func='sum',
    axis='obs',
    layer=None
)

log.info(f"cell_adata is {cell_adata}")

#compute cell centers
log.info(f"computing cell center coordinates")
spatial_df = pd.DataFrame({table: adata.obs_names,
                           0:adata.obsm['spatial'][:,0],
                           1:adata.obsm['spatial'][:,1],
                           'cell_id':adata.obs['cell_id'].values}).set_index(table)
center_coords = spatial_df.groupby('cell_id').mean()
n_bins_per_cell = spatial_df.groupby('cell_id').size()
cell_adata.obs['bins_per_cell'] = n_bins_per_cell.loc[cell_adata.obs_names].values

#estimate new cell diameter based on average bin count per cell
spot_diameter = adata.uns['spatial'][f"{sample}_hires_image"]['scalefactors']['spot_diameter_fullres']
avg_bins_per_cell = n_bins_per_cell.mean()
new_spot_diameter = spot_diameter * 2 * np.sqrt(avg_bins_per_cell/3.14)
adata.uns['spatial'][f"{sample}_hires_image"]['scalefactors']['spot_diameter_fullres'] = new_spot_diameter
adata.uns['spatial'][f"{sample}_hires_image"]['scalefactors']['orig_bin_spot_diameter'] = spot_diameter

log.info(f"writing mean cell diameter {new_spot_diameter:.2f} (ex {spot_diameter:.2f}) to adata.uns")

assert center_coords.shape[0] == cell_adata.n_obs, "number of cells do not match after aggregation"
assert all(center_coords.index == cell_adata.obs_names), "cell ids do not match after aggregation"

#compose a new adata with cell-level gene data
adata = ad.AnnData(
    X=cell_adata.layers['sum'],
    obs=cell_adata.obs.copy(),
    var=adata.var.copy(),
    uns=adata.uns.copy(),
    obsm={'spatial': np.array(center_coords)}
)

#save
outname = os.path.join(sample, f"{visium_hd}.h5ad")
adata.write_h5ad(filename=outname)

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))