from spatialdata_io import visium, visium_hd
import pandas as pd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq

#load test data
sdata = visium_hd('hdvisium', dataset_id='sample', bin_size=2)
adata = to_legacy_anndata(sdata, coordinate_system='sample')
print(f'got andata with {adata.n_obs} obs and {adata.n_vars} vars')
print(f'mean counts per spot: {adata.X.mean():.5f}')

#get neighbors info
sq.gr.spatial_neighbors(adata, n_rings=1)
neighbors = adata.obsp['spatial_connectivities'].toarray()

#compute synthetic cells from neighbors and therir neighbors
adata.obs['cell_id'] = None
for i in range(neighbors.shape[0]):
    neighbor_idxs = neighbors[i,:].nonzero()[0].tolist()
    neighbor_idxs.extend([j for nidx in neighbor_idxs for j in neighbors[nidx,:].nonzero()[0].tolist()])  #include neighbors of neighbors
    neighbor_idxs.append(i)  #include self
    idx = adata.obs.iloc[neighbor_idxs,:].index
    if adata.obs.loc[idx,'cell_id'].isnull().all():
        adata.obs.loc[idx,'cell_id'] = f'cell_{i}'

print(f'computed {adata.obs["cell_id"].nunique()} synthetic cells')

#create dummy barcode_mappings.parquet from adata.obs with following columns:
columns = ['square_002um', 'square_008um', 'square_016um', 'cell_id', 'in_nucleus','in_cell']
barcode_mappings = adata.obs[adata.obs.columns.intersection(columns)].copy()
barcode_mappings.reset_index(inplace=True)
barcode_mappings.rename(columns={'index':'square_002um'}, inplace=True)
barcode_mappings['in_nucleus'] = barcode_mappings['cell_id'].notnull()
barcode_mappings['in_cell'] = barcode_mappings['cell_id'].notnull()

#save to parquet
barcode_mappings.to_parquet('barcode_mappings.parquet')
