from spatialdata_io import visium_hd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq
import scanpy as sc
import geopandas as gpd
import h5py
import numpy as np
from scipy.sparse import csr_matrix
from pathlib import Path
import random
import json
from anndata import AnnData
from geopandas import GeoDataFrame
from typing import Any
from shapely.geometry import Polygon


#helper to check geometries extraction by spatialdata-io
def _extract_geometries_from_geojson(
    adata: AnnData,
    geojson_path: Path,
) -> GeoDataFrame:
    """Extract geometries and create a GeoDataFrame from a GeoJSON features map.

    Parameters
    ----------
    adata
        AnnData object containing cell data.
    geojson_path
        Path to the GeoJSON file containing cell segmentation geometries.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing cell IDs and their corresponding geometries.
    """
    with open(geojson_path) as f:
        geojson_data = json.load(f)
    geojson_features_map: dict[str, Any] = {
        f"cellid_{feature['properties']['cell_id']:09d}-1": feature for feature in geojson_data["features"]
    }

    geometries = []
    cell_ids_ordered = []

    for obs_index_str in adata.obs.index:
        feature = geojson_features_map.get(obs_index_str)
        print(obs_index_str, feature is not None)
        if feature:
            print(f"Extracting geometry for {obs_index_str}")
            polygon_coords = np.array(feature["geometry"]["coordinates"][0])
            geometries.append(Polygon(polygon_coords))
            cell_ids_ordered.append(obs_index_str)

    return GeoDataFrame({"cell_id": cell_ids_ordered, "geometry": geometries}, index=cell_ids_ordered)

#helper to create a 10X h5 file
def write_10X_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    
    Adapted from https://github.com/scverse/anndata/issues/595#issuecomment-1824376236
    
    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """

    with h5py.File(file, 'w') as w:
        grp = w.create_group("matrix")
        grp.create_dataset("barcodes", data=np.array(adata.obs_names))
        grp.create_dataset("data", data=np.array(adata.X.data), compression="gzip", compression_opts=9)
        ftrs = grp.create_group("features")
        # this group will lack the following keys:
        # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
        ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types))
        ftrs.create_dataset("genome", data=np.array(adata.var.genome))
        ftrs.create_dataset("id", data=np.array(adata.var.gene_ids))
        ftrs.create_dataset("name", data=np.array(adata.var.index))
        grp.create_dataset("indices", data=np.array(adata.X.indices))
        grp.create_dataset("indptr", data=np.array(adata.X.indptr))
        grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1]))


#load test data
sdata = visium_hd('hdvisium', dataset_id='sample')
adata = to_legacy_anndata(sdata, coordinate_system='sample', include_images=True, table_name='square_002um')
print(f'got andata with {adata.n_obs} obs and {adata.n_vars} vars')
print(f'mean counts per spot: {adata.X.mean():.5f}')

#get neighbors info
sq.gr.spatial_neighbors(adata, n_rings=1)
neighbors = adata.obsp['spatial_connectivities'].toarray()

#compute synthetic cells from neighbors and their neighbors
#pick some random cells and grow cells around them
random.seed(42)
cells = random.sample(range(neighbors.shape[0]), 1500)
adata.obs['cell_id'] = None
for i, cell in enumerate(cells):
    #get neighbors of cell
    nbrs = np.where(neighbors[cell,:] > 0)[0]
    #get neighbors of neighbors
    nbrs2 = []
    for nbr in nbrs:
        nbrs2.extend(np.where(neighbors[nbr,:] > 0)[0])
    nbrs2 = list(set(nbrs2))
    #assign cell_id to cell and its neighbors and neighbors of neighbors
    for spot in [cell] + list(nbrs) + nbrs2:
        adata.obs.at[adata.obs_names[spot], 'cell_id'] = i

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

#create polygons for each cell_id from coordinates of spots in that cell in adata.obsm['spatial']
barcode_mappings['y'] = adata.obsm['spatial'][:,1]
barcode_mappings['x'] = adata.obsm['spatial'][:,0]
barcode_mappings['points'] = barcode_mappings.apply(lambda row: gpd.points_from_xy([row['x']], [row['y']])[0], axis=1)
barcode_mappings = barcode_mappings[~barcode_mappings['cell_id'].isnull()]

gdf = gpd.GeoDataFrame(
    barcode_mappings[['cell_id','points']],
    geometry='points'
)


cell_gdf = gdf.dissolve(by='cell_id', as_index=True).convex_hull
#keep only polygon geometry
geo_types = cell_gdf.geometry.type
cell_gdf = cell_gdf[geo_types.isin(['Polygon', 'MultiPolygon', 'Circles'])]


# save the aggregated counts matrix to filtered_feature_cell_matrix.h5, adjust names to 'cellid_NNNNNNNNN-1' format where N is a digit
adata = adata[~adata.obs['cell_id'].isnull(), :]
agg_adata = sc.get.aggregate(adata, by='cell_id', func='sum', axis='obs', layer=None).copy()
agg_adata.obs_names = [f'cellid_{int(cid):09d}-1' for cid in agg_adata.obs_names]
agg_adata.X = csr_matrix(agg_adata.layers['sum'].copy())
del agg_adata.layers['sum']

cell_gdf.to_file('cell_segmentations.geojson')
test = _extract_geometries_from_geojson(agg_adata, Path('cell_segmentations.geojson'))
assert len(test) > 0
write_10X_h5(agg_adata, 'filtered_feature_cell_matrix.h5')