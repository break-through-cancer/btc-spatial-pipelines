import squidpy as sq
import scanpy as sc
import anndata as ad
import pandas as pd

import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata


#visiums SD
sdata = sd.visium("data/visium", dataset_id="testdata")
adata = to_legacy_anndata(sdata, coordinate_system="testdata")


def get_dummy_cell_types(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    ct = pd.DataFrame(adata.obsm["X_pca"][:,:5])
    ct.columns = [f"cell_type_{i}" for i in range(ct.shape[1])]
    ct = abs(ct).apply(lambda x: x/(sum(x)), axis=1)
    ct['barcode'] = adata.obs.index

    return ct

visium_cellular_compositions = get_dummy_cell_types(adata)
visium_cellular_compositions.to_csv("data/visium/visium_cellular_compositions.csv", index=False)

#visium HD
sdata_hd = sd.visium_hd("data/hdvisium", dataset_id="testdata_hd")
adata_hd = to_legacy_anndata(sdata_hd, coordinate_system="testdata_hd")
hdvisium_leden_cellular_compositions = get_dummy_cell_types(adata_hd)
hdvisium_leden_cellular_compositions.to_csv("data/hdvisium/hdvisium_cellular_compositions.csv", index=False)

