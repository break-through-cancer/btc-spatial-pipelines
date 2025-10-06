process ATLAS_MATCH {
    //return adata_sc with gene index matching adata_st by gene name or gene id
    tag "$meta.id"
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    //adata_sc adata_sc
    //adata_st adata_st
    input:
    tuple val(meta), path(adata_sc), path(adata_st)
    output:
    tuple val(meta), path("${prefix}/adata_matched.h5ad"), emit: adata_matched
    path "versions.yml",                                   emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3
import os
import anndata as ad
import numpy as np

print("Reading adata_sc in the backed mode")
adata_sc = ad.read_h5ad("$adata_sc", backed='r')
print("adata_sc:")
print(adata_sc)

print("Reading adata_st")
adata_st = ad.read_h5ad("$adata_st")
print("adata_st:")
print(adata_st)

os.makedirs("${prefix}", exist_ok=True)

#look for matching indices
matching_index = adata_sc.var.index.intersection(adata_st.var.index)
print(f"Found {len(matching_index)} matching genes in var.index")

#look for adata_sc.index in var["gene_ids"] of adata_st
if 'gene_ids' in adata_st.var.columns:
    matching_gene_ids = adata_sc.var.index.intersection(adata_st.var["gene_ids"])
    print(f"Found {len(matching_gene_ids)} matching genes in var[gene_ids]")
else:
    matching_gene_ids = []

#look for adata_st.index in adata_sc.var["feature_names"]
if 'feature_name' in adata_sc.var.columns:
    matching_feature_names = adata_st.var.index.intersection(adata_sc.var["feature_name"])
    print(f"Found {len(matching_feature_names)} matching genes in var[feature_name]")
else:
    matching_feature_names = []

#find largest matching case
matching_lengths = [len(x) for x in [matching_index, matching_gene_ids, matching_feature_names]]
which_matching = np.argmax(matching_lengths)

if matching_lengths[which_matching] == 0:
    raise RuntimeError("no matching genes found")

if which_matching == 0:
    print("Matching by index")
    matching = matching_index
    adata_st[:, matching].write_h5ad("${prefix}/adata_matched.h5ad")
    print(f"Saved adata_st with {len(matching)} matching genes")
elif which_matching == 1:
    print("Matching by gene_ids")
    matching = matching_gene_ids
    adata_st.var.reset_index(drop=False, inplace=True)
    adata_st.var.set_index("gene_ids", inplace=True)
    adata_st.var.index = adata_st.var.index.astype('object')
    adata_st[:, matching].write_h5ad("${prefix}/adata_matched.h5ad")
    print(f"Saved adata_st with {len(matching)} matching genes")
elif which_matching == 2:
    print("Matching by feature_name")
    matching = matching_feature_names
    m = {value: key for key, value in zip(adata_sc.var.index, adata_sc.var["feature_name"])}
    adata_st.var["name_matched"] = adata_st.var.index.map(m)
    adata_st.var.dropna(subset=["name_matched"], inplace=True)
    adata_st.var.reset_index(drop=False, inplace=True)
    adata_st.var.set_index("name_matched", inplace=True)
    adata_st.var.index = adata_st.var.index.astype('object')
    adata_st[:, adata_st.var.index].write_h5ad("${prefix}/adata_matched.h5ad")
    print(f"Saved adata_st with {len(matching)} matching genes")
else:
    raise RuntimeError("More cases than expected")

adata_sc.file.close()
adata_st.file.close()

with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    numpy: {}\\n".format(np.__version__))
"""
}

process ATLAS_GET {
    //download an atlas anndata file from a url
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
        val(url)
    output:
        path("*.h5ad"),         emit: atlas
        path("versions.yml"),   emit: versions

    script:
    prefix = task.ext.prefix
"""
#!/usr/bin/env python3
import os
import requests
from urllib.parse import urlparse
import boto3

myurl = "${url}"

if not(myurl.endswith(".h5ad")):
    raise ValueError("URL must end with .h5ad")

parsed_url = urlparse(myurl)
file_key = parsed_url.path.lstrip('/')

if myurl.startswith("s3://"):
    print("Downloading from S3")
    bucket_name = parsed_url.netloc
    s3 = boto3.client('s3')
    s3.download_file(bucket_name, file_key, os.path.basename(file_key))
else:
    print("Downloading from http")
    r = requests.get(myurl)
    r.raise_for_status()
    with open(os.path.basename(file_key), "wb") as f:
        f.write(r.content)
print(f"Downloaded atlas from {myurl}")

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    requests: {}\\n".format(requests.__version__))
    f.write("    boto3: {}\\n".format(boto3.__version__))
"""
}

process QC {
    //generate a simple report of the atlas adata
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
        tuple val(meta), path(adata), val(report_name)
    output:
        path("*report.csv"),                     emit: report
        path("versions.yml"),                    emit: versions

    script:
"""
#!/usr/bin/env python3
import os
import anndata as ad
import pandas as pd
import scanpy as sc


adata_path = "$adata"
outname = "$report_name"
adata = ad.read_h5ad(adata_path)

qc = sc.pp.calculate_qc_metrics(adata, inplace=False)

report = pd.DataFrame({
    "Sample": ["${meta.id}"],
    "n_genes": qc[0].shape[0],
    "n_cells": qc[1].shape[0],
    "mean_genes_by_counts": qc[0]["n_genes_by_counts"].mean(),
    "mean_cells_by_counts": qc[1]["n_cells_by_counts"].mean(),
    "mean_total_nnz_counts": adata.X[adata.X.nonzero()].mean(),
})

report.to_csv(f"{outname}_report.csv", index=False)

adata.file.close()

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    pandas: {}\\n".format(pd.__version__))
"""
}

process ADATA_FROM_VISIUM_HD {
    //convert vhd file to h5ad
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("${prefix}/${params.visium_hd}.h5ad"),   emit: adata
        path("versions.yml"),                                          emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq

sample = "${prefix}"
data = "${data}"
table = "${params.visium_hd}"
os.makedirs(sample, exist_ok=True)

#read visium_hd dataset
ds = sd.visium_hd(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
adata = to_legacy_anndata(ds, coordinate_system=sample,
                          table_name=table, include_images=True)
adata.var_names_make_unique()

#make compatible with BayesTME (uses an older, scanpy notation)
adata.X = adata.X.astype(int)
adata.uns['layout'] = 'IRREGULAR'
sq.gr.spatial_neighbors(adata)
adata.obsp['connectivities'] = adata.obsp['spatial_connectivities'].astype(bool)

#save
outname = os.path.join(sample, f"{table}.h5ad")
adata.write_h5ad(filename=outname)

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))
"""
}

process ADATA_FROM_VISIUM {
    //convert visium dir to h5ad
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("${prefix}/visium.h5ad"),         emit: adata
        path("versions.yml"),                                   emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3

import os
import spatialdata_io as sd
from spatialdata_io.experimental import to_legacy_anndata
import squidpy as sq

sample = "${prefix}"
data = "${data}"

os.makedirs(sample, exist_ok=True)

#read visium dataset
ds = sd.visium(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
adata = to_legacy_anndata(ds, coordinate_system=sample,
                          include_images=True)
adata.var_names_make_unique()

#make compatible with BayesTME (uses an older, scanpy notation)
adata.X = adata.X.astype(int)
adata.uns['layout'] = 'IRREGULAR'
sq.gr.spatial_neighbors(adata)
adata.obsp['connectivities'] = adata.obsp['spatial_connectivities'].astype(bool)

#save
outname = os.path.join(sample, "visium.h5ad")
adata.write_h5ad(filename=outname)

#versions
with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    spatialdata_io: {}\\n".format(sd.__version__))
    f.write("    squidpy: {}\\n".format(sq.__version__))
"""
}
