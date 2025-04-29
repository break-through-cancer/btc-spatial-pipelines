process ATLAS_MATCH {
    //return adata1 with gene index matching adata2 by gene name or gene id
    tag "$meta.id"
    label "process_high_memory"
    container "ghcr.io/break-through-cancer/btc-containers/squidpy:main"

    //adata_sc adata1
    //adata_st adata2
    input:
    tuple val(meta), path(adata1), path(adata2)
    output:
    path("${prefix}/adata_matched.h5ad"),                  emit: adata_matched
    path "versions.yml",                                   emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3
import os
import anndata as ad

print("Reading adata1 in the backed mode")
adata1 = ad.read_h5ad("$adata1", backed='r')
print("adata1:")
print(adata1)

print("Reading adata2")
adata2 = ad.read_h5ad("$adata2")
print("adata2:")
print(adata2)

os.makedirs("${prefix}", exist_ok=True)

matching = adata2.var.index.intersection(adata1.var.index)
print(f"Found {len(matching)} matching genes in var.index")

if (len(matching) > 0):
    adata2[:, matching].write_h5ad("${prefix}/adata_matched.h5ad")
    print(f"Saved adata2 with {len(matching)} matching genes")
else:
    print("Trying to match by feature_name")
    matching = adata2.var_names.intersection(adata1.var["feature_name"])
    if (len(matching) > 0):
        print(f"Found {len(matching)} matching genes in var[feature_name], resetting index.")
        m = {value: key for key, value in zip(adata1.var.index, adata1.var["feature_name"])}
        adata2.var["name_matched"] = adata2.var.index.map(m)
        adata2.var.dropna(subset=["name_matched"], inplace=True)
        adata2.var.reset_index(drop=False, inplace=True)
        adata2.var.set_index("name_matched", inplace=True)
        adata2.var.index = adata2.var.index.astype('object')
        adata2[:, adata2.var.index].write_h5ad("${prefix}/adata_matched.h5ad")
        print(f"Saved adata2 with {len(matching)} matching genes")
    else:
        print("No matching genes found")

adata1.file.close()
adata2.file.close()

with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    anndata: {}\\n".format(ad.__version__))
"""
}

process ATLAS_GET {
    //download an atlas anndata file from a url
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/squidpy:main"

    input:
        val(url)
    output:
        path("*.h5ad"),  emit: atlas

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
    print("URL does not end with .h5ad")
    exit(1)

parsed_url = urlparse(myurl)
file_key = parsed_url.path.lstrip('/')

if (myurl.startswith("s3://")):
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
"""
}

process VHD_TO_H5AD {
    //convert vhd file to h5ad
    label "process_medium"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
        tuple val(meta), path(data)
    output:
        tuple val(meta), path("**.h5ad"),  emit: adata

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3

import os
from spatialdata_io import visium_hd
from spatialdata_io.experimental import to_legacy_anndata
import dask

sample = "${prefix}"
data = "${data}"
table = "${params.hd}"
os.makedirs(sample, exist_ok=True)

#read visium_hd dataset
ds = visium_hd(data, dataset_id=sample, var_names_make_unique=True)

#convert to anndata
adata = to_legacy_anndata(ds, coordinate_system=sample,
                          table_name=table, include_images=True)
#save
outname = os.path.join(sample, f"{table}.h5ad")
adata.write_h5ad(filename=outname)
"""
}