process ATLAS_MATCH {
    //return adata_sc with gene index matching adata_st by gene name or gene id
    tag "$meta.id"
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/squidpy:main"

    input:
    tuple val(meta), path(adata_sc), path(adata_st)
    output:
    tuple val(meta), path("${prefix}/adata_sc_matched.h5ad"), emit: adata_sc_matched
    path "versions.yml",                                      emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3
import os
import anndata as ad

adata_st = ad.read_h5ad("$adata_st")
adata_sc = ad.read_h5ad("$adata_sc")

os.makedirs("${prefix}", exist_ok=True)

matching = adata_st.var.index.intersection(adata_sc.var.index)
print(f"Found {len(matching)} matching genes in var.index")

if (len(matching) > 0):
    adata_sc = adata_sc[:, matching]
    adata_sc.write_h5ad("${prefix}/adata_sc_matched.h5ad")
    print(f"Saved adata_sc with {len(matching)} matching genes")
else:
    print("Trying to match by feature_name")
    matching = adata_st.var.index.intersection(adata_sc.var["feature_name"])
    if (len(matching) > 0):
        print(f"Found {len(matching)} matching genes in var[feature_name], resetting index.")
        adata_sc.var.set_index("feature_name", inplace=True)
        adata_sc.var.index = adata_sc.var.index.astype('object')
        adata_sc = adata_sc[:, matching]
        adata_sc.write_h5ad("${prefix}/adata_sc_matched.h5ad")
        print(f"Saved adata_sc with {len(matching)} matching genes")
    else:
        print("No matching genes found")

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
        path("atlas.h5ad"),  emit: atlas

    script:
    prefix = task.ext.prefix
"""
#!/usr/bin/env python3
import os
import anndata as ad
import requests
import validators

if not(validators.url("$url")):
    print("Invalid URL")
    exit(1)
if not("$url".endswith(".h5ad")):
    print("URL does not end with .h5ad")
    exit(1)

r = requests.get("$url")
r.raise_for_status()
with open("atlas.h5ad", "wb") as f:
    f.write(r.content)
print("Downloaded atlas from $url")
"""
}