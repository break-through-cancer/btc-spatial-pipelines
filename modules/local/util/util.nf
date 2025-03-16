process MATCH_ADATAS {
    //return adata_sc with gene index matching adata_st by gene name or gene id
    tag "$meta.id"
    label "match_adatas"
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