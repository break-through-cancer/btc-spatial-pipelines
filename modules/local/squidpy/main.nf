process SQUIDPY {
    tag "$meta.id"
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/squidpy:main"

    input:
    tuple val(meta), path(data)
    output:
    tuple val(meta), path("${prefix}/spatially_variable_genes.csv"), emit: svgs
    path "versions.yml",                                             emit: versions

    stub:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3
import os
os.makedirs("${prefix}", exist_ok=True)
os.touch("spatially_variable_genes.csv")
with open ("versions.yml", "w") as f:
    f.write(f"${task.process}:\\n")
    f.write(f"    squidpy: {sq.__version__}\\n")
    f.write(f"    anndata: {ad.__version__}\\n")
    """

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env python3
import os
import anndata as ad
import squidpy as sq

os.makedirs("${prefix}", exist_ok=True)
adata = ad.read_h5ad("$data")
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode="moran")

svgs = adata.uns["moranI"]
svgs = svgs[svgs["pval_norm_fdr_bh"] < 0.1]
svgs = svgs[svgs["I"].notnull()]
svgs.to_csv("${prefix}/spatially_variable_genes.csv")

with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
"""
}