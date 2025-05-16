process SQUIDPY_MORANS_I {
    tag "$meta.id"
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
    tuple val(meta), path(adata)
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
adata = ad.read_h5ad("$adata")
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode="moran")

svgs = adata.uns["moranI"]
svgs = svgs[svgs["pval_norm"] < 0.05]
svgs = svgs[svgs["I"].notnull()]
svgs.to_csv("${prefix}/spatially_variable_genes.csv")

with open ("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
"""
}


process SQUIDPY_SPATIAL_PLOTS {
    tag "$meta.id"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse:main'

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("figures/spatial_scatter_${prefix}.png"),    emit: spatial_scatter_plot
    tuple val(meta), path("figures/interaction_matrix_${prefix}.png"), emit: interaction_matrix_plot
    tuple val(meta), path("figures/co_occurrence_${prefix}.png"),      emit: co_occurrence_plot
    tuple val(meta), path("figures/nhood_enrichment_${prefix}.png"),   emit: nhood_enrichment_plot

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'plot.py'
}

process SQUIDPY_LIGREC_ANALYSIS { //WIP
    tag "$meta.id"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse:main'

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("figures/ligrec_interactions_${prefix}.png"),  emit: ligrec_plot
    tuple val(meta), path("${prefix}/ligrec_interactions.csv"),         emit: ligrec_csv
    path "versions.yml",                                                emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'ligrec.py'
}
