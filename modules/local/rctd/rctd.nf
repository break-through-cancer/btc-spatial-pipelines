process RCTD {
    tag "$meta.id"
    label 'process_high'
    container 'ghcr.io/break-through-cancer/btc-containers/rctd:main'

    input:
    tuple val(meta), path(adata_sc), path(adata_st), val(n_top_genes)
    output:
    tuple val(meta), path("${prefix}/rctd_cell_types.csv"), emit: rctd_cell_types
    tuple val(meta), path("${prefix}/rctd.h5ad"),           emit: rctd_adata
    path "versions.yml",                                    emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'rctd.r'
}

process RCTD_PLOTS {
    tag "$meta.id"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse:main'

    input:
    tuple val(meta), path(rctd_adata)
    output:
    tuple val(meta), path("figures/${prefix}_spatial_scatter.png"),    emit: spatial_scatter_plot
    tuple val(meta), path("figures/${prefix}_interaction_matrix.png"), emit: interaction_matrix_plot

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'plot.py'
}