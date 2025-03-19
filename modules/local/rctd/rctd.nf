process RCTD {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/rctd:main'

    input:
    tuple val(meta), path(adata_sc), path(adata_st)
    output:
    tuple val(meta), path("${prefix}/rctd_cell_types.csv"), emit: rtcd_cell_types

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'rctd.r'
}