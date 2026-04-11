process RCTD {
    tag "$meta.id"
    label 'process_high'
    container 'ghcr.io/break-through-cancer/btc-containers/rctd@sha256:145facedacf5198ff97e7fdca314e24e3890d290b3fb9015ce44b710f5ecd95e'

    input:
    tuple val(meta), path(adata_sc), path(adata_st)
    output:
    tuple val(meta), path("${prefix}/rctd_cell_types.csv"), emit: rctd_cell_types
    tuple val(meta), path("${prefix}/rctd.rds"),            emit: rctd_object
    path "versions.yml",                                    emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'rctd.r'
}
