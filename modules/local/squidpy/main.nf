process SQUIDPY_MORANS_I {
    tag "$meta.id"
    label "process_low"
    container "ghcr.io/break-through-cancer/btc-containers/scverse:main"

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("${prefix}/spatially_variable_genes.csv"), emit: svgs
    path "versions.yml",                                             emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'moran.py'
}


process SQUIDPY_SPATIAL_PLOTS {
    tag "$meta.id"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse:main'

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("${prefix}/figures/spatial_scatter.png"),     emit: spatial_scatter_plot
    tuple val(meta), path("${prefix}/figures/interaction_matrix.png"),  emit: interaction_matrix_plot
    tuple val(meta), path("${prefix}/figures/co_occurrence*.png"),      emit: co_occurrence_plot
    tuple val(meta), path("${prefix}/figures/nhood_enrichment.png"),    emit: nhood_enrichment_plot
    tuple val(meta), path("${prefix}/figures/centrality_scores.png"),   emit: centrality_scores_plot

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
    tuple val(meta), path("${prefix}/ligrec_interactions.csv"),          emit: ligrec_csv
    path "versions.yml",                                                 emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'ligrec.py'
}
