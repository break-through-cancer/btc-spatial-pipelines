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
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse:main'

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("${prefix}/figures/spatial_scatter.png"),     emit: spatial_scatter_plot
    tuple val(meta), path("${prefix}/figures/interaction_matrix.png"),  emit: interaction_matrix_plot
    tuple val(meta), path("${prefix}/figures/co_occurrence*.png"),      emit: co_occurrence_plot
    tuple val(meta), path("${prefix}/figures/nhood_enrichment.png"),    emit: nhood_enrichment_plot
    tuple val(meta), path("${prefix}/figures/centrality_scores.png"),   emit: centrality_scores_plot
    tuple val(meta), path("${prefix}/squidpy.h5ad"),                    emit: adata


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
    tuple val(meta), path("**.png"),                                               emit: ligrec_plots
    tuple val(meta), path("${prefix}/ligrec_interactions.pickle"),                 emit: ligrec_interactions
    tuple val(meta), path("${prefix}/ligrec_metadata.csv"),                        emit: ligrec_metadata
    path "versions.yml",                                                           emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'ligrec.py'
}
