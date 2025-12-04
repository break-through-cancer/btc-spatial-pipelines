
process SQUIDPY_SPATIAL_PLOTS {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:0471909d51c29a5a4cb391ac86f5cf58dad448da7f6862577d206ae8eb831216'

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("${prefix}/figures/spatial_scatter.png"),     emit: spatial_scatter_plot
    tuple val(meta), path("${prefix}/figures/interaction_matrix.png"),  emit: interaction_matrix_plot
    tuple val(meta), path("${prefix}/figures/co_occurrence*.png"),      emit: co_occurrence_plot
    tuple val(meta), path("${prefix}/figures/nhood_enrichment.png"),    emit: nhood_enrichment_plot
    tuple val(meta), path("${prefix}/figures/centrality_scores.png"),   emit: centrality_scores_plot
    tuple val(meta), path("${prefix}/squidpy.h5ad"),                    emit: adata
    path "versions.yml",                                                emit: versions


    script:
    // we may want to run for multiple deconvolution results and same sample
    source = adata.simpleName
    sample = "${meta.id}"
    prefix = task.ext.prefix ?: "${sample}/${source}"
    template 'spatial.py'
}

process SQUIDPY_LIGREC_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:0471909d51c29a5a4cb391ac86f5cf58dad448da7f6862577d206ae8eb831216'

    input:
    tuple val(meta), path(adata)
    output:
    tuple val(meta), path("${prefix}/**.png"),                            emit: ligrec_plots, optional: true
    tuple val(meta), path("${prefix}/**ligrec_interactions.pickle"),      emit: ligrec_interactions
    tuple val(meta), path("${prefix}/**ligrec_metadata.csv"),             emit: ligrec_metadata
    path "versions.yml",                                                  emit: versions

    script:
    // we may want to run for multiple deconvolution results and same sample
    source = adata.simpleName
    sample = "${meta.id}"
    prefix = task.ext.prefix ?: "${sample}/${source}"
    template 'ligrec.py'
}
