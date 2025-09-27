
include { COGAPS_ADATA2DGC;
          COGAPS as COGAPS_MAIN;
          COGAPS_PREPROCESS } from '../../../modules/local/cogaps/main'

workflow COGAPS {
    take:
        ch_adata // Channel of [meta, adata]

    main:

    def patterns = params.npatterns.split(',').collect { it.toInteger() }
    ch_patterns = Channel.from(patterns)
    versions = Channel.empty()

    //example channel with cparams
    ch_fixed_params = Channel.of([niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1])

    ch_cparams = ch_patterns
      .combine(ch_fixed_params)
      .map { tuple([id:it[0].toString(), npatterns:it[0], niterations:it[1].niterations, sparse:it[1].sparse, distributed:it[1].distributed, nsets:it[1].nsets, nthreads:it[1].nthreads]) }

    // convert adata to dgCMatrix
    COGAPS_ADATA2DGC(ch_adata)
    ch_versions = versions.mix(COGAPS_ADATA2DGC.out.versions)


    // preprocess dgCMatrix
    ch_preprocess = COGAPS_ADATA2DGC.out.dgCMatrix
      .map { tuple(it[0], it[1]) }

    COGAPS_PREPROCESS(ch_preprocess)
    ch_versions = ch_versions.mix(COGAPS_PREPROCESS.out.versions)

    ch_input = COGAPS_PREPROCESS.out.dgCMatrix
      .map { tuple(it[0], it[1]) }
      .combine(ch_cparams)

    COGAPS_MAIN(ch_input)
    ch_versions = ch_versions.mix(COGAPS_MAIN.out.versions)
    ch_deconvolved = COGAPS_MAIN.out.cogapsResult.map { tuple(it[0], it[1]) }


    emit:
        ch_deconvolved
        versions
}