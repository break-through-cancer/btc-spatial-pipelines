include { SPACEMARKERS } from './analyze_spacemarkers.nf'
include { SQUIDPY } from './analyze_squidpy_ligrec.nf'
include { SQUIDPY_SPATIAL_PLOTS } from '../../../modules/local/squidpy/main'


workflow ANALYZE {

    take: 
        ch_sm_inputs   // from DECONVOLVE
        ch_squidpy
    main:

    versions = channel.empty()
    squidpy_ligrec = channel.empty()
    spacemarkers_ligrec = channel.empty()

    // ligrec - spacemarkers if requested
    if (params.analyze.spacemarkers){
        SPACEMARKERS(ch_sm_inputs)
        versions = versions.mix(SPACEMARKERS.out.versions)
        spacemarkers_ligrec = spacemarkers_ligrec.mix(SPACEMARKERS.out.spacemarkers)
    }

    // ligrec - squidpy if requested
    if (params.analyze.squidpy){
        SQUIDPY( ch_squidpy )
        versions = versions.mix(SQUIDPY.out.versions)
        squidpy_ligrec = squidpy_ligrec.mix(SQUIDPY.out.ligrec)
    }

    // do basic analysis anyway
    SQUIDPY_SPATIAL_PLOTS( ch_squidpy )
    versions = versions.mix(SQUIDPY_SPATIAL_PLOTS.out.versions)
    adata = SQUIDPY_SPATIAL_PLOTS.out.adata

    // wrap up - collect results from tools and save


    emit:
        versions
        adata
}