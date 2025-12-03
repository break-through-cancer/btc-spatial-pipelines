include { SPACEMARKERS } from './analyze_spacemarkers.nf'
include { SQUIDPY } from './analyze_squidpy_ligrec.nf'
include { SQUIDPY_MORANS_I;
          SQUIDPY_SPATIAL_PLOTS } from '../../../modules/local/squidpy/main'


workflow ANALYZE {

    take: 
        ch_sm_inputs   // from DECONVOLVE
        ch_squidpy
    main:

    versions = channel.empty()
    ligrec = channel.empty()

    // do basic analysis anyway
    SQUIDPY_MORANS_I( ch_squidpy )
    versions = versions.mix(SQUIDPY_MORANS_I.out.versions)

    SQUIDPY_SPATIAL_PLOTS( ch_squidpy )
    versions = versions.mix(SQUIDPY_SPATIAL_PLOTS.out.versions)
    adata = SQUIDPY_SPATIAL_PLOTS.out.adata

    // ligrec - spacemarkers if requested
    if (params.analyze.spacemarkers){
        SPACEMARKERS(ch_sm_inputs)
        versions = versions.mix(SPACEMARKERS.out.versions)
        ligrec = ligrec.mix(SPACEMARKERS.out.spacemarkers)
    }

    // ligrec - squidpy if requested
    if (params.analyze.squidpy){
        SQUIDPY( ch_squidpy )
        versions = versions.mix(SQUIDPY.out.versions)
        ligrec = ligrec.mix(SQUIDPY.out.ligrec)
    }


    emit:
        versions
        ligrec
        adata
}