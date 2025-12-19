include { SPACEMARKERS } from './analyze_spacemarkers.nf'
include { SQUIDPY } from './analyze_squidpy_ligrec.nf'
include { SQUIDPY_SPATIAL_PLOTS } from '../../../modules/local/squidpy/main'
include { STAPLE_ATTACH_LIGREC } from '../../../modules/local/util/main'
include { STAPLE_ATTACH_LIGREC as STAPLE_ATTACH_IMSCORES } from '../../../modules/local/util/main'
include { STAPLE_ATTACH_LIGREC as STAPLE_ATTACH_LRSCORES } from '../../../modules/local/util/main'


workflow ANALYZE {

    take: 
        ch_sm_inputs   // from DECONVOLVE
        ch_squidpy
    main:

    versions = channel.empty()
    ligrec = channel.empty()
    imscores = channel.empty()
    lrscores = channel.empty()


    // ligrec - spacemarkers if requested
    if (params.analyze.spacemarkers){
        SPACEMARKERS(ch_sm_inputs)
        versions = versions.mix(SPACEMARKERS.out.versions)
        // pass ligrec results along
        imscores = imscores.mix(SPACEMARKERS.out.imscores)
        lrscores = lrscores.mix(SPACEMARKERS.out.lrscores)
    }

    // ligrec - squidpy if requested
    if (params.analyze.squidpy){
        SQUIDPY( ch_squidpy )
        versions = versions.mix(SQUIDPY.out.versions)
        ligrec = ligrec.mix(SQUIDPY.out.ligrec)
    }

    // do basic analysis anyway
    SQUIDPY_SPATIAL_PLOTS( ch_squidpy )
    versions = versions.mix(SQUIDPY_SPATIAL_PLOTS.out.versions)

    // wrap up - collect results from tools and save
    // TODO: make every subworkflow write ./staple/staple.h5ad
    STAPLE_ATTACH_LIGREC(SQUIDPY_SPATIAL_PLOTS.out.adata.join(ligrec))
    STAPLE_ATTACH_IMSCORES(STAPLE_ATTACH_LIGREC.out.adata.join(imscores))
    STAPLE_ATTACH_LRSCORES(STAPLE_ATTACH_IMSCORES.out.adata.join(lrscores))

    adata = STAPLE_ATTACH_LRSCORES.out.adata

    emit:
        versions
        adata
}