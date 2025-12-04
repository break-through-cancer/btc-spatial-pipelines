include { SPACEMARKERS_HD } from '../../../modules/local/spacemarkers/nextflow/main_hd'
include { SPACEMARKERS as SPACEMARKERS_SD } from '../../../modules/local/spacemarkers/nextflow/main'

workflow SPACEMARKERS {

    take:
        ch_sm_inputs   // from DECONVOLVE, csv if available or (BayesTME, CoGAPS) obj - SpaceMarkers knows how to handle both
    main:
        versions = channel.empty()
        if(params.visium_hd) {
            SPACEMARKERS_HD( ch_sm_inputs.map {it -> [it[0], it[1], it[2]+"/binned_outputs/${params.visium_hd}" ]} )
            versions = versions.mix(SPACEMARKERS_HD.out.versions)
            spacemarkers = SPACEMARKERS_HD.out.LRscores

        } else {
            SPACEMARKERS_SD( ch_sm_inputs )
            versions = versions.mix(SPACEMARKERS_SD.out.versions)
            spacemarkers = SPACEMARKERS_SD.out.spaceMarkersScores
        }
    emit:
        versions
        spacemarkers

}