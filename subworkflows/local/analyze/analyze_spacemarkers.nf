include { SPACEMARKERS_HD } from '../../../modules/local/spacemarkers/nextflow/main_hd'
include { SPACEMARKERS as SPACEMARKERS_SD } from '../../../modules/local/spacemarkers/nextflow/main'
include { SPACEMARKERS_HARMONIZE as SPACEMARKERS_HARMONIZE_IMSCORES } from '../../../modules/local/util/'
include { SPACEMARKERS_HARMONIZE as SPACEMARKERS_HARMONIZE_LRSCORES } from '../../../modules/local/util/'


workflow SPACEMARKERS {

    take:
        ch_sm_inputs   // from DECONVOLVE, csv if available or (BayesTME, CoGAPS) obj - SpaceMarkers knows how to handle both
    main:
        versions = channel.empty()
        lrscores_raw = channel.empty() // because HD outputs both imscores and lrscores but SD only imscores

        if(params.visium_hd) {
            SPACEMARKERS_HD( ch_sm_inputs.map {it -> [it[0], it[1], it[2]+"/binned_outputs/${params.visium_hd}" ]} )
            // IMScores: IMScores.rds with row names holding gene names,
            // followed by cell_type1_near_cell_typeN columns, values are IMScores
            // LRscores: LRscores.rds with row names holding ligand-receptor pair names,
            // followed by cell_type1_near_cell_typeN columns, values are LRscores
            versions = versions.mix(SPACEMARKERS_HD.out.versions)
            imscores_raw = SPACEMARKERS_HD.out.IMscores
            lrscores_raw = SPACEMARKERS_HD.out.LRscores

        } else {
            SPACEMARKERS_SD( ch_sm_inputs )
            // IMScores: spacemarkers.csv first column is Gene with gene name, 
            // followed by cell_type1_cell_typeN columns, values are spacemarkers
            versions = versions.mix(SPACEMARKERS_SD.out.versions)
            imscores_raw = SPACEMARKERS_SD.out.spaceMarkersScores
        }

    SPACEMARKERS_HARMONIZE_IMSCORES( imscores_raw )
    SPACEMARKERS_HARMONIZE_LRSCORES ( lrscores_raw )

    imscores = SPACEMARKERS_HARMONIZE_IMSCORES.out.spacemarkers
    lrscores = SPACEMARKERS_HARMONIZE_LRSCORES.out.spacemarkers

    emit:
        versions
        imscores
        lrscores

}