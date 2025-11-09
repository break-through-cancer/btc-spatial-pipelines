include { ATLAS_GET;
          ATLAS_MATCH;
          ADATA_FROM_VISIUM;
          ADATA_FROM_VISIUM_HD;
          ADATA_FROM_SEGMENTED_VISIUM} from '../../modules/local/util/'


workflow LOAD_DATASET {
    take:
        ch_input // channel of tuples: (meta, data_directory, expression_profiles, ...)

    main:
        versions = channel.empty() // channel to collect versions of the tools used

        // Load visium HD or standard data
        if(params.visium_hd) {

            if(params.visium_hd == 'segmented'){
            //use Spaceranger4.0 segmented outputs
                ADATA_FROM_SEGMENTED_VISIUM( ch_input.map {it -> tuple(it.meta, it.data_directory) } )
                ch_adata = ADATA_FROM_SEGMENTED_VISIUM.out.adata
                data_directory = ch_input.map{ it -> tuple(it.meta, it.data_directory) }
                versions = versions.mix(ADATA_FROM_SEGMENTED_VISIUM.out.versions)

            } else {
                ADATA_FROM_VISIUM_HD( ch_input.map {it -> tuple(it.meta, it.data_directory) } )
                ch_adata = ADATA_FROM_VISIUM_HD.out.adata
                data_directory = ch_input.map{ it -> tuple(it.meta, it.data_directory + "/binned_outputs/${params.visium_hd}") }
                versions = versions.mix(ADATA_FROM_VISIUM_HD.out.versions)
            }

        } else {
            ADATA_FROM_VISIUM( ch_input.map {it -> tuple(it.meta, it.data_directory) } )
            ch_adata = ADATA_FROM_VISIUM.out.adata
            data_directory = ch_input.map{ it -> tuple(it.meta, it.data_directory) }
            versions = versions.mix(ADATA_FROM_VISIUM.out.versions)
        }

        // If an atlas has been provided download and prepare it
        if (params.ref_scrna) {
            ATLAS_GET(params.ref_scrna)
            ch_scrna = ch_input.map{ it -> tuple(it.meta)}
                .combine(ATLAS_GET.out.atlas)
            versions = versions.mix(ATLAS_GET.out.versions)
        } else{
        // Use the matched scrna per sample
            expression_profiles = ch_input.map { it -> tuple(it.meta, it.expression_profile) }
            ch_scrna = expression_profiles.filter { it -> it[1]!=null && it[1]!='' && it[1]!=[] }
        }

        // match scRNA atlas to spatial data
        ATLAS_MATCH(ch_scrna.join( ch_adata ))
        ch_matched_adata = ATLAS_MATCH.out.adata_matched
        versions = versions.mix(ATLAS_MATCH.out.versions)


    emit:
        ch_adata
        ch_matched_adata
        ch_scrna
        data_directory
        versions
}