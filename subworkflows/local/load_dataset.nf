include { ATLAS_GET;
          ATLAS_MATCH;
          ADATA_FROM_VISIUM;
          ADATA_FROM_VISIUM_HD; } from '../../modules/local/util/util'


workflow LOAD_DATASET {
    take:
        ch_input // Channel of tuples: (meta, data_directory, expression_profiles, find_annotations)

    main:
        data_directory = ch_input.map { it -> tuple(it[0], it[1]) }      // meta, data_directory column from samplesheet
        expression_profiles = ch_input.map { it -> tuple(it[0], it[2]) } // expression_profiles column from samplesheet
        find_annotations = ch_input.map { it -> tuple(it[0], it[3]) }    // find_annotations flag from samplesheet

        // Load visium HD or standard data
        if(params.hd) {
            ADATA_FROM_VISIUM_HD( ch_input.map { tuple(it[0], it[1]) } )
            ch_adata = ADATA_FROM_VISIUM_HD.out.adata
            data_directory = ch_input.map{ it -> tuple(it[0], it[1] + "/binned_outputs/${params.hd}") }
        } else {
            ADATA_FROM_VISIUM( ch_input.map { tuple(it[0], it[1]) } )
            ch_adata = ADATA_FROM_VISIUM.out.adata
            data_directory = ch_input.map{ it -> tuple(it[0], it[1]) }
        }


        // If an atlas has been provided download and prepare it
        if (params.reference_scrna) {
            ATLAS_GET(params.reference_scrna)
            ch_scrna = ch_input.map { tuple(it[0], it[1]) }
                .combine(ATLAS_GET.out.atlas)
                .map { tuple(it[0], it[2]) } // meta, adata_sc
        } else{
        // Look for matched scRNA deconvolution files using file mask
        // from expression_profiles column of the samplesheet
            ch_scrna = ch_input.map { tuple(it[0], it[1]) }
                .join(expression_profiles)
                .filter { it -> it[2].size()>0 } // only run if a profile is provided
                .flatMap { item -> 
                    def meta = item[0]
                    def data_path = item[1]
                    def seach_expr = item[2]
                    def scrna_files = []
                    data_path.eachFileRecurse { file ->
                        if (file.name.endsWith(seach_expr)) {
                            scrna_files.add(file)
                        }
                    }
                    scrna_files.collect { file -> [meta, file] } // meta, adata_sc
                }
        }

        // match scRNA atlas to spatial data
        ATLAS_MATCH(ch_scrna.join( ch_adata ))
        ch_matched_adata = ch_input.join(ATLAS_MATCH.out.adata_matched)
            .map(it -> tuple(it[0], it[-1])) // meta, adata_sc

        // CODA annotation channel - or any other external csv annotaion
        ch_coda = ch_input.map { tuple(it[0], it[1]) }
            .join(find_annotations)
            .filter { it -> it[2]==true} // only run if find_annotations is true
            .flatMap { item -> 
                def meta = item[0]
                def data_path = item[1]
                def coda_files = []
                data_path.eachFileRecurse { file ->
                    if (file.name.endsWith('tissue_positions_cellular_compositions.csv')) {
                        coda_files.add(file)
                    }
                }
                coda_files.collect { file -> [meta: meta, coda: file] }
            }

    emit:
        ch_adata
        ch_matched_adata  // TODO: think if this needs to get returned in ch_adata instead
        ch_scrna
        ch_coda
        data_directory
}