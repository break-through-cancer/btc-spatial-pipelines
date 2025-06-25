include { BAYESTME_FILTER_GENES;
          BAYESTME_BLEEDING_CORRECTION;
          BAYESTME_DECONVOLUTION;
          BAYESTME_LOAD_SPACERANGER;
        } from '../../modules/local/bayestme/nextflow/subworkflows/bayestme/bayestme_basic_visium_analysis/main'


workflow BAYESTME {
    take:
        ch_input // Channel of tuples: (meta, adata_sc, adata_spatial, coda_annotations)

    main:

        ch_input.map { tuple(it[0], it[1]) }.tap { data_directory }
        ch_input.map { tuple(it[0], it[2]) }.tap { n_cell_types }
        ch_input.map { tuple(it[0], it[3]) }.tap { should_run_bleeding_correction }
        ch_input.map { tuple(it[0], it[4]) }.tap { expression_profiles }
        ch_input.map { tuple(it[0], it[5]) }.tap { run_bayestme }
        ch_input.map { tuple(it[0], it[6]) }.tap { run_cogaps }
        ch_input.map { tuple(it[0], it[7]) }.tap { n_top_genes }
        ch_input.map { tuple(it[0], it[8]) }.tap { spatial_transcriptional_programs }
        ch_input.map { tuple(it[0], it[9]) }.tap { run_spacemarkers }
        ch_input.map { tuple(it[0], it[10]) }.tap { find_annotations }

        ch_btme = ch_input.map { tuple(it[0], it[1], it[2]) } // meta, data_directory
                          .join( run_bayestme )
                          .filter { it -> it[-1] == true }   // run_bayestme

        BAYESTME_LOAD_SPACERANGER( ch_input.map { tuple(it[0], it[1]) })
        ch_adata = BAYESTME_LOAD_SPACERANGER.out.adata

        filter_genes_input = ch_adata
        .join( n_top_genes )
        .map { tuple(
            it[0],               // sample_name
            it[1],               // adata
            true,                // filter_ribosomal_genes
            it[2],               // n_top_genes
            0.9,                 // spot_threshold
            [])                  // disabled ch_scrna, see #65
        }

        BAYESTME_FILTER_GENES( filter_genes_input )
        ch_versions = BAYESTME_FILTER_GENES.out.versions

        BAYESTME_FILTER_GENES.out.adata_filtered
            .join( should_run_bleeding_correction )
            .filter { it[-1] == true }
            .map { tuple(it[0], it[1]) }
            .tap { bleeding_correction_input }

        BAYESTME_FILTER_GENES.out.adata_filtered
            .join( should_run_bleeding_correction )
            .filter { it[-1] == false }
            .map { tuple(it[0], it[1]) }
            .join( ch_input.map { tuple(it[0], it[2]) } )
            .map { tuple(it[0], it[1], it[2], params.bayestme_spatial_smoothing_parameter) }
            .tap { not_bleed_corrected_deconvolution_input }

        BAYESTME_BLEEDING_CORRECTION( bleeding_correction_input )
        ch_versions = ch_versions.mix(BAYESTME_BLEEDING_CORRECTION.out.versions)

        deconvolution_input = BAYESTME_BLEEDING_CORRECTION.out.adata_corrected
            .join( ch_input.map { tuple(it[0], it[2]) } )
            .map { tuple(it[0], it[1], it[2], params.bayestme_spatial_smoothing_parameter) }
            .concat( not_bleed_corrected_deconvolution_input )
            .map { tuple(it[0], //meta
                        it[1], //dataset_filtered
                        it[2], //n_cell_types
                        it[3], //smoothing_parameter
                        [])    //expression truth placeholder, see #65
                }
        BAYESTME_DECONVOLUTION( deconvolution_input )
        ch_versions = ch_versions.mix(BAYESTME_DECONVOLUTION.out.versions)

        // Output the deconvolved data
        ch_deconvolved = BAYESTME_DECONVOLUTION.out.adata_deconvolved

    emit:
        ch_deconvolved
        ch_versions
}