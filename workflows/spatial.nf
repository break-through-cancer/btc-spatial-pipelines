/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSpatial.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { BAYESTME_LOAD_SPACERANGER;
          BAYESTME_FILTER_GENES;
          BAYESTME_BLEEDING_CORRECTION;
          BAYESTME_DECONVOLUTION;
          BAYESTME_SPATIAL_TRANSCRIPTIONAL_PROGRAMS;
        } from '../modules/local/bayestme/nextflow/subworkflows/bayestme/bayestme_basic_visium_analysis/main'
        
include { SPACEMARKERS; 
          SPACEMARKERS_MQC;
          SPACEMARKERS_PLOTS;
        } from '../modules/local/spacemarkers/nextflow/main'

include { COGAPS;
          COGAPS_ADATA2DGC; } from '../modules/local/cogaps/nextflow/main'

include { SQUIDPY } from '../modules/local/squidpy/main'

include { MATCH_ADATAS } from '../modules/local/util/match_adatas.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MULTIQC } from '../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SPATIAL {

    //gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    ch_versions = Channel.empty()

    // Optional inputs to SpaceMarkers
    ch_sm_inputs = Channel.empty()

    INPUT_CHECK (
        file(params.input)
    )
    
    // NOTE: append to the list to avoid other indices being off
    ch_input = INPUT_CHECK.out.datasets.map { tuple(
        id:it.sample_name,
        it.data_directory,
        it.n_cell_types,
        it.bleeding_correction,
        it.expression_profile,
        it.run_bayestme,
        it.run_cogaps,
        it.n_top_genes,
        it.spatial_transcriptional_programs,
        it.run_spacemarkers,
        it.find_annotations
    ) }

    ch_input.map { tuple(it[0], it[1]) }.tap { data_directory }
    ch_input.map { tuple(it[0], it[2]) }.tap { n_cell_types }
    ch_input.map { tuple(it[0], it[3]) }.tap { should_run_bleeding_correction }
    ch_input.map { tuple(it[0], it[4]) }.tap { expression_profiles }
    ch_input.map { tuple(it[0], it[5]) }.tap { run_bayestme }
    ch_input.map { tuple(it[0], it[7]) }.tap { n_top_genes }
    ch_input.map { tuple(it[0], it[8]) }.tap { spatial_transcriptional_programs }
    ch_input.map { tuple(it[0], it[9]) }.tap { run_spacemarkers }
    ch_input.map { tuple(it[0], it[10]) }.tap { find_annotations }

    // A new channel that contains *.html spaceranger reports for multiqc
    ch_sr_reports = data_directory.flatMap { item ->
        def meta = item[0]
        def data_path = item[1]
        def html_files = file(data_path).listFiles().findAll { it.name.endsWith('.html') }
        html_files.collect { file -> [meta: meta, sr_report: file] }
    }


    // CODA annotation channel
    ch_coda = data_directory
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

    // SCRNA reference data channel
    ch_scrna = data_directory
        .join(expression_profiles)
        .filter { it -> it[2].length()>0 } // only run if a profile is provided
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
            scrna_files.collect { file -> [meta, file] }
        }
    
    ch_sm_inputs = ch_sm_inputs.mix(ch_coda.map { coda -> tuple(coda.meta, coda.coda) })
        .join(data_directory)

    ch_multiqc_files = ch_multiqc_files.mix(ch_sr_reports.map { it.sr_report })

    ch_btme = ch_input.map { tuple(it[0], it[1]) }
    BAYESTME_LOAD_SPACERANGER( ch_btme )
    ch_versions = ch_versions.mix(BAYESTME_LOAD_SPACERANGER.out.versions)

    //Match the gene indices of the spatial data with the scRNA data
    ch_scrna = ch_scrna
        .join(BAYESTME_LOAD_SPACERANGER.out.adata)
    MATCH_ADATAS( ch_scrna.map { tuple(it[0], it[1], it[2]) } )

    filter_genes_input = BAYESTME_LOAD_SPACERANGER.out.adata
        .join( n_top_genes )
        .map { tuple(
            it[0],               // sample_name
            it[1],               // adata
            true,                // filter_ribosomal_genes
            it[2],               // n_top_genes
            0.9)                 // spot_threshold
    }.join(MATCH_ADATAS.out.adata_sc_matched)

    BAYESTME_FILTER_GENES( filter_genes_input )
    ch_versions = ch_versions.mix(BAYESTME_FILTER_GENES.out.versions)

    BAYESTME_FILTER_GENES.out.adata_filtered
        .join( should_run_bleeding_correction )
        .filter { it[2] == true }
        .map { tuple(it[0], it[1]) }
        .tap { bleeding_correction_input }

    BAYESTME_FILTER_GENES.out.adata_filtered
        .join( should_run_bleeding_correction )
        .filter { it[2] == false }
        .map { tuple(it[0], it[1]) }
        .join( ch_input.map { tuple(it[0], it[2]) } )
        .map { tuple(it[0], it[1], it[2], params.bayestme_spatial_smoothing_parameter) }
        .join(MATCH_ADATAS.out.adata_sc_matched)
        .tap { not_bleed_corrected_deconvolution_input }

    BAYESTME_BLEEDING_CORRECTION( bleeding_correction_input )
    ch_versions = ch_versions.mix(BAYESTME_BLEEDING_CORRECTION.out.versions)

    deconvolution_input = BAYESTME_BLEEDING_CORRECTION.out.adata_corrected
        .join( ch_input.map { tuple(it[0], it[2]) } )
        .map { tuple(it[0], it[1], it[2], params.bayestme_spatial_smoothing_parameter) }
        .join(MATCH_ADATAS.out.adata_sc_matched)
        .concat( not_bleed_corrected_deconvolution_input )
        .combine( run_bayestme )
        .filter { it -> it[-1] == true }   // run_bayestme
        .map { tuple(it[0], it[1], it[2], it[3], it[4]) }

    BAYESTME_DECONVOLUTION( deconvolution_input )
    ch_versions = ch_versions.mix(BAYESTME_DECONVOLUTION.out.versions)

    BAYESTME_DECONVOLUTION.out.adata_deconvolved
        .join(BAYESTME_DECONVOLUTION.out.deconvolution_samples)
        .join( spatial_transcriptional_programs )
        .filter { it -> it[3] == true }                           // spatial_transcriptional_programs bool
        .map { tuple(it[0], it[1], it[2], []) }
        .tap { stp_input }
    ch_versions = ch_versions.mix(BAYESTME_DECONVOLUTION.out.versions)
    ch_sm_inputs = ch_sm_inputs.mix(BAYESTME_DECONVOLUTION.out.adata_deconvolved.map { tuple(it[0], it[1]) }
        .join(data_directory))


    BAYESTME_SPATIAL_TRANSCRIPTIONAL_PROGRAMS( stp_input )
    ch_versions = ch_versions.mix(BAYESTME_SPATIAL_TRANSCRIPTIONAL_PROGRAMS.out.versions)

    // squidpy - spatially variable genes
    ch_svgs = BAYESTME_BLEEDING_CORRECTION.out.adata_corrected
        .concat( not_bleed_corrected_deconvolution_input )
        .map { tuple(it[0], it[1]) }
    SQUIDPY( ch_svgs )
    ch_versions = ch_versions.mix(SQUIDPY.out.versions)

    //cogaps - make use of the BTME preprocessing
    ch_samplesheet = ch_input.map { tuple(it[0], it[1]) }

    ch_cogaps_input = BAYESTME_BLEEDING_CORRECTION.out.adata_corrected
        .concat( not_bleed_corrected_deconvolution_input )
        .map( it -> tuple(it[0], it[1]) )

    COGAPS_ADATA2DGC( ch_cogaps_input )
    ch_versions = ch_versions.mix(COGAPS_ADATA2DGC.out.versions)

    ch_gaps = INPUT_CHECK.out.datasets
        .filter { it -> it.run_cogaps == true }
        .map { tuple([id:it.sample_name], [niterations:params.cogaps_niterations,
                                           npatterns:it.n_cell_types,
                                           sparse:params.cogaps_sparse,
                                           distributed:params.cogaps_distributed,
                                           nsets:params.cogaps_nsets,
                                           nthreads:params.cogaps_nthreads]) }
        .join(COGAPS_ADATA2DGC.out.dgCMatrix.map { tuple(it[0], it[1]) })
        .map { tuple(it[0], it[2], it[1]) }                          // reorder to match cogaps input

    COGAPS(ch_gaps)
    ch_versions = ch_versions.mix(COGAPS.out.versions)

    ch_sm_inputs = ch_sm_inputs.mix(COGAPS.out.cogapsResult.map { tuple(it[0], it[1]) }.join(ch_samplesheet))
    ch_sm_inputs = ch_sm_inputs.combine(run_spacemarkers, by:0)
        .filter { it -> it[3] == true }                             // make spacemarkers optional
        .map { tuple(it[0], it[1], it[2]) }

    //spacemarkers - main
    SPACEMARKERS( ch_sm_inputs )
    ch_versions = ch_versions.mix(SPACEMARKERS.out.versions)


    //spacemarkers - plots
    ch_plotting_input = SPACEMARKERS.out.spaceMarkersScores
        .map { tuple(it[0], it[1]) }
    ch_plotting_input = ch_plotting_input.join(SPACEMARKERS.out.overlapScores)
        .map { tuple(it[0], it[1], it[2], it[3]) }
    
    SPACEMARKERS_PLOTS( ch_plotting_input)
    ch_versions = ch_versions.mix(SPACEMARKERS_PLOTS.out.versions)

    //spacemarkers - mqc
    SPACEMARKERS_MQC( SPACEMARKERS.out.spaceMarkers.map { tuple(it[0], it[1], it[2]) } )
    ch_versions = ch_versions.mix(SPACEMARKERS_MQC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SPACEMARKERS_MQC.out.spacemarkers_mqc.map { it[1] })


    //collate versions
    version_yaml = Channel.empty()
    version_yaml = softwareVersionsToYAML(ch_versions)
                   .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
    

    // MultiQC
    // NOTE - will fail to find spaceranger reports unless the full path is provided
    MULTIQC (
            ch_multiqc_files.collect().ifEmpty([]),[],[],[],[],[]
            )
    multiqc_report = MULTIQC.out.report.toList()


    emit:
    multiqc_report // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
