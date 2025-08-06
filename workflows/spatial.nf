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

ch_multiqc_config          = Channel.fromPath(params.multiqc_config, checkIfExists: true)
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

include { BAYESTME;} from '../subworkflows/local/deconvolve_bayestme'

include { SPACEMARKERS; 
          SPACEMARKERS_MQC;
          SPACEMARKERS_PLOTS;
        } from '../modules/local/spacemarkers/nextflow/main'

include { COGAPS;
          COGAPS_ADATA2DGC; } from '../modules/local/cogaps/main'

include { SQUIDPY_MORANS_I;
          SQUIDPY_SPATIAL_PLOTS;
          SQUIDPY_LIGREC_ANALYSIS } from '../modules/local/squidpy/main'

include { RCTD } from '../modules/local/rctd/rctd'

include { LOAD_DATASET } from '../subworkflows/local/load_dataset'


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

    // Squidpy analysis ch stub
    ch_squidpy = Channel.empty()

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
    ch_input.map { tuple(it[0], it[6]) }.tap { run_cogaps }
    ch_input.map { tuple(it[0], it[7]) }.tap { n_top_genes }
    ch_input.map { tuple(it[0], it[8]) }.tap { spatial_transcriptional_programs }
    ch_input.map { tuple(it[0], it[9]) }.tap { run_spacemarkers }
    ch_input.map { tuple(it[0], it[10]) }.tap { find_annotations }

    // Grab datasets
    LOAD_DATASET(ch_input.map { tuple(it[0], it[1], it[4], it[10]) }) //[meta, data_directory, expression_profiles, find_annotations]

    ch_adata = LOAD_DATASET.out.ch_adata
    ch_scrna = LOAD_DATASET.out.ch_scrna
    ch_coda = LOAD_DATASET.out.ch_coda
    data_directory = LOAD_DATASET.out.data_directory
    ch_matched_adata = LOAD_DATASET.out.ch_matched_adata
    ch_versions = ch_versions.mix(LOAD_DATASET.out.versions)

    ch_sm_inputs = ch_sm_inputs.mix(ch_coda.map { coda -> tuple(coda.meta, coda.coda) })
        .join(data_directory)

    // BayestME deconvolution and plots, run only if not hd as the tool does not support it
    if(!params.hd) {
        BAYESTME(ch_input)
        ch_sm_inputs = ch_sm_inputs.mix(BAYESTME.out.ch_deconvolved.map { tuple(it[0], it[1]) }
                                   .join(data_directory))
        ch_squidpy = ch_squidpy.mix(BAYESTME.out.ch_deconvolved)
        .map { tuple(it[0], it[1]) }
    }

    // RCTD reference-based deconvolution and plots
    // plots are temporary as there is the idea to plot
    // all deconvolution stats with the same process/workflow
    ch_rctd_input = data_directory
        .join( ch_scrna )
        .join( ch_matched_adata )
        .map { tuple(it[0], it[2], it[-1]) }
        .join( n_top_genes)
    
    RCTD( ch_rctd_input )
    ch_versions = ch_versions.mix(RCTD.out.versions)
    ch_sm_inputs = ch_sm_inputs.mix(RCTD.out.rctd_cell_types.map { tuple(it[0], it[1]) }
        .join(data_directory))

    //CoGAPS
    ch_convert_adata = ch_adata
        .join(run_cogaps)
        .filter { it -> it[-1] == true }
        .map( it -> tuple(it[0], it[1]) )
    
    COGAPS_ADATA2DGC( ch_convert_adata )
    ch_versions = ch_versions.mix(COGAPS_ADATA2DGC.out.versions)

    ch_cogaps = COGAPS_ADATA2DGC.out.dgCMatrix.map { tuple(it[0], it[1]) }
        .join(n_cell_types)
        .map { tuple(it[0], it[1], [niterations:params.niterations,
                                           npatterns:it[-1],
                                           sparse:params.sparse,
                                           distributed:params.distributed,
                                           nsets:params.nsets,
                                           nthreads:params.nthreads]) }

    COGAPS(ch_cogaps)
    ch_versions = ch_versions.mix(COGAPS.out.versions)
    ch_sm_inputs = ch_sm_inputs.mix(COGAPS.out.cogapsResult.map { tuple(it[0], it[1]) }.join(data_directory))
    ch_sm_inputs = ch_sm_inputs.combine(run_spacemarkers, by:0)
        .filter { it -> it[-1] == true }                             // make spacemarkers optional
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
    //ch_multiqc_files = ch_multiqc_files.mix(SPACEMARKERS_MQC.out.spacemarkers_mqc.map { it[1] })

    // squidpy analysis 
    ch_squidpy = RCTD.out.rctd_adata
        .map { tuple(it[0], it[1]) }

    SQUIDPY_MORANS_I( ch_squidpy )
    ch_versions = ch_versions.mix(SQUIDPY_MORANS_I.out.versions)

    SQUIDPY_SPATIAL_PLOTS( ch_squidpy )
    ch_versions = ch_versions.mix(SQUIDPY_SPATIAL_PLOTS.out.versions)

    SQUIDPY_LIGREC_ANALYSIS( ch_squidpy )
    ch_versions = ch_versions.mix(SQUIDPY_LIGREC_ANALYSIS.out.versions)

    ch_ligrec_output = SQUIDPY_LIGREC_ANALYSIS.out.ligrec_interactions.collect()

    //collate versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'versions.yml',
            sort: true,
            newLine: true
        ).set { version_yaml }
    ch_multiqc_files = ch_multiqc_files.mix(version_yaml)

    // MultiQC
    // NOTE - will fail to find spaceranger reports unless the full path is provided
    // multiqc does not find spaceranger report for VisiumHD, address with #24
    // ch_multiqc_files.view()
    // MULTIQC (
    //          ch_multiqc_files.collect().ifEmpty([]),
    //          ch_multiqc_config,
    //          [],
    //          [],
    //          [],
    //          []
    //          )
    // multiqc_report = MULTIQC.out.report.toList()


    // emit:
    // multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
    //   versions       = ch_versions                 // channel: [ versions.yml ]


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
