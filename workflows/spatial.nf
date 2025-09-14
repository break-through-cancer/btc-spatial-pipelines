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

include { SPACEMARKERS_HD;  // temp - allow spacemarkers to run on dev
          SPACEMARKERS_HD_PLOTS;
        } from '../modules/local/spacemarkers/nextflow/main_hd'

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


    // Load input paths and metadata
    INPUT_CHECK (
        file(params.input)
    )

    // Grab datasets
    ch_datasets = INPUT_CHECK.out.datasets
    LOAD_DATASET(ch_datasets)

    ch_adata = LOAD_DATASET.out.ch_adata
    ch_scrna = LOAD_DATASET.out.ch_scrna
    ch_coda = LOAD_DATASET.out.ch_coda
    data_directory = LOAD_DATASET.out.data_directory
    ch_matched_adata = LOAD_DATASET.out.ch_matched_adata
    ch_versions = ch_versions.mix(LOAD_DATASET.out.versions)

    ch_sm_inputs = ch_sm_inputs.mix(ch_coda.map { coda -> tuple(coda.meta, coda.coda) })
        .join(data_directory)

    // BayestME deconvolution and plots, run only if not hd as the tool does not support it
    if(!params.visium_hd && params.deconvolve.bayestme) {
        BAYESTME(ch_datasets)
        ch_sm_inputs = ch_sm_inputs.mix(BAYESTME.out.ch_deconvolved.map { tuple(it[0], it[1]) }
                                   .join(data_directory))
        ch_squidpy = ch_squidpy.mix(BAYESTME.out.ch_deconvolved)
        .map { tuple(it[0], it[1]) }
    }

    // RCTD reference-based deconvolution and plots
    // plots are temporary as there is the idea to plot
    // all deconvolution stats with the same process/workflow
    if(params.deconvolve.rctd) {
        ch_rctd_input = data_directory
            .join( ch_scrna )
            .join( ch_matched_adata )
            .map { tuple(it[0], it[2], it[-1]) }

        RCTD( ch_rctd_input )
        ch_versions = ch_versions.mix(RCTD.out.versions)
        ch_sm_inputs = ch_sm_inputs.mix(RCTD.out.rctd_cell_types.map { tuple(it[0], it[1]) }
            .join(data_directory))
        
        ch_squidpy = ch_squidpy.mix(RCTD.out.rctd_adata)
            .map { tuple(it[0], it[1]) }
    }

    //CoGAPS reference-free spatially unaware deconvolution
    if(params.deconvolve.cogaps) {
        ch_convert_adata = ch_adata
    
        COGAPS_ADATA2DGC( ch_convert_adata )
        ch_versions = ch_versions.mix(COGAPS_ADATA2DGC.out.versions)

        ch_cogaps = COGAPS_ADATA2DGC.out.dgCMatrix.map { tuple(it[0], it[1]) }
            .map { tuple(it[0], it[1], [niterations:params.niterations,
                                            npatterns:it[-1],
                                            sparse:params.sparse,
                                            distributed:params.distributed,
                                            nsets:params.nsets,
                                            nthreads:params.nthreads]) }

        COGAPS(ch_cogaps)
        ch_versions = ch_versions.mix(COGAPS.out.versions)
        ch_sm_inputs = ch_sm_inputs.mix(COGAPS.out.cogapsResult.map { tuple(it[0], it[1]) }.join(data_directory))
    }

    //spacemarkers - main
    if (params.analyze.spacemarkers){
        if(params.visium_hd) {

            SPACEMARKERS_HD( ch_sm_inputs )   //temp - allow spacemarkers to run on dev

        } else {

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

        }
    }


    // squidpy analysis
    if (params.analyze.squidpy){
        SQUIDPY_MORANS_I( ch_squidpy )
        ch_versions = ch_versions.mix(SQUIDPY_MORANS_I.out.versions)

        SQUIDPY_SPATIAL_PLOTS( ch_squidpy )
        ch_versions = ch_versions.mix(SQUIDPY_SPATIAL_PLOTS.out.versions)

        SQUIDPY_LIGREC_ANALYSIS( ch_squidpy )
        ch_versions = ch_versions.mix(SQUIDPY_LIGREC_ANALYSIS.out.versions)

        ch_ligrec_output = SQUIDPY_LIGREC_ANALYSIS.out.ligrec_interactions.collect()
    }

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
