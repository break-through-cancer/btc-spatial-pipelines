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
        
include { SPACEMARKERS } from '../modules/local/jhu-spatial/modules/local/spacemarkers'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SPATIAL {
    ch_versions = Channel.empty()

    INPUT_CHECK (
        file(params.input)
    )

    ch_input = INPUT_CHECK.out.datasets.map { tuple(
        [id:it.sample_name, single_end: false],
        it.data_directory,
        it.n_cell_types,
        it.bleeding_correction,
        it.expression_profile
    ) }

    ch_input.map { tuple(it[0], it[3]) }.tap { should_run_bleeding_correction }
    ch_input.map { tuple(it[0], it[4]) }.tap { expression_profiles }
    ch_input.map { tuple(it[0], it[1]) }.tap { data_directory }

    BAYESTME_LOAD_SPACERANGER( ch_input.map { tuple(it[0], it[1]) } )
    ch_versions = ch_versions.mix(BAYESTME_LOAD_SPACERANGER.out.versions)

    filter_genes_input = BAYESTME_LOAD_SPACERANGER.out.adata.map { tuple(
        it[0],
        it[1],
        true,
        1000,
        0.9)
    }.join(expression_profiles)


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
        .join(expression_profiles)
        .tap { not_bleed_corrected_deconvolution_input }


    BAYESTME_BLEEDING_CORRECTION( bleeding_correction_input )
    ch_versions = ch_versions.mix(BAYESTME_BLEEDING_CORRECTION.out.versions)

    deconvolution_input = BAYESTME_BLEEDING_CORRECTION.out.adata_corrected
        .join( ch_input.map { tuple(it[0], it[2]) } )
        .map { tuple(it[0], it[1], it[2], params.bayestme_spatial_smoothing_parameter) }
        .join(expression_profiles)
        .concat( not_bleed_corrected_deconvolution_input )

    BAYESTME_DECONVOLUTION( deconvolution_input )
    ch_versions = ch_versions.mix(BAYESTME_DECONVOLUTION.out.versions)

    SPACEMARKERS( BAYESTME_DECONVOLUTION.out.adata_deconvolved.map { tuple(it[0], it[1]) }.join(data_directory) )
    ch_versions = ch_versions.mix(SPACEMARKERS.out.versions)

    BAYESTME_DECONVOLUTION.out.adata_deconvolved.join(BAYESTME_DECONVOLUTION.out.deconvolution_samples)
        .map { tuple(it[0], it[1], it[2], []) }
        .tap { stp_input }
    ch_versions = ch_versions.mix(BAYESTME_DECONVOLUTION.out.versions)

    BAYESTME_SPATIAL_TRANSCRIPTIONAL_PROGRAMS( stp_input )
    ch_versions = ch_versions.mix(BAYESTME_SPATIAL_TRANSCRIPTIONAL_PROGRAMS.out.versions)

    //collate versions
    version_yaml = Channel.empty()
    version_yaml = softwareVersionsToYAML(ch_versions)
            .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
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
