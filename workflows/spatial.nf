/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { DECONVOLVE } from '../subworkflows/local/deconvolve'

include { SPACEMARKERS; 
          SPACEMARKERS_MQC;
          SPACEMARKERS_PLOTS;
        } from '../modules/local/spacemarkers/nextflow/main'

include { SPACEMARKERS_HD;  // temp - allow spacemarkers to run on dev
          SPACEMARKERS_HD_PLOTS;
        } from '../modules/local/spacemarkers/nextflow/main_hd'

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
    //multiqc_report   = Channel.empty()
    ch_versions = Channel.empty()

    // Load input paths and metadata
    INPUT_CHECK (
        file(params.input)
    )

    // Grab datasets
    ch_datasets = INPUT_CHECK.out.datasets
    LOAD_DATASET(ch_datasets)

    ch_adata = LOAD_DATASET.out.ch_adata
    ch_scrna = LOAD_DATASET.out.ch_scrna
    ch_matched_adata = LOAD_DATASET.out.ch_matched_adata
    ch_versions = ch_versions.mix(LOAD_DATASET.out.versions)

    // Deconvolve / cell type / use external
    DECONVOLVE (ch_datasets, ch_matched_adata, ch_scrna, ch_adata)
    ch_versions = ch_versions.mix(DECONVOLVE.out.ch_versions)
  
    ch_squidpy = Channel.empty()
    //Analyze - spacemarkers
    if (params.analyze.spacemarkers){

        //csv if available, otherwise deconvolution object - SpaceMarkers knows how to handle both
        ch_sm_inputs = DECONVOLVE.out.ch_deconvolved
                        .map { [it.meta, it.cell_probs?:it.obj] }
                        .combine( ch_datasets.map { [it.meta, it.data_directory] }, by:0 )

        if(params.visium_hd) {
        
            SPACEMARKERS_HD( ch_sm_inputs.map {[it[0], it[1], it[2]+"/binned_outputs/${params.visium_hd}" ]} )   //temp - allow spacemarkers to run on dev

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
        // spatially variable genes do not depend on deconvolution
        SQUIDPY_MORANS_I( ch_adata )
        ch_versions = ch_versions.mix(SQUIDPY_MORANS_I.out.versions)

        // squidpy now anndata to plot spatial plots and ligrec
        ch_squidpy = DECONVOLVE.out.ch_deconvolved
                            .filter { it -> it.obj != null }
                            .filter { it -> it.obj.name.endsWith('.h5ad') }
                            .map { [it.meta, it.obj] }

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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
