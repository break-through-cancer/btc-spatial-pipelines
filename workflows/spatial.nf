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
        } from '../modules/local/spacemarkers/nextflow/main'

include  { QC } from '../modules/local/util/util'

include { SPACEMARKERS_HD;  // temp - spacemarkers hd runs on dev
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
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    versions = Channel.empty()

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
    versions = versions.mix(LOAD_DATASET.out.versions)

    // Report read data to MultiQC
    ch_report = ch_scrna.map {it -> [it[0], it[1], 'atlas_input']} 
    ch_report = ch_report.mix(ch_adata.map {it -> [it[0], it[1], 'adata_input']})


    // Deconvolve / cell type / use external
    DECONVOLVE (ch_datasets, ch_matched_adata, ch_scrna, ch_adata)
    versions = versions.mix(DECONVOLVE.out.versions)

    // Report deconvolution results to MultiQC
    ch_report = ch_report.mix( DECONVOLVE.out.ch_deconvolved.map {it -> [it.meta, it.obj, 'data_output']} )

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
            versions = versions.mix(SPACEMARKERS.out.versions)

        }
    }


    // squidpy analysis
    if (params.analyze.squidpy){
        // spatially variable genes do not depend on deconvolution
        SQUIDPY_MORANS_I( ch_adata )
        versions = versions.mix(SQUIDPY_MORANS_I.out.versions)

        // squidpy now anndata to plot spatial plots and ligrec
        ch_squidpy = DECONVOLVE.out.ch_deconvolved
                            .filter { it -> it.obj != null }
                            .filter { it -> it.obj.name.endsWith('.h5ad') }
                            .map { [it.meta, it.obj] }

        SQUIDPY_SPATIAL_PLOTS( ch_squidpy )
        versions = versions.mix(SQUIDPY_SPATIAL_PLOTS.out.versions)

        SQUIDPY_LIGREC_ANALYSIS( ch_squidpy )
        versions = versions.mix(SQUIDPY_LIGREC_ANALYSIS.out.versions)

        ch_ligrec_output = SQUIDPY_LIGREC_ANALYSIS.out.ligrec_interactions.collect()
    }

    //make reports
    QC( ch_report )
    ch_multiqc_files = ch_multiqc_files.mix(QC.out.report)
    versions = versions.mix(QC.out.versions)

    //collate versions
    softwareVersionsToYAML(versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'all_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set{ all_versions }
    ch_multiqc_files = ch_multiqc_files.mix(all_versions)

    // MultiQC
    // NOTE - will fail to find spaceranger reports unless the full path is provided
    // multiqc does not find spaceranger report for VisiumHD, address with #24
     MULTIQC (
             ch_multiqc_files.collect().ifEmpty([]),
             file(params.multiqc_config).exists() ? file(params.multiqc_config) : null,
             [],
             [],
             [],
             []
             )
    multiqc_report = MULTIQC.out.report.toList()


    //emit:
      multiqc_report = multiqc_report           // channel: /path/to/multiqc_report.html
      versions       = versions                 // channel: [ versions.yml ]


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
