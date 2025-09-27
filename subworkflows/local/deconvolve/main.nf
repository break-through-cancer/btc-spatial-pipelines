
include { BAYESTME } from './deconvolve_bayestme'
include { COGAPS } from './deconvolve_cogaps'
include { RCTD } from '../../../modules/local/rctd/rctd'

workflow DECONVOLVE {
    take:
        ch_datasets          // from LOAD_DATASET, for BayesTME
        ch_matched_adata     // from LOAD_DATASET, for RCTD
        ch_scrna             // from LOAD_DATASET, for RCTD
        ch_adata             // from LOAD_DATASET, for CoGAPS
    main:

    versions = Channel.empty()
    ch_deconvolved = Channel.empty()  //outputs: meta, optional cell type probs, optional deconvolution object

    // Grab external deconvolution results
    // CODA annotation channel - or any other external csv annotaion
        if (params.deconvolve.external){
            ch_deconvolved = ch_datasets.map { tuple(it.meta, it.data_directory) }
                .flatMap { item -> 
                    def meta = item[0]
                    def data_path = item[1]
                    def coda_files = []
                    data_path.eachFileRecurse { file ->
                        if (file.name.endsWith('cellular_compositions.csv')) {
                            coda_files.add(file)
                        }
                    }
                    coda_files.collect { file -> [meta: meta, cell_probs:file, obj: null] }
                }
        }

    // BayestME deconvolution and plots, run only if not hd as the tool does not support it
    if(!params.visium_hd && params.deconvolve.bayestme) {
        BAYESTME(ch_datasets)
        ch_deconvolved = ch_deconvolved.mix(BAYESTME.out.ch_deconvolved.map { [meta:it[0], cell_probs:null, obj:it[1]] })
        versions = versions.mix(BAYESTME.out.versions)
    }

    // RCTD reference-based deconvolution and plots
    // plots are temporary as there is the idea to plot
    // all deconvolution stats with the same process/workflow
    if(params.deconvolve.rctd) {
        ch_rctd_input = ch_scrna.join(ch_matched_adata)
        RCTD( ch_rctd_input )
        versions = versions.mix(RCTD.out.versions)
        ch_rctd_output = RCTD.out.rctd_cell_types
                            .join(RCTD.out.rctd_adata)
                            .map { [meta:it[0], cell_probs:it[1], obj:it[2]] }
        ch_deconvolved = ch_deconvolved.mix(ch_rctd_output)
        versions = versions.mix(RCTD.out.versions)
    }

    //CoGAPS reference-free spatially unaware deconvolution
    if(params.deconvolve.cogaps) {
        COGAPS( ch_adata )
        versions = versions.mix(COGAPS.out.ch_versions)
        ch_cogaps_output = COGAPS.out.ch_deconvolved
                                .map { [meta:it[0], cell_probs:null, obj:it[1]] }
        ch_deconvolved = ch_deconvolved.mix(ch_cogaps_output) //need to add csv with cell type probs
    }

    emit:
        ch_deconvolved
        versions
}