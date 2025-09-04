include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK(samplesheet)
            .csv
            .splitCsv(header:true, sep:',')
            .map{ row -> [
                id: row.sample,
                data_directory: file(row.data_directory),
                n_cell_types: row.n_cell_types,
                bleeding_correction: row.bleeding_correction.toBoolean(),
                expression_profile: (row.expression_profile == null || row.expression_profile == "") ? [] : row.expression_profile,
                run_bayestme: row.run_bayestme.toBoolean(),
                run_cogaps: row.run_cogaps.toBoolean(),
                cogaps_niterations: row.cogaps_niterations,
                n_top_genes: row.n_top_genes,
                run_spacemarkers: row.run_spacemarkers.toBoolean(),
                find_annotations: row.find_annotations.toBoolean()
            ] }
            .set{ datasets }


    emit:
        datasets
}
