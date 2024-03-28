include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK(samplesheet)
            .csv
            .splitCsv(header:true, sep:',')
            .map{ row -> [
                sample_name: row.sample,
                data_directory: row.data_directory,
                n_cell_types: row.n_cell_types,
                bleeding_correction: row.bleeding_correction.toBoolean(),
                spatial_transcriptional_programs: row.spatial_transcriptional_programs.toBoolean(),
                expression_profile: (row.expression_profile == null || row.expression_profile == "") ? [] : row.expression_profile
            ] }
            .set{ datasets }


    emit:
        datasets
}
