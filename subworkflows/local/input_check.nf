include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK(samplesheet)
            .csv
            .splitCsv(header:true, sep:',')
            .map{ row -> [sample_name: row.sample, data_directory: row.data_directory, n_components: row.n_components] }
            .set{ datasets }


    emit:
        datasets
}
