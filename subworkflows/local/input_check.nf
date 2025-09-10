include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK(samplesheet)
            .csv
            .splitCsv(header:true, sep:',')
            .map{it + [id: it.sample]}
            .map{it + [data_directory: file(it.data_directory)]}
            .map{it + [expression_profile: (it.expression_profile == null || it.expression_profile == "") ? [] : file(it.expression_profile)]}
            .set{ datasets }

    emit:
        datasets
}
