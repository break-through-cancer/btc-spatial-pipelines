include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK(samplesheet)
            .csv
            .splitCsv(header:true, sep:',')
            .map{ row -> tuple row.sample, row.data_directory }
            .set{ reads }

    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
        versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
