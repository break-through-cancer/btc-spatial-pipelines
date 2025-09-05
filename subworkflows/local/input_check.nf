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
            .map{it + [bleeding_correction: it.bleeding_correction.toBoolean()]}
            .map{it + [run_bayestme: it.run_bayestme.toBoolean()]}
            .map{it + [run_cogaps: it.run_cogaps.toBoolean()]}
            .map{it + [run_spacemarkers: it.run_spacemarkers.toBoolean()]}
            .map{it + [find_annotations: it.find_annotations.toBoolean()]}
            .set{ datasets }

    emit:
        datasets
}
