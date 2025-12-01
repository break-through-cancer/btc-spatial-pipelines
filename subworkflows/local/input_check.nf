include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK(samplesheet)
            .csv
            .splitCsv(header:true, sep:',')
            .map{ it -> it + [id:it.sample]}                   //place everything except data_directory and expression_profile into meta map
            .map{ it -> [meta: it.findAll{ k, _v -> !(k in ['sample','data_directory','expression_profile']) },
                  data_directory: file(it.data_directory),
                  expression_profile: (it.expression_profile == null || it.expression_profile == "") ? [] : file(it.expression_profile)
                  ]}
            .set{ datasets }
        samplesheet_valid = SAMPLESHEET_CHECK.out.csv
    emit:
        datasets
        samplesheet_valid
}
