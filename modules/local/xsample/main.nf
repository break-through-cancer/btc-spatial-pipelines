process STAPLE_XSAMPLE {
    tag "cross-sample"
    container "ghcr.io/break-through-cancer/btc-containers/scverse@sha256:0471909d51c29a5a4cb391ac86f5cf58dad448da7f6862577d206ae8eb831216"


    input:
    path collected_items, stageAs: "?/*" // stage numbers the files as 0,1,2...

    output:
    path "reports/*.csv",               emit: reports,          optional: true //full csv reports
    path "reports/mqc/*mqc*",           emit: multiqc_files,    optional: true //multiqc reports for viewing and ai
    path "versions.yml",                emit: versions

    // collect multiple ligrec results into a single file

    script:
    template 'xsample.py'
}