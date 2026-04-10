process CREATE_TEST_ADATA {
    tag "create_test_adata"
    label 'process_single'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:0471909d51c29a5a4cb391ac86f5cf58dad448da7f6862577d206ae8eb831216'

    input:
        val num_adatas

    output:
        path "*adata.h5ad", emit: adata

    script:
    template 'create_adata.py'

}