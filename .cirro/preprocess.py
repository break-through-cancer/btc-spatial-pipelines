#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
import validators

SAMPLESHEET_REQUIRED_COLUMNS = ("sample", 
                                "data_directory", 
                                "n_cell_types", 
                                "bleeding_correction", 
                                "expression_profile",
                                "run_bayestme",
                                "run_cogaps",
                                "n_top_genes",
                                "spatial_transcriptional_programs",
                                "run_spacemarkers",
                                "find_annotations"
                                )


def set_params_as_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    ds.logger.info([ds.params])
    
    # If the reference_scrna is a URL, we assume it is a file mask string
    # to look for in the data directory downstream
    if not validators.url(ds.params['reference_scrna']):
        ds.params['expression_profile'] = ds.params['reference_scrna']
    
    samplesheet = df_from_params(ds.params)

    for colname in SAMPLESHEET_REQUIRED_COLUMNS:
        if colname not in samplesheet.columns:
            samplesheet[colname] = np.nan

    for colname in samplesheet.columns:
        if colname not in SAMPLESHEET_REQUIRED_COLUMNS:
            del samplesheet[colname]

    # Save to a file
    samplesheet.to_csv("samplesheet.csv", index=None)

    # Clear params that we wrote to the samplesheet
    # cleared params will not overload the nextflow.params
    to_remove = []
    for k in ds.params:
        if k not in ["outdir","reference_scrna","input"]:
            to_remove.append(k)

    for k in to_remove:
        ds.remove_param(k)

    ds.add_param("input", "samplesheet.csv")

    # Log the samplesheet
    ds.logger.info(samplesheet.to_dict())


def df_from_params(params):
    pipeline_param_names = [c for c in SAMPLESHEET_REQUIRED_COLUMNS]
    pipeline_params = { k: [params[k]] for k in pipeline_param_names if k in params.keys()}

    data_params = pd.DataFrame({
        'sample':[x['name'] for x in params['cirro_input']],
        'data_directory': [x['s3']+'/data' for x in params['cirro_input']]
        })
    
    samplesheet = data_params.join(pd.DataFrame(pipeline_params), how='cross')

    return samplesheet

def main():
    ds = PreprocessDataset.from_running()

    set_params_as_samplesheet(ds)

    # log
    ds.logger.info(ds.params)
    print(ds.params)


if __name__ == "__main__":
    main()
