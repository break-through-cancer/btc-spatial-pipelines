#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
from urllib.parse import urlparse
from pathlib import Path

SAMPLESHEET_REQUIRED_COLUMNS = ("sample", 
                                "data_directory", 
                                "n_cell_types", 
                                "bleeding_correction", 
                                "expression_profile",
                                "run_bayestme",
                                "run_cogaps",
                                "n_top_genes",
                                "run_spacemarkers",
                                "find_annotations",
                                "response"
                                )

# Helper function to check if a string is a URL
def is_url(string):
    try:
        result = urlparse(string)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


def set_params_as_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    ds.logger.info([ds.params])
    
    # If the reference_scrna is not a URL, we assume it is a file mask string
    # to look for in the data directory downstream
    if 'reference_scrna' in ds.params and not is_url(ds.params['reference_scrna']):
        ds.params['expression_profile'] = ds.params['reference_scrna']
    
    samplesheet = df_from_params(ds.params, ds)

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
        if k in SAMPLESHEET_REQUIRED_COLUMNS:
            to_remove.append(k)

    for k in to_remove:
        ds.remove_param(k)

    ds.add_param("input", "samplesheet.csv")

    # Log the samplesheet
    ds.logger.info(samplesheet.to_dict())


def df_from_params(params, ds):
    pipeline_param_names = [c for c in SAMPLESHEET_REQUIRED_COLUMNS]
    pipeline_params = { k: [params[k]] for k in pipeline_param_names if k in params.keys()}

    files = ds.files

    # Assumes samplesheet associates sample with a file in the sample's root directory
    files['data_directory'] = files['file'].apply(lambda x: str(Path(x).parent))
    files = files[['sample','data_directory']]

    data_params = pd.merge(ds.samplesheet,files,on='sample')
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
