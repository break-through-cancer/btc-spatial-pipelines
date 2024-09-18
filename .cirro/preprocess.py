#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np

SAMPLESHEET_REQUIRED_COLUMNS = ("sample", "data_directory", "n_cell_types", "bleeding_correction", "expression_profile")


def set_params_as_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    ds.logger.info([ds.params])
    samplesheet = pd.DataFrame([ds.params]).explode("data_directory")

    for colname in SAMPLESHEET_REQUIRED_COLUMNS:
        if colname not in samplesheet.columns:
            samplesheet[colname] = np.nan

    for colname in samplesheet.columns:
        if colname not in SAMPLESHEET_REQUIRED_COLUMNS:
            del samplesheet[colname]

    # Save to a file
    samplesheet.to_csv("samplesheet.csv", index=None)

    # Clear all nextflow params other than --outdir and --input
    # since the input samplesheet now contains all the information we need.
    to_remove = []
    for k in ds.params:
        if k != "outdir":
            to_remove.append(k)

    for k in to_remove:
        ds.remove_param(k)

    ds.add_param("input", "samplesheet.csv")

    # Log the samplesheet
    ds.logger.info(samplesheet.to_dict())


def main():
    ds = PreprocessDataset.from_running()

    set_params_as_samplesheet(ds)

    # log
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
