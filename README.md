## Introduction

**btc/spatial** is a bioinformatics pipeline for 10X Visium and Visium HD spatial data. Being a simple preprocessing-deconvolution-interaction pipeline it features multiple tools, including reference-based deconvolution. In case of reference-based deconvolution, atlas may be specified as an URL on an S3 location (e.g. CellxGene), or matched-scRNAmade made available in the sample folder.
![image info](assets/btc-visium.svg)

## Usage

>[!note]
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. 

First, prepare a samplesheet with your input data that looks as follows: [samplesheet.csv](samplesheet.csv):

Each row represents a spatial transcriptomics sample and configuration parameters specific to that sample.

`sample`: A unique identifier for the sample.

`data_directory`: The path to the directory containing the output of the spaceranger pipeline.

`expression_profile`: Optional reference expression profiles (leave blank if not using) if you have known cell types in your Visium data from matched scRNA (different local path to each scRNA atlas) or a local scRNA atals (same path for each sample). If using a remotely stored atlas (such as CellXGene), rather pass `params.ref_scrna` and the atlas will be downloaded from the web.


>[!IMPORTANT]
`export NXF_SINGULARITY_HOME_MOUNT=true` in order to allow matplotlib to write its logs (and avoid related error) if using singularity.

Simple Visium SD run with reference-free deconvolution:

```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --deconvolve.bayestme \
```

Run on Visium SD with local reference atlas specified in the samplesheet:

```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --deconvolve.rctd
   --ref_scrna_type_col cell_type_column_name
```


Run with Visium HD support and remote reference atlas annotation
```
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \ 
   --outdir out \
   --visium_hd 'square_016um' \
   --reference_scrna https://datasets.cellxgene.cziscience.com/d1d90d18-2109-412f-8dc0-e014e8abb338.h5ad
```

>[!WARNING]
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
