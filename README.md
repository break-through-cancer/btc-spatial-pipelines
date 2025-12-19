## Introduction

**STAPLE** is a bioinformatics pipeline for 10X Visium and Visium HD spatial data. Being a simple preprocessing-deconvolution-interaction pipeline it features multiple tools for reference-free and reference-based deconvolution (cell typing) and cell-cell interaction analysis. In case of reference-based deconvolution, atlas may be fetched from a URL on an S3 location (e.g. CellxGene) or a local file, or matched-scRNA made available.


![image info](assets/btc-visium.svg)

## Usage

>[!note]
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. 

First, prepare a samplesheet with your input data that looks as follows: [samplesheet.csv](samplesheet.csv), where each row represents a spatial transcriptomics sample.

The default named columns are following:

  * `sample` required, a unique identifier for the sample

  * `data_directory` required, path to the 10x spaceranger `outs` directory

  * `expression_profile`: optional, (leave blank if not using), reference expression profiles from matched scRNA (different local path to each scRNA atlas) or a local scRNA atals (same path for each sample). If using a remotely stored atlas (such as CellXGene), rather pass `params.ref_scrna` and the atlas will be downloaded from the web.

Any extra columns will be treated as metadata and copied into the `meta` map and the resulting `.h5ad` object.


>[!IMPORTANT]
`export NXF_SINGULARITY_HOME_MOUNT=true` in order to allow matplotlib to write its logs (and avoid related error) if using singularity.

Run on segmented Visium HD with RCTD for cell typing and squipy ligand-receptor analysis (default) using remote atlas annotation. In case of CellXGene atlas, the cell type column is always `cell_type`, so it does not need to be explicitly specified.
```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \ 
   --outdir <OUTDIR> \
   --visium_hd 'cell_segmentations' \
   --reference_scrna https://datasets.cellxgene.cziscience.com/d1d90d18-2109-412f-8dc0-e014e8abb338.h5ad
```
In order to specify a different Visium HD resolution, change `visium_hd` param to an existing table name such as `square_008um` or `square_016um`. 

Run on Visium SD or HD with local reference atlas specified in the samplesheet:

```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --ref_scrna_type_col cell_type_column_name
```

Run on Visium SD with reference-free deconvolution:

```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --deconvolve.bayestme \
```

## Limitations

Not all the tools support all the formats! Use these guidelines to pick parameters in case the fully functioning defaults (RCTD + Squidpy) are not desired.

| tool/format | Visium SD | Visium HD | HD segmented | MultiQC |
| ----------- | --------- | --------- | ------------ | ---------- |
| RCTD | OK | OK | OK | OK |
| Squidpy | OK | OK | OK | OK |
| CoGAPS | OK | reduce gene N | reduce gene N | OK |
| BayesTME | OK | | | |
| SpaceMarkers SD | OK | | | OK |
| SpaceMarkers HD | | OK | | OK |

Spacemarkers for SD reports IMscores for gene names and undirected cell type interactions (cell_type1 near cell_type2 is no different to cell_type2 near cell_type1)

SpaceMarkers for HD reports IMscores for gene names in a directed fashion (cell_type1 near cell_type2 is different to cell_type2 near cell_type1) but also reports LRscores, which are the interaction scores between genes listed in a database that SpaceMarkers uses (CellChat) by default.

Furthermore, not all the tools are fully featured in the cross-sample analysis. So, Spacemakers are not yet integrated into the MultuQC module, as well as since BayesTME and CoGAPS are reference-free, the synthetic cell type outputs they produce does not match across samples.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
