## Introduction

**btc/spatial** is a bioinformatics pipeline for 10X Visium and Visium HD spatial data. Being a simple preprocessing-deconvolution-interaction pipeline it features multiple tools, including reference-based deconvolution. In case of reference-based deconvolution, atlas may be specified as an URL on an S3 location (e.g. CellxGene), or matched-scRNAmade made available in the sample folder.
![image info](assets/btc-visium.svg)

## Usage

>[!note]
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.


First, prepare a samplesheet with your input data that looks as follows: [samplesheet.csv](samplesheet.csv):


Each row represents a spatial transcriptomics sample and configuration parameters specific to that sample.

`sample`: A unique identifier for the sample.

`data_directory`: The path to the directory containing the output of the spaceranger pipeline.

`n_cell_types`: This parameter controls how many cell types to deconvolve into. If you pass `expression_profile` this value will be ignored, we will determine number of cell types from the matched scRNA data.

`bleeding_correction`: set to `true` if you want to enable bleeding correction for that sample.

`expression_profile`: Optional reference expression profiles (leave blank if not using) if you have known cell types in your Visium data (perhaps from matched scRNA). See more info about how to generate this file here: https://bayestme.readthedocs.io/en/latest/fine_mapping_workflow.html

`run_bayestme`: boolean, whether to run BayesTME deconvolution.

`run_cogaps`: boolean, whether to run BayesTME deconvolution.

`run_spacemarkers`: boolean, whether to run SpaceMarkers

`find_annotations`: boolean, if `true`, CODA annotations will be sought by `*tissue_positions_cellular_compositions.csv` string and the data will be fed to SpaceMarkers. This has potential to run any external annotation from a `csv`.


>[!IMPORTANT]
`export NXF_SINGULARITY_HOME_MOUNT=true` in order to allow matplotlib to write its logs (and avoid related error) if using singularity.

Simple run without reference deconvolution:

```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```


Run with Visium HD support and reference atlas annotation
```
nextflow run break-through-cancer/btc-spatial-pipelines 
   -profile docker
   --input samplesheet.csv 
   --outdir out
   --hd 'square_016um' 
   --reference_scrna https://datasets.cellxgene.cziscience.com/d1d90d18-2109-412f-8dc0-e014e8abb338.h5ad
   --type_col_scrna Clusters
   -resume
```

>[!WARNING]
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

```
Zhang H, Hunter MV, Chou J, Quinn JF, Zhou M, White RM, Tansey W. BayesTME: An end-to-end method for multiscale spatial transcriptional profiling of the tissue microenvironment. Cell Syst. 2023 Jul 19;14(7):605-619.e7. doi: 10.1016/j.cels.2023.06.003. PMID: 37473731; PMCID: PMC10368078.
```

```
Deshpande A, Loth M, Sidiropoulos DN, Zhang S, Yuan L, Bell ATF, Zhu Q, Ho WJ, Santa-Maria C, Gilkes DM, Williams SR, Uytingco CR, Chew J, Hartnett A, Bent ZW, Favorov AV, Popel AS, Yarchoan M, Kiemen A, Wu PH, Fujikura K, Wirtz D, Wood LD, Zheng L, Jaffee EM, Anders RA, Danilova L, Stein-O'Brien G, Kagohara LT, Fertig EJ. Uncovering the spatial landscape of molecular interactions within the tumor microenvironment through latent spaces. Cell Syst. 2023 Apr 19;14(4):285-301.e4. doi: 10.1016/j.cels.2023.03.004. Erratum in: Cell Syst. 2023 Aug 16;14(8):722. PMID: 37080163; PMCID: PMC10236356.
```

```
Sherman TD, Gao T, Fertig EJ. CoGAPS 3: Bayesian non-negative matrix factorization for single-cell analysis with asynchronous updates and sparse data structures. BMC Bioinformatics. 2020 Oct 14;21(1):453. doi: 10.1186/s12859-020-03796-9. PMID: 33054706; PMCID: PMC7556974.
```

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
