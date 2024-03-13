## Introduction

**btc/spatial** is a bioinformatics pipeline for spatial data. Currently we only have support for spatial transcriptomics (10X Visium), but in the future support for other modalities
(Xenium, IMC, etc.) will be integrated into this pipeline as well.

## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,data_directory,n_cell_types,bleeding_correction,spatial_transcriptional_programs,expression_profile
"my_sample","/path/to/spaceranger/dir",5,false,false,"/path/to/profile.csv"
```

Each row represents a spatial transcriptomics sample and configuration parameters specific to that sample.

`sample`: A unique identifier for the sample.

`data_directory`: The path to the directory containing the output of the spaceranger pipeline.

`n_cell_types`: This parameter controls how many cell types to deconvolve into.

`bleeding_correction`: set to `true` if you want to enable bleeding correction for that sample.

`spatial_transcriptional_programs`: set to `true` if you want to enable spatial transcriptional programs for that sample.

`expression_profile`: Optional reference expression profiles (leave blank if not using) if you have known cell types in your Visium data (perhaps from matched scRNA). See more info about how to generate this file here: https://bayestme.readthedocs.io/en/latest/fine_mapping_workflow.html

Now, you can run the pipeline using:

```bash
nextflow run break-through-cancer/btc-spatial-pipelines \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

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
