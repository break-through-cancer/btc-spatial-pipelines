# btc/spatial: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - 2025-07-11

### `Added`
- optionally look for external (e.g. CODA) annotation in the input to run spaceMarkers on
- spatial overlap and interaction plots for spaceMarkers
- RCTD deconvolution with single cell atlas reference
- Visium HD support (RCTD)
- Squidpy module for plotting deconvolution results and spatial analysis

### `Fixed`
- refactor: separate subworkflow for data loading
- refactor: separate subworkflow for BayesTME to avoid it clashing with VisiumHD pipeline

### `Dependencies`

### `Deprecated`
- spatial transcriptional programs from the BayesTME package
- BAYESTME_LOAD_SPACERANGER deprecated to allow Visium HD (and possible more) input
- do not reuse BayesTME preprocessing for CoGAPS


## v1.1.0 - 2024-12-16

### `Added`
- multiqc report for Spaceranger
- multiqc report for SpaceMarkers
- standalone SpaceMarkers module to run after any/both deconvolutions
- CoGAPS deconvolution module
- major analysis steps can now be turned on/off
- number of top genes to keep as parameter for deconvolution
- radically run CoGAPS after BayesTME preprocessing

### `Fixed`
- Cirro not being able to run on multiple samples
- Cirro showing 'sample' instead of the proper sample name
- SpaceMarkers returning no interactions in BayesTME results
- SpaceMarkers not running multithreaded
- spatial transcriptional programs running if enable = False


## v1.0dev - [date]

Initial release of btc/spatial, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
