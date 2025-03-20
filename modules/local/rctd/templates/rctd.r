#!/usr/bin/env Rscript
reticulate::use_virtualenv('/app/venv')

ad <- reticulate::import('anndata')

#args:  adata_sc_path, adata_st_path, ncores
args <- commandArgs(trailingOnly = TRUE)
adata_sc_path <- '${adata_sc}'
adata_st_path <- '${adata_st}'
ncores <- ${task.cpus}
outdir <- '${prefix}'

message('prepare rctd ref object')
#1. counts - genes need to be in rows, cells in columns
adata_sc <- ad[['read_h5ad']](adata_sc_path)
counts_sc <- Matrix::t(as(adata_sc[['X']], "CsparseMatrix"))
rownames(counts_sc) <- as.character(adata_sc[['var_names']][['values']])
colnames(counts_sc) <- as.character(adata_sc[['obs_names']][['values']])
#2. celltypes must be named factor
celltypes_sc <- as.character(adata_sc[['obs']][['cell_type']])
names(celltypes_sc) <- adata_sc[['obs_names']][['values']]
celltypes_sc <- as.factor(celltypes_sc)
#3. create reference object
ref <- spacexr::Reference(cell_types=celltypes_sc, counts=counts_sc)

message('prepare query object')
adata_st <- ad[['read_h5ad']](adata_st_path)
#1. coords need to be present in anndata.obsm['spatial']
spatial <- as.data.frame(adata_st[['obsm']]['spatial'])
rownames(spatial) <- adata_st[['obs_names']][['values']]
#2. counts - genes need to be in rows, cells in columns
counts_st <- Matrix::t(as(adata_st[['X']], "CsparseMatrix"))
rownames(counts_st) <- as.character(adata_st[['var_names']][['values']])
colnames(counts_st) <- as.character(adata_st[['obs_names']][['values']])
query <- spacexr::SpatialRNA(coords=spatial, counts=counts_st)

message('run rctd')
rctd <- spacexr::create.RCTD(spatialRNA=query, reference=ref, max_cores = ncores)
rctd_res <- spacexr::run.RCTD(rctd, doublet_mode = 'full')

#cell type deconvolution results
message('getting cell type deconvolution results')
cell_types <- spacexr::normalize_weights(rctd_res@results[['weights']])
message(sprintf('saving results to %s/', outdir))
dir.create(outdir, showWarnings = FALSE)
write.csv(as.matrix(cell_types), file=file.path(outdir, 'rctd_cell_types.csv'))
