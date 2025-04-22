#!/usr/bin/env Rscript
library(reticulate)
library(spacexr)
library(Matrix)

reticulate::use_virtualenv('/app/venv')
ad <- reticulate::import('anndata')

#args:  adata_sc_path, adata_st_path, ncores, cell_type_col
adata_sc_path <- '${adata_sc}'
adata_st_path <- '${adata_st}'
ncores <- ${task.cpus}
outdir <- '${prefix}'
process <- '${task.process}'
cell_type_col <- '${params.type_col_scrna}'

#!/usr/bin/env Rscript
library(reticulate)
ad <- reticulate::import("anndata")

### prep spatial (query) object
#1. coords need to be present in anndata.obsm['spatial']
adata_st <- ad[["read_h5ad"]](adata_st_path)
spatial <- as.data.frame(adata_st[["obsm"]]["spatial"])
rownames(spatial) <- adata_st[["obs_names"]][["values"]]
#2. counts - genes need to be in rows, cells in columns
counts_st <- Matrix::t(as(adata_st[['X']], "CsparseMatrix"))
rownames(counts_st) <- as.character(adata_st[['var_names']][['values']])
colnames(counts_st) <- as.character(adata_st[['obs_names']][['values']])
query <- spacexr::SpatialRNA(coords=spatial, counts=counts_st)


### prep reference object
#1. counts - genes need to be in rows, cells in columns, opposite to anndata
adata_sc <- ad[["read_h5ad"]](adata_sc_path, backed = "r")

if (is.null(adata_sc[["raw"]])) {
  counts_sc <- Matrix::t(adata_sc[["X"]][])
  rownames(counts_sc) <- rownames(adata_sc[["var"]])
} else {
  counts_sc <- Matrix::t(adata_sc[["raw"]][["X"]][])
  rownames(counts_sc) <- rownames(adata_sc[["raw"]][["var"]])
}

colnames(counts_sc) <- as.character(adata_sc[["obs_names"]][["values"]])

#2. celltypes must be named factor
obs_cols <- colnames(adata_sc[['obs']])
message("found columns: ", paste(obs_cols, collapse = ', '))
if(!cell_type_col %in% obs_cols) {
  stop(sprintf('cell type column "%s" not found in adata_sc', cell_type_col))
}
celltypes_sc <- as.character(adata_sc[['obs']][[cell_type_col]])
names(celltypes_sc) <- adata_sc[['obs_names']][['values']]
celltypes_sc <- as.factor(celltypes_sc)

#3. drop non-matching genes from atlas object & convert to column orientation
counts_sc <- counts_sc[rownames(counts_sc) %in% rownames(counts_st), ]

#4. create reference object
ref <- spacexr::Reference(cell_types=celltypes_sc, counts=counts_sc)

message('run rctd') 
#this converts the reference object to a dense matrix
rctd <- spacexr::create.RCTD(spatialRNA=query, reference=ref, max_cores = ncores)
rctd_res <- spacexr::run.RCTD(rctd, doublet_mode = 'full')

#cell type deconvolution results
message('getting cell type deconvolution results')
cell_types <- spacexr::normalize_weights(rctd_res@results[['weights']])
message(sprintf('saving results to %s/', outdir))
dir.create(outdir, showWarnings = FALSE)
write.csv(as.matrix(cell_types), file=file.path(outdir, 'rctd_cell_types.csv'))

#versions
message("reading session info")
sinfo <- sessionInfo()
versions <- lapply(sinfo[["otherPkgs"]], function(x) {sprintf("  %s: %s",x[["Package"]],x[["Version"]])})
versions[['R']] <- sprintf("  R: %s\n",packageVersion("base"))
cat(paste0(process,":\n"), file="versions.yml")
cat(unlist(versions), file="versions.yml", append=TRUE, sep="\n")
