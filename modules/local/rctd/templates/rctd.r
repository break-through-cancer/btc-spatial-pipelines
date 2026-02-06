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

set.seed(${params.seed})

cell_type_col <- '${params.ref_scrna_type_col}'
n_top_genes <- as.numeric('${params.deconvolve.n_top_genes}')
doublet_mode <- '${params.doublet_mode}'

### prep spatial (query) object
#1. coords need to be present in anndata.obsm['spatial']
adata_st <- ad[["read_h5ad"]](adata_st_path)
spatial <- as.data.frame(adata_st[["obsm"]]["spatial"])
rownames(spatial) <- adata_st[["obs_names"]][["values"]]

#2. counts - genes need to be in rows, cells in columns
counts_st <- Matrix::t(as(adata_st[['X']], "CsparseMatrix"))
rownames(counts_st) <- as.character(adata_st[['var_names']][['values']])
colnames(counts_st) <- as.character(adata_st[['obs_names']][['values']])

#3. select top variable genes with less RAM footprint (as opposed to apply())
gene_vars <- pbapply::pbsapply(rownames(counts_st), function(x){var(counts_st[x,])})
gene_vars <- gene_vars[gene_vars > 0] #remove genes with zero variance
top_genes <- sort(gene_vars, decreasing = TRUE)[1:min(n_top_genes, length(gene_vars))]
counts_st <- counts_st[rownames(counts_st) %in% names(top_genes), ]

#4. prep query object
query <- spacexr::SpatialRNA(coords=spatial, counts=counts_st)


### prep reference object
#1. counts - genes need to be in rows, cells in columns, opposite to anndata
adata_sc <- ad[["read_h5ad"]](adata_sc_path, backed = "r")
adata_sc[["obs_names_make_unique"]]()

if (is.null(adata_sc[["raw"]])) {
  message('no raw data, using X from adata_sc')
  gene_names <-rownames(adata_sc[["var"]])
  select_genes <- which(gene_names %in% names(top_genes))
  counts_sc <- Matrix::t(adata_sc[["X"]][, select_genes - 1]) # -1 for 0-based index
  rownames(counts_sc) <- gene_names[select_genes]
} else {
  message('using raw.X from adata_sc')
  gene_names <- rownames(adata_sc[["raw"]][["var"]])
  select_genes <- which(gene_names %in% names(top_genes))
  counts_sc <- Matrix::t(adata_sc[["raw"]][["X"]][, select_genes - 1]) # -1 for 0-based index
  rownames(counts_sc) <- gene_names[select_genes]
}

#room for RAM improvement: read in at most 10K cells by index as RCTD limit
colnames(counts_sc) <- as.character(adata_sc[["obs_names"]][["values"]])

#2. celltypes must be named factor
obs_cols <- colnames(adata_sc[['obs']])
message("found columns: ", paste(obs_cols, collapse = ', '))
if(!cell_type_col %in% obs_cols) {
  stop(sprintf('cell type column "%s" not found in adata_sc', cell_type_col))
}
#read and clean up cell_types
celltypes_sc <- as.character(adata_sc[['obs']][[cell_type_col]])
celltypes_sc <- iconv(celltypes_sc,to="ASCII//TRANSLIT")                             #rm diacritics
celltypes_sc <- gsub(pattern = "[^[:alnum:] ]", replacement = " ", x = celltypes_sc) #rm punctuation
celltypes_sc <- gsub(pattern = " +", replacement = " ", x = celltypes_sc)            #rm extra spaces
names(celltypes_sc) <- adata_sc[['obs_names']][['values']]
celltypes_sc <- as.factor(celltypes_sc)

#3. drop rare cells & convert to column orientation
cell_stats <- table(celltypes_sc)
all_cells <- names(cell_stats)
rare_cells <- names(cell_stats[cell_stats < 25])
if (length(rare_cells) > 0) {
  message(sprintf("dropping %d rare cell types: %s", 
                  length(rare_cells), paste(rare_cells, collapse = ", ")))
} else {
  message("no rare cell types found, proceeding with all cells")
}
counts_sc <- counts_sc[, celltypes_sc %in% setdiff(all_cells, rare_cells)]
celltypes_sc <- celltypes_sc[celltypes_sc %in% setdiff(all_cells, rare_cells)]
celltypes_sc <- droplevels(celltypes_sc)
counts_sc <- as(counts_sc, "CsparseMatrix")

#4. create reference object and cleanup
ref <- spacexr::Reference(cell_types=celltypes_sc, counts=counts_sc)
counts_sc <- NULL
gc()

message('run rctd')
rctd_res <- tryCatch({
  spacexr::run.RCTD(spacexr::create.RCTD(spatialRNA=query, reference=ref, max_cores = ncores),
                    doublet_mode = doublet_mode)
}, error = function(e) {
  message('RCTD threw error: "',e[["message"]],'"')
  if(grep(pattern="UMI_min_sigma", x=e[["message"]])){
    message('RCTD error caught, retrying with UMI_min_sigma=1')
    spacexr::run.RCTD(spacexr::create.RCTD(spatialRNA=query, reference=ref, max_cores = ncores, UMI_min_sigma = 1),
                      doublet_mode = doublet_mode)
  } else {
    stop("Could not catch the error")
  }
})

#cell type deconvolution results
message('getting cell type deconvolution results')
cell_types <- spacexr::normalize_weights(rctd_res@results[['weights']])
message(sprintf('saving results to %s/', outdir))
dir.create(outdir, showWarnings = FALSE)
#cell types csv
write.csv(as.matrix(cell_types), file=file.path(outdir, 'rctd_cell_types.csv'))
#raw object
saveRDS(rctd_res, file=file.path(outdir, 'rctd.rds'))

#versions
message("reading session info")
sinfo <- sessionInfo()
versions <- lapply(sinfo[["otherPkgs"]], function(x) {sprintf("    %s: %s",x[["Package"]],x[["Version"]])})
versions[['R']] <- sprintf("    R: %s
",packageVersion("base"))
cat(paste0(process,":
"), file="versions.yml")
cat(unlist(versions), file="versions.yml", append=TRUE, sep="
")
