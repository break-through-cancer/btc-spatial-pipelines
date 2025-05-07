#!/usr/bin/env python3
import anndata as ad
import squidpy as sq
import logging
import os

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


adata_path = "${rctd_adata}"
sample = "${prefix}"

#squidpy insists on dir naming, not creating outdir as usually
log.info("loading {}".format(adata_path))
adata = ad.read_h5ad(adata_path, backed="r")
log.info("adata is {}".format(adata))


# Extract spatial coordinates from the AnnData object
spatial_coords = adata.obsm['spatial']

# Get the minimum and maximum x and y coordinates
min_x, min_y = spatial_coords[:, 0].min(), spatial_coords[:, 1].min()
max_x, max_y = spatial_coords[:, 0].max(), spatial_coords[:, 1].max()

# Plot the spatial scatter plot
#format selector
if 'spatial' in adata.uns:
    # Get the library id for plotting
    try:
        lib_id = [k for k in adata.uns["spatial"].keys() if "hires" in k][0]
    except IndexError:
        log.warning("no hires library found, using the first one")
        try:
            lib_id = adata.uns["spatial"].keys()[0]
        except:
            log.error("no library data found in adata.uns['spatial']")
            raise
    #plot
    sq.pl.spatial_scatter(adata, 
                        color=["cell_type"],
                        crop_coord=(min_x, min_y, max_x, max_y),
                        shape=None,
                        save="{}_spatial_scatter.png".format(sample),
                        title="{} Spatial Scatter Plot".format(sample),
                        dpi=300
                        )
else: #bayestme adata
    sq.pl.spatial_scatter(adata, 
                        color=["cell_type"],
                        shape=None,
                        crop_coord=(min_x, min_y, max_x, max_y),
                        save="{}_spatial_scatter.png".format(sample),
                        title="{} Spatial Scatter Plot".format(sample),
                        dpi=300
                        )


# Plot the interaction matrix
sq.gr.spatial_neighbors(adata)
sq.gr.interaction_matrix(adata, cluster_key="cell_type")
sq.pl.interaction_matrix(adata,
                        cluster_key="cell_type",
                        save="{}_interaction_matrix.png".format(sample),
                        title="{} Interaction Matrix".format(sample))
