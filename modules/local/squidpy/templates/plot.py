#!/usr/bin/env python3
import anndata as ad
import squidpy as sq
import logging
import os
import numpy as np

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

adata_path = "${adata}"
sample = "${sample}"
out = "${prefix}"
process = "${task.process}"
cell_type = "cell_type"

os.makedirs(out, exist_ok=True)

log.info("loading {}".format(adata_path))
adata = ad.read_h5ad(adata_path)
log.info("adata is {}".format(adata))

#get most abundant cell type from bayestme
if cell_type not in adata.obs.columns:
    log.info("cell_type {} not found in adata.obs, calculating from bayestme_cell_type_counts")
    most_abundant = np.argmax(adata.obsm['bayestme_cell_type_counts'], axis=1)
    adata.obs[cell_type] = most_abundant.astype('str')

#squidpy insists on dir naming, not creating outdir as usually
main_dir = os.getcwd()
os.chdir(out)

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
                        color=[cell_type],
                        crop_coord=(min_x, min_y, max_x, max_y),
                        library_id=lib_id,
                        save="spatial_scatter.png",
                        title="{} Spatial Scatter Plot".format(sample),
                        dpi=300
                        )
else: #bayestme adata, no image
    sq.pl.spatial_scatter(adata, 
                        color=[cell_type],
                        shape=None,
                        na_color=(1,1,1,0),
                        crop_coord=(min_x, min_y, max_x, max_y),
                        save="spatial_scatter.png",
                        title="{} Spatial Scatter Plot".format(sample),
                        dpi=300
                        )


# Plot the interaction matrix
sq.gr.spatial_neighbors(adata)
sq.gr.interaction_matrix(adata, cluster_key="cell_type")
sq.pl.interaction_matrix(adata,
                        cluster_key="cell_type",
                        save="interaction_matrix.png",
                        title="{} Interaction Matrix".format(sample))

#Plot the co-occurence, needs not NaN clusters to run
nona_adata = adata[~adata.obs[cell_type].isna()]
clusters = nona_adata.obs[cell_type].unique()
sq.gr.spatial_neighbors(nona_adata)
sq.gr.co_occurrence(nona_adata, cluster_key=cell_type)

for c in clusters:
    sq.pl.co_occurrence(nona_adata,
                        cluster_key=cell_type,
                        clusters=c,
                        save="co_occurrence_{}.png".format(c),
                        dpi=300
                        )


#Plot nhood enrichment
sq.gr.nhood_enrichment(nona_adata, cluster_key=cell_type)
sq.pl.nhood_enrichment(nona_adata,
                        cluster_key=cell_type,
                        save="nhood_enrichment.png",
                        title="{} Nhood Enrichment".format(sample),
                        dpi=300
                        )

#Plot centrality scores
sq.gr.centrality_scores(adata, cluster_key=cell_type)
sq.pl.centrality_scores(adata,
                        score="degree_centrality",
                        cluster_key=cell_type,
                        save="centrality_scores.png",
                        dpi=300
                        )

#save the object with newly calculated attributes
adata.write_h5ad("squidpy.h5ad", compression="gzip")

#versions
os.chdir(main_dir)
with open("versions.yml", "w") as f:
    f.write("{}:\\n".format(process))
    f.write("    squidpy: {}\\n".format(sq.__version__))
    f.write("    anndata: {}\\n".format(ad.__version__))
    f.write("    numpy: {}\\n".format(np.__version__))
    f.write("    python: {}\\n".format(os.sys.version.split()[0]))