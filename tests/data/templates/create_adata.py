#!/usr/bin/env python
import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
import scanpy as sc

#set random seed for reproducibility
np.random.seed(42)

num_adatas = "${num_adatas}"
with_metadata = "${with_metadata}".lower() == 'true'

def make_one_adata(n=25, m=1000, pct_mito=0.1, sample_id='sample', with_metadata=True):

    # create anndata object n obs by m vars with random counts
    adata = ad.AnnData(X=np.random.poisson(1, (n, m)), 
                    obs=pd.DataFrame(index=[f'obs_{i}' for i in range(n)]), 
                    var=pd.DataFrame(index=[f'var_{j}' for j in range(m)]))
    #make sparse 
    adata.X = adata.X.astype(np.float32)

    #mark a percentage of the genes as mitochondrial by name prefix
    mito_genes = np.random.choice(adata.var_names, size=int(pct_mito*m), replace=False)
    adata.var_names = ['MT-' + name if name in mito_genes else name for name in adata.var_names]

    # add spatial dimensions to represent a 5x5 grid
    adata.obsm['spatial'] = np.array([[i // 5, i % 5] for i in range(n)])
    library_id = 'spatial_data'
    adata.uns['spatial'] = {
        library_id: {
            'scalefactors': {
                'tissue_hires_scalef': 1.0,
                'spot_diameter_fullres': 1.0
            },
            'images': {}
        }
    }
    sq.gr.spatial_neighbors(adata, spatial_key='spatial')


    # add an assignment with all outside cells being stroma surrounding
    # tumor core and a couple randomly assigned other of type other
    adata.obs['cell_type'] = 'cancer'
    outside_indices = [i for i in range(n) if adata.obsm['spatial'][i, 0] in [0, 4] or adata.obsm['spatial'][i, 1] in [0, 4]]
    for i in outside_indices:
        adata.obs.at[f'obs_{i}', 'cell_type'] = 'stroma'
    adata.obs.at[f'obs_{np.random.randint(0,n-1)}', 'cell_type'] = 'other'  # random cell
    adata.obs.at[f'obs_{np.random.randint(0,n-1)}', 'cell_type'] = 'other'  # another random cell
    
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    
    # attach cell type interaction report
    sq.gr.interaction_matrix(adata, cluster_key='cell_type')
    
    # add a random cell type assignment of 3 cell types: tumor, stroma, other
    cell_types = ['tumor', 'stroma', 'other']
    adata.obs['cell_type'] = np.random.choice(cell_types, size=n)
    adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'], categories=cell_types)

    # add sample id for testing
    adata.obs['id'] = sample_id

    # add a random response variable for testing (no variation within adata)
    adata.obs['response'] = np.random.choice(['responder', 'non-responder'])

    # add a continuous variable for testing
    adata.obs['age'] = np.random.randint(20, 100)

    #simulate staple behavior of added metadata from samplesheet
    if with_metadata:
        adata.uns['staple_meta_fields'] = ['response', 'id', 'age']
    
    # make some genes differentially expressed in stroma vs tumor for testing
    stroma = adata.obs['cell_type'] == 'stroma'
    tumor = adata.obs['cell_type'] == 'tumor'
    adata.X[stroma, :50] += np.random.poisson(5, (stroma.sum(), 50))
    adata.X[tumor, :50] += np.random.poisson(1, (tumor.sum(), 50))

    # compute Moran's I (results stored in adata.uns['moranI']) for testing
    sq.gr.spatial_autocorr(adata, mode="moran", n_jobs=1)
    
    # mark some genes as spatially variable for testing
    adata.var['spatially_variable'] = adata.uns['moranI']['I'] > adata.uns['moranI']['I'].median()
    
    # add dummy ligand-receptor interactions for testing. ligand receptor data
    # is adata.uns['ligrec_means'], adata.uns['ligrec_pvalues']
    # with gene-pairs in row index, and celltype pairs in columns
    cell_types = adata.obs['cell_type'].cat.categories
    gene_pairs =  ['-'.join(['var_' + str(i), 'var_' + str(j)]) for i in range(10) for j in range(10)]
    celltype_pairs = ['-'.join([ct1, ct2]) for ct1 in cell_types for ct2 in cell_types]
    ligrec_means = pd.DataFrame(np.random.rand(len(gene_pairs), len(celltype_pairs)), index=gene_pairs, columns=celltype_pairs)
    ligrec_pvalues = pd.DataFrame(np.random.rand(len(gene_pairs), len(celltype_pairs)), index=gene_pairs, columns=celltype_pairs)
    adata.uns['ligrec_means'] = ligrec_means
    adata.uns['ligrec_pvalues'] = ligrec_pvalues

    # recompute interactions using the final cell_type categories so that
    # the stored interaction matrix and centrality scores stay in sync
    sq.gr.interaction_matrix(adata, cluster_key='cell_type')

    #add centrality measures
    sq.gr.centrality_scores(adata, cluster_key='cell_type')

    #set a small number of cell types to NA for testing NA handling
    na_indices = np.random.choice(adata.obs.shape[0], size=10, replace=False)
    adata.obs.loc[adata.obs.index[na_indices], 'cell_type'] = np.nan

    return adata


def make_many_adata(num_adatas=2, n=25, m=1000, pct_mito=0.1, with_metadata=True):
    adatas = []
    for i in range(num_adatas):
        adata = make_one_adata(n=n, m=m, pct_mito=pct_mito, sample_id=f'sample_{i}', with_metadata=with_metadata)
        adatas.append(adata)
    return adatas

if __name__ == "__main__":
    adatas = make_many_adata(num_adatas=int(num_adatas), n=25, m=1000, pct_mito=0.1, with_metadata=with_metadata)
    for i, adata in enumerate(adatas):
        adata.write_h5ad(f'{i}_adata.h5ad')