#!/usr/bin/env python3
# Cross-sample analysis template script
# Input is a space-delimited string with adatas

import os
import logging
import pandas as pd
import scipy as sp
import anndata as ad
import scanpy as sc
import numpy as np
import json
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def ligrec_from_adatas(adatas, type='ligrec_means', axis=1,
                            samples=None, spotlight=None):

    # extract stat from each adata
    ligrecs = [x.uns[type] for x in adatas if type in x.uns]

    # combine sample level
    combined = pd.concat(ligrecs, keys=samples, axis=axis)

    # move cell types to columns and perform the t-test
    res = combined.stack(future_stack=True)

    # if a cell type or pair has been specified, filter to that only
    if spotlight:
        sp_index = [np.all([y in x for y in spotlight ]) for x in res.index.get_level_values(-1)]
        res = res[sp_index]

    return res

def heatmap_report(adatas, spotlight=None, groups=None, show=100, filter=0.05, tool=None):
    samples = [a.obs['id'].unique()[0] for a in adatas]
    pvalues = None
    if tool =='squidpy_ligrec':
        ligrecs = ligrec_from_adatas(adatas, type='ligrec_means', spotlight=spotlight, samples=samples)
        pvalues = ligrec_from_adatas(adatas, type='ligrec_pvalues', spotlight=spotlight, samples=samples)
    elif tool =='spacemarkers_LRscores':
        ligrecs = ligrec_from_adatas(adatas, type='LRscores', spotlight=spotlight, samples=samples)
    elif tool =='Moran_I':
        #Moran's I is not per cell type pair, so no spotlighting (yet)
        ligrecs = ligrec_from_adatas(adatas, type='moranI', spotlight=None, samples=samples)
        #pick only the Moran's I value index
        ligrecs = ligrecs[ligrecs.index.get_level_values(-1) == 'I']

    if pvalues is not None:
        # filter sample ligrec pairs by p-value to reduce multiple testing
        sig_mask = (pvalues <= filter).any(axis=1)
        ligrecs = ligrecs[sig_mask]

    #join multiindex of squidpy ligrec into single index
    if(isinstance(ligrecs.index, pd.MultiIndex)):
        ligrecs.index = ['-'.join(map(str, idx)) for idx in ligrecs.index]

    if groups is not None:
        ligrecs_ttest = xsample_ttest(ligrecs, groups[0], groups[1])
        ligrecs_ttest_sig = ligrecs_ttest[ligrecs_ttest['pval_adj'] <= filter]

        #mark significant with a star in the index for plotting purposes
        ligrecs_ttest.index = [idx + '*' if idx in ligrecs_ttest_sig.index \
            else idx for idx in ligrecs_ttest.index]

        if(len(ligrecs_ttest_sig) == 0):
            log.warning(f"No significant differential interactions found for {tool}.")
        else:
            log.info(f"{len(ligrecs_ttest_sig)} significant differential \
                interactions found for {tool} with p_adj<={filter}.")

        res = ligrecs_ttest.sort_values('pval', ascending=True)
        memo = f"Top {show} differential interactions across samples. There are \
                {len(ligrecs_ttest_sig['pval_adj'])} significant differential \
                interactions (p_adj<={filter}) found between \
                {groups[0]} and {groups[1]}, marked with *."

    else:
        ligrecs['mean'] = ligrecs.mean(axis=1)
        res = ligrecs.sort_values('mean', ascending=False)
        memo = f"Top {show} mean interactions across samples shown as no groups \
                 were specified in the sample sheet."

    res_show = res[samples][:show]
    res_dict = res_show.to_dict()

    mqc_report = {
        "id": f"{tool}_interactions",
        "description": memo,
        "plot_type": "heatmap",
        "pconfig": {
            "ylab": "Sample",
            "ycats_samples": True,
            "xcats_samples": False,
            "xlab": "Interaction",
            "zlab": "Score",
            "title": "Interaction scores",
            "square": False
        },
        "data": res_dict
    }
    return mqc_report, res


def pseudobulk_adatas(adatas, vars=None, only_spatial=False):
    # collect pseudobulk expression for each adata
    # if only_spatial, select only genes with spatially variable expression
    pb_adatas = []
    for adata in adatas:
        adata = adata.to_memory()
        if only_spatial:
            adata = adata[:, adata.var['spatially_variable']==True]
        pb_adatas.append(sc.get.aggregate(adata, by=vars, func='sum'))
    pb_adata = ad.concat(pb_adatas, axis=0, join='outer')
    counts = pb_adata.layers['sum']
    # keep raw counts in X for pydeseq2, handling both dense and sparse matrices
    if sp.sparse.issparse(counts):
        counts = counts.tocsr().copy()
        # replace NaNs with zeros in the underlying sparse data array
        data = counts.data
        nan_mask = np.isnan(data)
        if np.any(nan_mask):
            data[nan_mask] = 0
    else:
        counts = np.asarray(counts).copy()
        counts[np.isnan(counts)] = 0
    pb_adata.X = counts
    del pb_adata.layers['sum']

    return pb_adata


def pydeseq_results(pb_adata, spotlight=None, cpus=1, design=None, contrast=None):
    # estimate size factors and dispersions with pydeseq2
    inference = DefaultInference(n_cpus=cpus)
    dds = DeseqDataSet(adata=pb_adata,
                       design=design,
                       inference=inference)
    dds.deseq2()
    de = DeseqStats(dds, contrast=contrast)

    return de

def de_report(de_dict, spotlight=None, filter=0.05, show=100, contrast=None, p='padj'):
    #plot logfoldchangs vs -log padj for all genes
    # multiqc want this format: "gene1": {"x": 1, "y": 2}
    # x is log2 fold change, y is -log10 adjusted p-value
    # show_p is mainly for testng purposes to allow plotting
    reports = []
    for de in de_dict.values():
        res_sig = de[de[p] <= filter]
        res_sig_show = res_sig.sort_values(p).head(show)
        res_sig_show.rename(columns={"log2FoldChange": "x", p: "y"}, inplace=True)
        res_sig_show['y'] = -np.log10(res_sig_show['y'])
        res_dict = res_sig_show.loc[:, ['x','y']].to_dict(orient='index')
        reports.append(res_dict)
    memo = f"DESeq2 results for {contrast[0]} variable with significant DE genes ({p}<={filter}) between {contrast[1]} and {contrast[2]}. Log2 fold change (X) and -log10 adjusted p-value (Y) shown."
    mqc_report = {
        "id": f"deseq2_{contrast[0]}",
        "description": memo,
        "plot_type": "scatter",
        "pconfig": {
            "xlab": "log2 fold change",
            "ylab": "-log10 adjusted p-value",
            "title": f"DESeq2 results for {contrast[0]}",
            "data_labels": list(de_dict.keys()),
            "ymin": 0,
            "y_lines": [
                {"value": -np.log10(filter), "color": "#ff0000", "width": 1, "dash": "dash", "label": "significance"}
            ]
        },
        "data": reports
    }
    return mqc_report



def neighbors_report(adatas, spotlight=None):
    if spotlight:
        cell_types = spotlight    
    else:
        cell_types = {}
        for adata in adatas:
            if "cell_type" in adata.obs:
                ct_counts = adata.obs['cell_type'].value_counts()
                for ct, count in ct_counts.items():
                    if ct in cell_types:
                        cell_types[ct] += count
                    else:
                        cell_types[ct] = count
        cell_types = cell_types.keys()
    
    reports = []
    for cell_type in cell_types:
        sample_dict = {}
        for adata in adatas:
            if "cell_type_interactions" in adata.uns:
                if cell_type not in adata.obs['cell_type'].cat.categories:
                    continue
                adata_cell_type_index = adata.obs['cell_type'].cat.categories.get_loc(cell_type)
                interactions = adata.uns["cell_type_interactions"][adata_cell_type_index]
                #construct dict of other cell types and interaction values
                other_cell_types = adata.obs['cell_type'].cat.categories.tolist()
                interaction_dict = {}
                for i, other_cell_type in enumerate(other_cell_types):
                    interaction_dict[other_cell_type] = interactions[i]
                sample_dict[adata.obs['id'].unique()[0]] = interaction_dict
        reports.append(sample_dict)

    mqc_report = {
        "id": "spatial_neighbors",
        "plot_type": "bar",
        "description": "Cell type immediate neighborhood across samples. \
        The rate of self-neighborhood indicates clustering of a cell type. \
        The rate of neighbors with other cell types (self omitted) indicates \
            how these clusters interact with each other.",
        "pconfig": {
            "title": "Cell type neighborhood across samples",
            "ylab": "Neighboring cell type share",
            "xlab": "Sample",
            "data_labels": list(cell_types)
        },
        "data": reports
    }
    return mqc_report

def centrality_reports(adatas, spotlight=None, scores=None, uns_key='cell_type_centrality_scores'):
    #separately for each score as multiqc does not support data_labels for heatmaps
    memos = {
        'closeness_centrality': "This graph measure reflects how close the group is to other nodes.",
        'average_clustering': "This graph measure shows the degree to which nodes cluster together.",
        'degree_centrality': "This graph measure is the fraction of non-group members connected to group members"
    }

    adatas_scores = [a.uns[uns_key].keys().tolist() for a in adatas if uns_key in a.uns]
    if not adatas_scores:
        return {}
    all_scores = set.union(*[set(s) for s in adatas_scores])
    if scores:
        scores = all_scores.intersection(set(scores))
    else:
        scores = all_scores
    reports = {}
    for score in scores:
        sample_dict = {}
        for adata in adatas:
            if uns_key not in adata.uns or score not in adata.uns[uns_key]:
                continue
            centrality_scores = adata.uns[uns_key][score]
            available_cell_types = adata.obs['cell_type'].cat.categories.tolist()
            if spotlight:
                cell_types = [cell_type for cell_type in available_cell_types if cell_type in spotlight]
            else:
                cell_types = available_cell_types
            centrality_dict = {}
            for i, cell_type in enumerate(available_cell_types):
                if cell_type in cell_types:
                    centrality_dict[cell_type] = centrality_scores.iloc[i]
            sample_dict[adata.obs['id'].unique()[0]] = centrality_dict

        mqc_report = {
            "id": f"{score}",
            "description": memos.get(score, f"Centrality score {score} across samples"),
            "plot_type": "table",
            "pconfig": {
                "ylab": "Sample",
                "xlab": "Cell type",
                "title": f"{score} centrality score across samples"
            },
            "data": sample_dict
        }
        reports[score] = mqc_report
    return reports


def plot_hist(df_pair, title=None, save=True):
    import matplotlib.pyplot as plt
    import seaborn as sns
    x = 8
    y = df_pair.shape[0] / 4

    plt.figure(figsize=(x, y))
    sns.heatmap(df_pair, cmap='coolwarm', fmt=".2f", linewidths=.5)
    plt.title(title)
    if save:
        plt.savefig(title.replace(" ", "_")+".png", dpi=300, bbox_inches='tight')

def xsample_ttest(df, group1, group2):
    res = df.copy()
    test = sp.stats.ttest_ind(res[group1], res[group2], axis=1)
    res['statistic'] = test.statistic
    res['pval'] = test.pvalue
    res.dropna(inplace=True)
    res['pval_adj'] = sp.stats.false_discovery_control(res['pval'], method='bh')
    res.sort_values('pval_adj', inplace=True)

    return res

def versions():
    with open ("versions.yml", "w") as f:
        f.write(f"{process}:\\n")
        f.write(f"    scipy: {sp.__version__}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
        f.write(f"    anndata: {ad.__version__}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    json: {json.__version__}\\n")


def get_vars(adatas, only=None):
    # find vars added by the staple pipeline
    vars = [x.uns['staple_meta_fields'].tolist() for x in adatas if 'staple_meta_fields' in x.uns]
    if (len(vars) == 0):
        log.warning("No added metadata fields found in any of the provided anndatas.")
        return None

    # drop id var from analysis vars
    for v in vars:
        if 'id' in v:
            v.remove('id')

    # find vars present in all samples
    common_vars = set(vars[0]).intersection(*vars[1:])
    log.info(f"Common added metadata fields across all samples: {common_vars}")
    if len(common_vars) == 0:
        log.warning("No common added metadata fields found across all provided anndatas.")

    res = common_vars

    if(only=='categorical'):
        # check that vars are categorical, same inside sample and differ across samples
        cats = {}
        for var in common_vars:
            is_categorical = all([isinstance(x.obs[var].dtype, pd.CategoricalDtype) for x in adatas])
            if not is_categorical:
                log.warning(f"Variable {var} is not categorical in all samples, skipping.")
                continue
            n_categories = [x.obs[var].nunique() for x in adatas]
            if len(set(n_categories)) != 1:
                log.warning(f"Variable {var} differs across individual samples .obs, skipping.")
                continue
            combined_categories = pd.concat([x.obs[[var,'id']] for x in adatas])
            levels = combined_categories.groupby(var, observed=True)['id'].unique()
            if (len(levels) != 2):
                log.warning(f"Can only contrast vars with 2 levels ({var} has {len(levels)}), skipping.")
                continue
            cats[var] = levels.to_dict()
        res = cats

    return res

def save_reports(mqc, res, name, mqc_reports="reports/mqc", reports="reports"):
    if mqc is not None:
        with open(f"{mqc_reports}/{name}_mqc.json","w") as f:
            json.dump(mqc, f, indent=4)
    if res is not None:
        res.to_csv(f"{reports}/{name}.csv")

if __name__ == '__main__':
    process       = "${task.process}"
    collected     = "${collected_items}"                 # these are whitespace separated paths to anndatas
    show          = int("${params.analyze.show_top}")    # how many top results to show
    cpus          = int("${task.cpus}")
    filter        = float("${params.analyze.filter}")    # p-value or adjusted p-value threshold for significance
    pb_vars       = "${params.analyze.pb_vars}"          # vars to use for pseudobulk grouping, comma-separated string
    only_spatial  = "${params.analyze.only_spatial}"\
                                    .lower() == 'true'   # only use spatially variable genes for ligrec and pseudobulk 

    adata_paths = collected.split(" ")
    adatas = [ad.read_h5ad(path, backed="r") for path in adata_paths]
    spotlight = "${params.analyze.spotlight}"  # this is a comma-separated string of cell type pairs to spotlight
    if ',' in spotlight:
        spotlight = spotlight.split(',')

    #place all mqc reports here
    reports_dir = "reports"
    os.makedirs(reports_dir, exist_ok=True)
    mqc_reports_dir = "reports/mqc"
    os.makedirs(mqc_reports_dir, exist_ok=True)

    # generate neighbors report
    try:
        log.info("Generating neighbors report.")
        neighbors = neighbors_report(adatas, spotlight=spotlight)
        with open(f"{mqc_reports_dir}/neighbors_mqc.json","w") as f:
            json.dump(neighbors, f, indent=4)
    except Exception as e:
        log.warning(f"Could not generate neighbors report: {e}")

    # generate centrality report - separately
    try:
        log.info("Generating centrality report.")
        centrality = centrality_reports(adatas, spotlight=spotlight)
        for score, report in centrality.items():
            with open(f"{mqc_reports_dir}/{score}_mqc.json","w") as f:
                json.dump(report, f, indent=4)
    except Exception as e:
        log.warning(f"Could not generate centrality report: {e}")

    # print versions now because later may be too late
    versions()

    # generate reports using added meta vars
    all_vars = get_vars(adatas)
    log.info(f"All variables added upstream for analysis: {all_vars}")
    cats = get_vars(adatas, only='categorical')
    log.info(f"Cat variables suitable for cross-sample contrasts: {cats}")
    
    # prepare pseudobulk adata for deseq2, with all variables as grouping variables
    if pb_vars:
        # split comma-separated variables, trim whitespace, and preserve order without duplicates
        raw_vars = pb_vars.split(",")
        pb_vars = []
        for v in raw_vars:
            v = v.strip()
            if v and v not in pb_vars:
                pb_vars.append(v)
        log.info(f"Using specified variables for pseudobulk grouping: {pb_vars}")
    else:
        # start from all_vars as a list, preserving order and removing duplicates
        pb_vars = list(dict.fromkeys(all_vars)) if all_vars is not None else []
        log.info(f"No specific variables for pseudobulk grouping specified, using all variables: {pb_vars}")
        # ensure required grouping variables are present while preserving order
        if 'cell_type' not in pb_vars:
            pb_vars.append('cell_type')  # add cell type as grouping variable
        if 'id' not in pb_vars:
            pb_vars.append('id')
    pb_adata = pseudobulk_adatas(adatas, vars=pb_vars, only_spatial=only_spatial)
    if(only_spatial):
        log.info(f"Used only spatially variable genes for pseudobulk, resulting in {pb_adata.shape[1]} genes.")
        pb_adata.write(f"{reports_dir}/svg_pseudobulk.h5ad")
    else:
        log.info(f"Used all genes for pseudobulk, resulting in {pb_adata.shape[1]} genes.")
        pb_adata.write(f"{reports_dir}/pseudobulk.h5ad")
    
    # make reports
    # if no cats found, just produce overall ligrec report
    if cats is None or len(cats) == 0:
        try:
            res_mqc, res = heatmap_report(adatas, spotlight=spotlight, show=show, tool='squidpy_ligrec', filter=filter)
            save_reports(res_mqc, res, "ligrec_overall")
        except Exception as e:
            log.warning(f"Could not generate overall ligand-receptor report: {e}")
        try:
            lrs_mqc, lrs = heatmap_report(adatas, spotlight=spotlight, show=show, tool='spacemarkers_LRscores', filter=filter)
            save_reports(lrs_mqc, lrs, "lrscores_overall")
        except Exception as e:
            log.warning(f"Could not generate overall LR scores report: {e}")
        try:
            moran_mqc, moran = heatmap_report(adatas, spotlight=spotlight, show=show, tool='Moran_I', filter=filter)
            save_reports(moran_mqc, moran, "moranI_overall")
        except Exception as e:
            log.warning(f"Could not generate overall Moran's I report: {e}")
    else:
        # for variables with 2 groups, perform appropriate tests
        for var in cats.keys():
            groups = [x for x in cats[var]]
            group1 = cats[var][groups[0]].tolist()
            group2 = cats[var][groups[1]].tolist()
            
            #squidpy ligrec
            try:
                res_mqc, res = heatmap_report(adatas, groups=[group1,group2],
                                              spotlight=spotlight, show=show, tool="squidpy_ligrec", filter=filter)
                save_reports(res_mqc, res, f"ligrec_diff_{var}_results")
            except Exception as e:
                log.warning(f"Could not generate ligand-receptor report for variable {var}: {e}")
            # SpaceMarkers LR scores
            try:
                lrs_mqc, lrs = heatmap_report(adatas, groups=[group1,group2],
                                              spotlight=spotlight, show=show, tool='spacemarkers_LRscores', filter=filter)
                save_reports(lrs_mqc, lrs, f"lrscores_diff_{var}_results")
            except Exception as e:
                log.warning(f"Could not generate LR scores report for variable {var}: {e}")

            # Moran's I
            try:
                moran_mqc, moran = heatmap_report(adatas, groups=[group1,group2],
                                                  spotlight=spotlight, show=show, tool='Moran_I', filter=filter)
                save_reports(moran_mqc, moran, f"Moran_I_diff_{var}_results")
            except Exception as e:
                log.warning(f"Could not generate Moran's I report for variable {var}: {e}")


            # deseq2 on pseudobulks split by cell type
            try:
                contrasts = [var]+[k for k in cats[var].keys()]
                stratify_var = 'cell_type' if 'cell_type' in pb_adata.obs else None
                strata = pb_adata.obs[stratify_var].unique() if stratify_var else ['unstratified']
                de_results = {}
                for ct in strata:
                    mask = pb_adata.obs[stratify_var] == ct if stratify_var else np.array([True]*pb_adata.shape[0])
                    adata = pb_adata[mask].copy()
                    if adata.shape[0] < 2:
                        log.warning(f"Not enough samples for DESeq2 analysis for cell type {ct} with variable {var}, skipping.")
                        continue
                    de = pydeseq_results(adata, 
                                        spotlight=spotlight, 
                                        cpus=cpus, 
                                        design=f"~{var}", 
                                        contrast= contrasts)
                    de.summary()
                    deseq_res = de.results_df
                    deseq_res = deseq_res[deseq_res['padj'] <= filter]
                    if(deseq_res.empty):
                        log.warning(f"No significant DE genes found for {ct} cell type with variable {var}.")
                    else:
                        log.info(f"Deseq2 results for {ct} cell type with variable {var}:{deseq_res.shape[0]} significant genes found.")
                        de_results[ct] = deseq_res
                        deseq_res.to_csv(f"{reports_dir}/deseq2_diff_{var}_results_{ct}.csv")
                mqc_report = de_report(de_results, spotlight=spotlight, filter=filter, show=show, contrast=contrasts)
                save_reports(mqc_report, None, f"deseq2_diff_{var}_results",
                                mqc_reports_dir, reports_dir)

            except Exception as e:
                log.warning(f"Could not perform DESeq2 analysis for variable {var}: {e}")

    #wrapup
    for adata in adatas:
        adata.file.close()
