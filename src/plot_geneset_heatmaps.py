import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

#Generates a geneset*geneset intersection proportion plot
#denominator is based on index, not column
def gs_gs_plot(genes, geneset_df, output_prefix='eGeneEnrich_plots'):
    """
    genes: genes to look at for comparisons
    geneset_df: a pandas DataFrame containing a row for every gene in geneset
    geneset column should be labeled 'geneset', gene column should be labeled 
    'ensembl_gene_id'
    
    output_prefix: prefix for file to be written. Should contain path
    chunk_size: memory manager for large number of genes. 
    """
    #get size for each geneset
    geneset_sizes = geneset_df.groupby('gene_set')['gene_set'].size()

    #dict for each geneset
    gs_dict = pd.Series(geneset_df.ensembl_gene_id)
    gs_dict.index = geneset_df.gene_set
    gs_dict = gs_dict.groupby(level=0).agg(set).to_dict()

    #Get combination values
    gs_series = pd.Series({a:b for a,b in zip(itertools.combinations_with_replacement(gs_dict.keys(),2),
                 [len(c&d) for c,d in itertools.combinations_with_replacement(gs_dict.values(),2)])})

    #unstack series
    gs_df = gs_series.unstack()
    proper_index = gs_series.index.get_level_values(0).unique()
    gs_df = gs_df.reindex(proper_index)[proper_index.values].copy()

    #reflect matrix
    i_lower = np.tril_indices(np.array(len(gs_df.columns)), -1)
    gs_matrix = gs_df.values
    gs_matrix[i_lower] = gs_matrix.T[i_lower]

    #Turn each value into the proportion of the index gs rather than a count
    for column in gs_df.columns:
        gs_df[column] = ((gs_df[column])/
                               (float(geneset_sizes[column])))

    gs_df = gs_df.T.copy()
    #Generate plots
    ax = sns.clustermap(gs_df,cmap='YlGnBu', rasterized=True, vmin=0, vmax=1,
                       figsize=(15, 15), annot=False)
    ax.savefig(output_prefix+'_gs_gs_intersection.pdf')

#Generates a geneset*geneset intersection count plot
#denominator is based on index, not column
def gs_gs_count_plot(genes, geneset_df, output_prefix='eGeneEnrich_plots'):
    """
    genes: genes to look at for comparisons
    geneset_df: a pandas DataFrame containing a row for every gene in geneset
    geneset column should be labeled 'geneset', gene column should be labeled 
    'ensembl_gene_id'
    
    output_prefix: prefix for file to be written. Should contain path
    chunk_size: memory manager for large number of genes. 
    """
    #get size for each geneset
    geneset_sizes = geneset_df.groupby('gene_set')['gene_set'].size()

    #dict for each geneset
    gs_dict = pd.Series(geneset_df.ensembl_gene_id)
    gs_dict.index = geneset_df.gene_set
    gs_dict = gs_dict.groupby(level=0).agg(set).to_dict()

    #Get combination values
    gs_series = pd.Series({a:b for a,b in zip(itertools.combinations_with_replacement(gs_dict.keys(),2),
                 [len(c&d) for c,d in itertools.combinations_with_replacement(gs_dict.values(),2)])})

    #unstack series
    gs_df = gs_series.unstack()
    proper_index = gs_series.index.get_level_values(0).unique()
    gs_df = gs_df.reindex(proper_index)[proper_index.values].copy()

    #reflect matrix
    i_lower = np.tril_indices(np.array(len(gs_df.columns)), -1)
    gs_matrix = gs_df.values
    gs_matrix[i_lower] = gs_matrix.T[i_lower]

    gs_df = gs_df.T.copy()
    #Generate plots
    ax = sns.clustermap(gs_df,cmap='YlGnBu', rasterized=True, annot=False,
                        robust=True, figsize=(15,15))
    hm = ax.ax_heatmap.get_position()

    ax.savefig(output_prefix+'_gs_gs_intersection.pdf')


def gs_gene_plot(genes, sig_egenes, geneset_df, output_prefix='eGeneEnrich_plots'):
    """
    genes: genes to look at for comparison
    sig_egenes: genes considered significant
    geneset_df: a pandas DataFrame containing a row for every gene in geneset
    geneset column should be labeled 'geneset', gene column should be labeled 
    'ensembl_gene_id'

    output_prefix: prefix for file to be written. Should contain path
    """
    genes = [gene for gene in genes if gene in list(geneset_df.ensembl_gene_id)]
    sig_egenes = [gene for gene in sig_egenes if gene in list(geneset_df.ensembl_gene_id)]

    #dict for genes in each geneset
    gs_dict = geneset_df.groupby('gene_set')['ensembl_gene_id'].apply(set).to_dict()

    #table for geneset x gene
    plot_df = pd.DataFrame({}, index = genes)

    for key, val in gs_dict.items():
        plot_df[key] = pd.Series(np.ones(len(val)), index=list(val))

    plot_df = plot_df.fillna(0).astype(int)
    plot_df['significant'] = 0
    plot_df.loc[plot_df.index.isin(sig_egenes), 'significant'] = 2
    for col in plot_df.columns:
        if col != 'significant':
            plot_df[col] = plot_df[col] + plot_df['significant']
    plot_df.drop('significant', axis=1, inplace=True)
    #replace 2s with 0. and 3s with 2 Values will be 0, 1, 2
    plot_df = plot_df.replace(2,0).copy()
    plot_df = plot_df.replace(3,2).copy()
    #generate plots
#    ax = sns.clustermap(plot_df,cmap='YlGnBu', rasterized=True, vmin=0, vmax=2)
#    ax.savefig(output_prefix+'gs_gene_intersection.pdf')
    ax = sns.clustermap(plot_df,cmap='YlGnBu', rasterized=True, vmin=0, vmax=2)
    ax.savefig(output_prefix+'_gs_gene_intersection.pdf')


def gs_gene_all_labels_plot(genes, sig_egenes, geneset_df, output_prefix='eGeneEnrich_plots'):
    """
    genes: genes to look at for comparison
    sig_egenes: genes considered significant
    geneset_df: a pandas DataFrame containing a row for every gene in geneset
    geneset column should be labeled 'geneset', gene column should be labeled 
    'ensembl_gene_id'

    output_prefix: prefix for file to be written. Should contain path
    """
    genes = [gene for gene in genes if gene in list(geneset_df.ensembl_gene_id)]
    sig_egenes = [gene for gene in sig_egenes if gene in list(geneset_df.ensembl_gene_id)]

    #dict for genes in each geneset
    gs_dict = geneset_df.groupby('gene_set')['ensembl_gene_id'].apply(set).to_dict()

    #table for geneset x gene
    plot_df = pd.DataFrame({}, index = genes)

    for key, val in gs_dict.items():
        plot_df[key] = pd.Series(np.ones(len(val)), index=list(val))

    plot_df = plot_df.fillna(0).astype(int)
    plot_df['significant'] = 0
    plot_df.loc[plot_df.index.isin(sig_egenes), 'significant'] = 2
    for col in plot_df.columns:
        if col != 'significant':
            plot_df[col] = plot_df[col] + plot_df['significant']
    plot_df.drop('significant', axis=1, inplace=True)
    #replace 2s with 0. and 3s with 2 Values will be 0, 1, 2
    plot_df = plot_df.replace(2,0).copy()
    plot_df = plot_df.replace(3,2).copy()

    ####Plotting parameters to make all genes appear
    fontsize_pt = 6
    dpi = 90
    height = len(plot_df.index)
    width = len(plot_df.columns)

    #Height in pixels and inches
    matrix_height_pt = fontsize_pt * height
    matrix_height_in = matrix_height_pt / dpi
    
    matrix_width_pt = fontsize_pt * width
    matrix_width_in = matrix_width_pt / dpi

    cm = sns.clustermap(plot_df,cmap='YlGnBu', rasterized=True, vmin=0, vmax=2,
                        yticklabels=True, figsize=(matrix_width_in*5, matrix_height_in))
    hm = cm.ax_heatmap.get_position()
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
    plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), fontsize=6)
    cm.cax.yaxis.label.set_size(6)
    
    cm.ax_heatmap.set_position([0.1, 0.01, 0.9, 0.99])
    col = cm.ax_col_dendrogram.get_position()
    cm.ax_col_dendrogram.set_position([0.1, 1, 0.9, 0.01])    
    row = cm.ax_row_dendrogram.get_position()
    cm.ax_row_dendrogram.set_position([0.05, 0.01, 0.05, 0.99])
    cbar = cm.cax.get_position()
    cm.cax.set_position([0, 0.95, 0.01, 0.005])
    cm.savefig(output_prefix+'_gs_gene_intersection.pdf')

    
    
