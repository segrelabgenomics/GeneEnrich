#!/usr/bin/env python3

"""
GeneEnrich was written under the supervision of Ayellet Segre
Department of Ophthalmology and Ocular Genomics Institute, Massachusetts Eye and Ear, Harvard Medical School

Authors: John Rouhana, Andrew Hamel

Code written in Python version 3.6

GeneEnrich assesses enrichment of a list of genes of interst in prespecified gene sets or pathways
with the option to correct for tissue expression levels.

"""

import pandas as pd
import numpy as np
import time
import sys
import argparse
from functools import partial
from decimal import Decimal
from tools import str2bool, padjust_bh, df_to_delim, get_bins
from parse_genes_list import *
import scipy.stats
import re
import gzip
from plot_geneset_heatmaps import gs_gs_plot, gs_gene_plot
import datetime


def parse_args():
    """
    arguments for GeneEnrich
    """
    parser = argparse.ArgumentParser(description='Run GeneEnrich')
    parser.add_argument('--min_permutations', type=int, default=1000,
                        help='Minimum number of permutations to run with.')
    parser.add_argument('--max_permutations', type=int, default=10000,
                        help='Maximum number of permutations to run with.')
    parser.add_argument('--gene_set_size_min', type=int, default=10,
                        help='Minimum number of genes in geneset expressed in tissue.')
    parser.add_argument('--gene_set_size_max', type=int, default=1000,
                        help='Maximum number of genes in geneset expressed in tissue.')
    parser.add_argument('--significant_genes', type=str, required=True,
                        help='Path to file with new-line separated variants for your significant set.')
    parser.add_argument('--null_genes', type=str, required=True,
                        help='Path to file with new-line separated variants for your null set.')
    parser.add_argument('--gene_set_file', type=str, required=True,
                        help='Path to file with genesets.')
    parser.add_argument('--prefix', type=str, default='GeneEnrich_output',
                        help='Prefix for file output name.')
    parser.add_argument('--output_dir', type=str, default='',
                        help='Directory to write to. Defaults to working directory.')
    parser.add_argument('--fast_permute', type=str2bool, default=True,
                        help='If False, each geneset gets its own unique samples per permutation. Faster if True.')
    parser.add_argument('--HLA_remove', type=str2bool, default=False,
                        help='If True, removes genes that intersect with HLA region from Significant & Null set.')
    parser.add_argument('--HLA_file', type=str, default='../data/HLA_region_file_hg38.tsv',
                        help='Path to file indicating HLA region. Default file for HG38 provided.')
    parser.add_argument('--bed_remove', type=str,
                        help='Path to simple bed file with 3 unlabeled columns: chr, start, end. Significant '
                             '& Null genes that intersect with any of these intervals will be removed.')
    parser.add_argument('--gtf_file', type=str, required=True,
                        help='Path to gtf file to use for gene start & end definitions. Also used to find '
                             'gene names. Be certain gtf adheres to desired HG build version.')
    parser.add_argument('--p_val_cutoff', type=float, default=0.05,
                        help='p-value cutoff at which genesets are considered significant.')
    parser.add_argument('--q_val_cutoff', type=float, default=0.1,
                        help='q-value cutoff at which genesets are considered significant.')
    parser.add_argument('--fdr_cutoff', type=float, default=0.1,
                        help='fdr cutoff at which genesets are considered significant.')
    parser.add_argument('--genes_to_mark_table', type=str,
                        help=('Path to table with genes to mark. At the minimum, one column '
                              'should be labeled "gene". This column should have truncated '
                              'ensembl gene IDs. All columns will be annotated in output table.'))
    parser.add_argument('--mark_columns', type=str,
                        help=('Comma separated column names from genes_to_mark_table file to '
                              'include in main output table.'))
    parser.add_argument('--gene_expression_levels', type=str,
                        help=('Path to a file that indicates gene expression levels in the tissue '
                              'of interest. If this argument is passed, GeneEnrich will use '
                              'expression level as a confounding factor in the analysis.'))
    parser.add_argument('--tissue_of_interest', type=str, required=('-gene_expression_levels' in sys.argv),
                        help=('Tissue of interest to parse. GeneEnrich will ignore other tissue columns.'))
    parser.add_argument('--development', type=str2bool, default=False,
                        help=("Whether this is a development testing run. Determines whether to print hypergeometric null samples table."))
    parser.add_argument("--restrict_genes_database", action="store_true", help="Option to restrict gene-set enrichment analysis considering only genes in the input significant and null list present in the gene set database")
    parser.add_argument("--resource_name", type=str, default="",
                        help="Name of database resource, e.g. KEGG")
    parser.add_argument("--null_set", action="store_true",
                        help="Does not include sampling of significant genes.")
    parser.add_argument("--sig_gs_cutoff", action="store_true",
                        help="If false, use nominal gene set enrichment significance to select gene sets to be plotted in heatmap and genes/gene sets to be printed in gene centric table (as defined by --p_val_cutoff). If true, significant gene sets to be determined by BH FDR cutoff (as defined by --q_val_cutoff).")

    return parser.parse_args()

def count_in_rowwise(matr,vec,assume_sorted=False):
    """
    permutation algorithm. 
    matr: matrix of null values by permutation
    vec: geneset to check rows against
    """
    if assume_sorted==1:
        sorted_vec = vec
    else:
        sorted_vec = np.sort(vec)
    idx = np.searchsorted(sorted_vec,matr)
    idx[idx==len(sorted_vec)] = 0
    return (sorted_vec[idx] == matr).sum(1)

def permutation_manager(keys, sig_egenes, genesets, null_egenes, sample_size, min_permutations, 
                        max_permutations, geneset_sig_counts, fast_permute, gene_levels):
    """
    keys: list of geneset names
    sig_egenes: list of significant genes, for expr_level bins
    genesets: a groupby object of pandas dfs of genesets
    null_egenes: a list of all expressed genes (which genes will be sampled from). **NOT** only null genes
    min_permutations: number of permutations for this instance
    max_permutations: maximum number of permutations before breaking
    geneset_sig_counts: df containing intersection with significant set per geneset
    fast_permute: whether fast_permute is enabled
    gene_levels: a pd.Series with expression levels for all genes
    """
    print("Final number of working genesets: "+str(len(keys)))
    if len(keys) == 0:
        print("No workable genesets. Test cannot be run.")
        exit()
    #Check if we need to control on gene_levels
    if len(gene_levels) > 0:
        #put together an expression dataframe; index = ensembl_gene_id
        expr_df = pd.DataFrame({'expression_levels':gene_levels, 'bin':get_bins(gene_levels)},
                               index=gene_levels.index)
        expr_df['bin'] = expr_df['bin'].astype(str)
        #Figure out significant gene bins
        sig_df = expr_df.loc[expr_df.index.isin(sig_egenes)].copy()
        sig_bin_counts = sig_df[['bin']].groupby('bin').size().to_dict()
        #Remove significants from expr_df. no longer necessary
        expr_df = expr_df.loc[~(expr_df.index.isin(sig_egenes))].copy()
        #Pass a groupby object
        expr_df = expr_df.groupby('bin')
    else:
        print("No expression data passed. Sampling will not use gene expression data as confounding factor.")
        expr_df = pd.DataFrame([])
        sig_bin_counts = {}
    if not fast_permute: 
        row_list = []
        print("Sampling and intersecting with genesets...")
        # setup toolbar
        toolbar_width = 40
        tabs = np.ceil(len(keys)/toolbar_width)
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

        for idx, key in enumerate(keys):
            if idx % tabs == 0:
                sys.stdout.write("-")
                sys.stdout.flush()
            row = []
            #Get df for geneset
            group = genesets.get_group(key)
            #Sample from the nulls 
            egenes = list(set(group['ensembl_gene_id']))
            #Retrieve observed count
            obs_count = list(geneset_sig_counts.loc[geneset_sig_counts['gene_set']==key, 'number_significant_genes'])[0]
            row = multiple_random_sampling_for_egene(egenes, null_egenes, sample_size,
                                                      min_permutations, max_permutations,
                                                      obs_count, expr_df, sig_bin_counts)
            row_list.append(row)
        sys.stdout.write("\n")
        #change row indices
        df = pd.DataFrame(row_list)
        df.index=keys
        return(df)
    else:
        row_list = []
        print("Sampling and intersecting with genesets...")
        #Retrieve list of list of null samples, one list for each permutation
        initial_sample = single_random_sampling_for_egenes(null_egenes, sample_size, min_permutations, expr_df, sig_bin_counts)
        #Setup null_samples for intersection with genesets
        u, ids = np.unique(initial_sample, return_inverse=True)
        #Store in dictionary for retrieval
        null_samples_dict = {}
        null_samples_us = {}
        null_samples_ids = {}
        null_samples_dict[str(min_permutations)]=initial_sample
        null_samples_us[str(min_permutations)]=u
        null_samples_ids[str(min_permutations)]=ids
        #Setup toolbar
        toolbar_width = 40
        tabs = np.ceil(len(keys)/toolbar_width)
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
        for idx, key in enumerate(keys):
            if idx % tabs == 0:
                sys.stdout.write("-")
                sys.stdout.flush()
            #get df for geneset
            group = genesets.get_group(key)
            egenes = list(set(group['ensembl_gene_id']))
            #Retrieve observed count
            obs_count = list(geneset_sig_counts.loc[geneset_sig_counts['gene_set']==key, 'number_significant_genes'])[0]
            #Get the null data that we need
            row, null_samples_dict, null_samples_us, null_samples_ids  = iterative_overlap_with_geneset(null_samples_dict,
                                                                             null_samples_us, null_samples_ids, null_egenes,
                                                                             egenes, sample_size, min_permutations, max_permutations, 
                                                                             obs_count, expr_df, sig_bin_counts)
            row_list.append(row)
        sys.stdout.write("\n")
        #change row indices
        df = pd.DataFrame(row_list)
        df.index=keys
        return(df)

def iterative_overlap_with_geneset(null_samps_dict, null_samps_us, null_samps_ids, null_gene_list,
                                   gene_list, sample_size, num_samples, max_samples, obs_count, expr_df, sig_bin_counts):
    """
    dictionaries: stuff to retrieve initial null gene list
    gene_list: genes in the geneset
    sample_size: size of each sample
    num_samples: number of samples to take
    max_samples: max number of samples to take
    obs_count: observed intersection between geneset and significant genes
    expr_df: dataframe with expression levels. Can be empty if not provided
    sig_bin_counts: Dict containing number of draws to make from each bin
    outputs: 
    list of number of intersections between null samples and genesets
    updated inputs
    """
    null_samples_dict = null_samps_dict
    null_samples_us = null_samps_us
    null_samples_ids = null_samps_ids
    #If sample number already exists, don't sample again
    if str(num_samples) in null_samples_dict.keys():
        null_samples = null_samples_dict[str(num_samples)]
        u = null_samples_us[str(num_samples)]
        ids = null_samples_ids[str(num_samples)]
    else:
        fewer_samples = num_samples//10
        new_samples = num_samples - fewer_samples
        #Retrieve list of list of null samples, one list for each permutation
        null_samples1 = single_random_sampling_for_egenes(null_gene_list, sample_size, new_samples, expr_df, sig_bin_counts)
        #Setup null_samples for intersection with genesets
        u1, ids1 = np.unique(null_samples1, return_inverse=True)
        null_samples2 = null_samples_dict[str(fewer_samples)]
        u2 = null_samples_us[str(fewer_samples)]
        ids2 = null_samples_ids[str(fewer_samples)]
        null_samples = np.append(null_samples2, null_samples1, axis=0)
        u = np.append(u2, u1, axis=0)
        ids = np.append(ids2, ids1, axis=0)
        #Add values to dicts
        null_samples_dict[str(num_samples)] = null_samples
        null_samples_us[str(num_samples)] = u
        null_samples_ids[str(num_samples)] = ids
    #Get number of intersections between null and geneset
    samp_tf = count_in_rowwise(ids.reshape(null_samples.shape), ids[np.searchsorted(u, gene_list)])
    #Decide if more permutations are needed
    num_greater = sum([x > obs_count for x in samp_tf])
    if num_greater >= 5:
        return(samp_tf, null_samples_dict, null_samples_us, null_samples_ids)
    else:
        if num_samples >= max_samples:
            return(samp_tf, null_samples_dict, null_samples_us, null_samples_ids)
        else:
            samples = num_samples*10
            return(iterative_overlap_with_geneset(null_samples_dict, null_samples_us, null_samples_ids, null_gene_list,
                                                  gene_list, sample_size, samples, max_samples, obs_count, expr_df, sig_bin_counts))


def single_random_sampling_for_egenes(null_gene_list, sample_size, permutations, expr_df, sig_bin_counts):
    """
    null_gene_list: genes to sample from
    sample size: size of each sample
    permutations: number of permutations being done
    expr_df: dataframe with expression levels. Can be empty if not provided
    sig_bin_counts: dict containing number of draws to make from each bin
    """
    if len(expr_df) == 0:
        #Take random sample from uniform distribution, get indices by sorted numbers
        idx = np.random.rand(permutations, len(null_gene_list)).argsort(axis=1)[:,:sample_size]
        #Turn indices to genes
        samp_genes = np.take(null_gene_list, idx)
        return(samp_genes)
    else: #If controlling on expression levels
        return_genes = None
        for key in sig_bin_counts.keys():
            null_gene_list = list(expr_df.get_group(key).index)
            #Take random sample from uniform distribution, get indices by sorted numbers
            idx = np.random.rand(permutations, len(null_gene_list)).argsort(axis=1)[:,:sig_bin_counts[key]]
            #Turn indices to genes
            samp_genes = np.take(null_gene_list, idx)
            try:
                #Append each list by permutation number
                return_genes = np.append(return_genes, samp_genes, axis=1)
            except:
                return_genes = samp_genes 
        return(return_genes)
      
def multiple_random_sampling_for_egene(gene_list, null_gene_list, sample_size, 
                                       num_samples, max_samples, obs_count, expr_df, sig_bin_counts):
    """
    gene_list: genes in the geneset
    null_gene_list: genes to sample from
    sample_size: size of each sample
    num_samples: number of samples to take
    max_samples: max number of samples to take
    obs_count: observed intersection between geneset and significant genes
    expr_df: dataframe with expression levels. Can be empty if not provided
    sig_bin_counts: dict containing number of draws to make from each bin
    outputs: 
    list of number of intersections between null samples and genesets
    """
    if len(expr_df) == 0:
        #Take random sample from uniform distribution, get indices by sorted numbers
        idx = np.random.rand(num_samples, len(null_gene_list)).argsort(axis=1)[:,:sample_size]
        #Turn indices to genes
        samp_genes = np.take(null_gene_list, idx)
    else:
        samp_genes = None
        for key in sig_bin_counts.keys():
            null_gene_list = list(expr_df.get_group(key).index)
            #Take random sample from uniform distribution, get indices by sorted numbers
            idx = np.random.rand(num_samples, len(null_gene_list)).argsort(axis=1)[:,:sig_bin_counts[key]]
            #Turn indices to genes
            sub_samp_genes = np.take(null_gene_list, idx)
            try: 
                #Append each list by permutation number
                samp_genes = np.append(samp_genes, sub_samp_genes, axis=1)
            except:
                samp_genes = sub_samp_genes
    #Prepare for intersections
    u, ids = np.unique(samp_genes, return_inverse=True)
    #Get number of intersections between null and geneset
    samp_tf = count_in_rowwise(ids.reshape(samp_genes.shape), ids[np.searchsorted(u, gene_list)])
    #Decide if more permutations are needed
    num_greater = sum([x > obs_count for x in samp_tf])
    if num_greater >= 5:
        return(samp_tf)
    else:
        if num_samples >= max_samples:
            return(samp_tf)
        else:
            samples = num_samples*10
            return(multiple_random_sampling_for_egene(gene_list, null_gene_list, sample_size, samples, 
                                                      max_samples, obs_count, expr_df, sig_bin_counts))
 
def hypergeometric(x, M, df, col_num):
    """
    x: Hypergeometric population numerator (num sig genes)
    M: hypergeometric population denominator (total number of expressed genes)
    df: dataframe to perform hypergeometric test on. Should have a
    column 'size', which will be used for sample denominator. All
    other columns will be used as sample numerators.
    col_num: number of columns to consider. 
    """
    size_series = pd.Series(df['size'], index=df.index)
    df_use = df.drop('size', axis=1).copy()
    df_use = df_use.iloc[:,:col_num].copy()
    partial_row = partial(hypergeometric_by_row, x=x, M=M, size_series = size_series)
    df_hyper_p = df_use.apply(partial_row, axis=1)
    return(df_hyper_p)

def hypergeometric_by_row(row, x, M, size_series):
    """
    inner function for hypergeometric
    """
    N = size_series.loc[size_series.index == row.name].values[0] #geneset size
    p_values = row.apply(lambda v: 1-(scipy.stats.hypergeom.cdf((v-1),M,x,N)))
    return(p_values)

def find_empirical_p_value(row, sig_val_df):
    """
    row: row from permutation df with normalized values
    sig_val_df: df containing significant normalized values for genesets, or float
    """
    sig_val = sig_val_df.loc[sig_val_df.index == (row.name)].values[0]
    list_row = np.array(list(row))
    list_row = np.append(list_row, sig_val) #We added the significant value to the row
    no_nan_row = list_row[~np.isnan(list_row)] #Remove any nan values
    len_no_nan = float(len(no_nan_row) - 1) #Number of non-nan values has to exclude sig value
    no_nan_row = -np.sort(-no_nan_row) #Sort in descending order
    match_idx = np.where(no_nan_row==sig_val) #where values = sig
    num_great_eq = float(max(match_idx[0])) #How many are >= sig
    empirical_p_val = num_great_eq/len_no_nan
    return(empirical_p_val)

#find_fdr and compare_within_fdr function separately from above functions
def find_fdr(df, empirical):
    """
    df: df from permutation with hypergeometric p-values
    empirical: Series containing hypergeometric p-values for genesets
    """
    fdr = empirical.copy()
    all_df_vals = df.values.flatten()
    num_denum = float(len(all_df_vals)) #numerator's denominator
    denum_denum = float(len(empirical.index)) #denominator's denominator
    sig_df_vals = np.array(empirical)
    sig_array_sorted = np.sort(sig_df_vals)
    #Make the whole thing one long list, sort in ascending order
    no_nan_vals = np.append(all_df_vals, sig_df_vals)
    no_nan_vals = np.sort(no_nan_vals)
    #compute fdr
    for idx in empirical.index:
        p_val = empirical.loc[empirical.index==idx].values[0]
        geneset_den_num = (float(max(np.where(sig_array_sorted==p_val)[0]))+1.0)
        geneset_num_num = ((float(max(np.where(no_nan_vals==p_val)[0])))-geneset_den_num)
        fdr.loc[fdr.index==idx] = ((geneset_num_num/num_denum)/(geneset_den_num/denum_denum))
    return(fdr)


def annotate_sig_genes(geneset_df, p_sig, all_genesets, all_genes_df, sig_genes, output_df, gene_dict, 
                       gene_colname='name', add_gene_col=True,
                       sig_colname='name', add_sig_col=True, mark_df=None,
                       mark_columns=None):
    """
    Adds a column indicating what genes are significant to the output df
    Returns the output df with the column added, and a separate table
    with unique gene per row
    Recording whether a geneset passes a certain significance cutoff (T/F)
    geneset_df: mostly the inputted geneset database resource table, with 
    column to indicate gene significance.
    p_sig: genesets that are significant in this case
    all_genesets: all genesets
    all_genes_df: input file with all genes
    sig_genes_df: input file with significant genes
    output_df: output df that's ready besides this function
    gene_dict: dictionary with gene objects inside
    gene_colname: What you want the new gene column to be named
    sig_colname: What you want the new significant T/F column to be named
    add_gene_col: whether you want to add column to include all significant
    genes in the geneset. 
    add_sig_col: whether you want to add a significant T/F column
    mark_df: data frame with genes you want annotated. Should have, at least,
    column called 'gene'. Overrides significant gene column writing.
    mark_columns: columns you want marked in the main output table
    """
    sig_genes_df = all_genes_df.drop_duplicates()
    #Get significant genes from input gene files from significant genesets
    #p_sig_df = geneset_df.loc[geneset_df['gene_set'].isin(set(all_genesets))].copy()#&
                           #   (geneset_df['significant']==True))].copy()
    gene_name_dict = {key:value.gene_name for key,value in gene_dict.items()}
    #Annotate with any other information from input
    p_sig_df = []
    if add_sig_col:
        output_df[sig_colname] = False
    #Annotate output_df
    for geneset in set(all_genesets):
        p_sig_df_gs = geneset_df.loc[(geneset_df['gene_set']==geneset) & (geneset_df['significant']==True)].copy()
        p_sig_df_gs.drop('significant', axis=1, inplace=True)
        sig_info = sig_genes_df.loc[sig_genes_df['gene'].isin(set(p_sig_df_gs['ensembl_gene_id']))].copy()
        sig_info = sig_info.loc[sig_info[['gene']].drop_duplicates().index].copy()
        if mark_df is not None:
            sig_info = mark_df.loc[mark_df['gene'].isin(set(sig_info['gene']))].copy().merge(sig_info, how='left', on='gene').copy()
            sig_info.dropna(inplace=True)
            #Get columns to mark
            if sig_info is None:
                continue #If there is nothing to annotate, skip
            if mark_columns is not None:
                sig_info = sig_info[mark_columns].drop_duplicates().copy()
                p_sig_df_gs = p_sig_df_gs.merge(sig_info, how='right',
                                                left_on='ensembl_gene_id', right_on='gene').copy()
                p_sig_df_gs = p_sig_df_gs.merge(sig_genes_df, how='left',
                                                left_on='ensembl_gene_id', right_on='gene').copy()
        if add_gene_col:
            #turn gene ids to gene names, then add column
            sig_info['gene'] = sig_info['gene'].map(gene_name_dict)
            if len(sig_info.gene) > 0:
                output_df.loc[output_df.index == geneset, gene_colname] = df_to_delim(sig_info)
            else:
                output_df.loc[output_df.index == geneset, gene_colname] = 'None'
            #count number of marked genes
            if mark_columns is not None:
                output_df.loc[output_df.index == geneset, 'marked_significant_gene_count'] = str(len(set(sig_info.gene))) 
        if add_sig_col:
            if geneset in p_sig:
                output_df.loc[output_df.index == geneset, sig_colname] = True
        p_sig_df.append(p_sig_df_gs)
    if len(p_sig_df) > 1:
        p_sig_df = pd.concat(p_sig_df)
    elif len(p_sig_df) == 1:
        p_sig_df = p_sig_df[0]
    if len(p_sig_df) > 0:
        p_sig_df['significant_gene'] = False
        if 'ensembl_gene_id' in p_sig_df.columns: 
            p_sig_df['gene_name'] = p_sig_df['ensembl_gene_id'].map(gene_name_dict) 
            p_sig_df.loc[p_sig_df['ensembl_gene_id'].isin(list(set(sig_genes.gene))), 'significant_gene'] = True
        elif 'gene' in p_sig_df.columns: 
            p_sig_df['gene_name'] = p_sig_df['gene'].map(gene_name_dict) 
            p_sig_df.loc[p_sig_df['gene'].isin(list(set(sig_genes.gene))), 'significant_gene'] = True
        p_sig_df.dropna(inplace=True)
    return(output_df, p_sig_df)

def retrieve_resource(geneset_file, resource_name):
    """
    retrieves resource name from file
    first index and sets to uppercase

    ARH 29 January 2021
    """
    if not resource_name:
        resource = os.path.basename(geneset_file)
        resource = resource.split("_")[0]
        return(resource.upper())
    return(resource_name.upper())

def retrieve_key_from_value(dict, value):
    """
    returns key from a value in dictionary
    ARH 29 January 2021
    """ 
    for key, val in dict.items():
        if val == value:
            return(key)
    return('')

def retrieve_variant_from_gene(dict, sig_genes_df, gene_name):
    """
    writes df to dictionary {gene: variant}

    returns variant value from gene key

    ARH 29 Jan 2021
    """
    #TODO:
    # return NA if variant not found
    gene_variant_dict = pd.Series(sig_genes_df.variant.values,index=sig_genes_df.gene).to_dict()
    gene_symbol = retrieve_key_from_value(dict, gene_name)
    if gene_symbol:
        try:
            variant = gene_variant_dict[gene_symbol]
            return(variant)
        except KeyError:
            return('NA')
    return('NA')

def extract_genesets_lists(df, geneset_list, col_name):
    """
    selects genesets that pass significance cutoff
    if not, sets to NA
    if only NA in string, returns one NA
    """
    #geneset_string = ';'.join(set([x if list(df.loc[df['gene_set'] == x, col_name])[0] == True else 'NA' for x in geneset_list.split(';')]))
    geneset_string = ';'.join(set([x for x in geneset_list.split(';') if list(df.loc[df['gene_set'] == x, col_name])[0] == True]))

    # [x if list(df.loc[df['gene_set'] == x, col_name])[0] == True else 'NA' for x in geneset_list.split(';')]
    #make sure unique genesets
    #double check 5 Mar 2021
    # if empty list: return NA
    geneset_string_split = geneset_string.split(';')
    if not geneset_string_split:
        return('NA')
    if not geneset_string:
        return('NA')
    return(geneset_string)

def format_directory(dir):
    """
    adds '/' if not at end of dir

    ARH 5 Mar 2021
    """
    if not dir.endswith('/'):
        if not dir == '':
            return(dir + '/')
        return(dir)
    return(dir)

def set_mark_columns(mark_columns):
    """
    split mark columns 

    ARH 5 Mar 2021
    """
    #Get all annotations to mark
    if mark_columns is not None:
        if ',' in mark_columns:
            return(mark_columns.split(','))
        return(mark_columns)
    return(mark_columns)

def parse_gencode_gtf(gtf_file):
    """
    manages parse_gtf function

    ARH 5 Mar 2021
    """
    if gtf_file.endswith('gz'):
        gene_dict = parse_gtf(gtf_file, 'gz')
        return(gene_dict)
    elif gtf_file.endswith('gtf'):
        gene_dict = parse_gtf(gtf_file)
        return(gene_dict)
    else:
        raise Exception('GTF file passed cannot be parsed. Please pass a .gtf or .gz.')

def mark_genes(mark_genes_file):
    """
    processes mark_genes_file if inputted
    """
    if mark_genes_file:
        mark_genes_df = pd.read_csv(mark_genes_file, sep='\t')
        if '.' in list(mark_genes_df['gene'])[0]:
            mark_genes_df['gene'] = mark_genes_df['gene'].apply(lambda x: x.split('.')[0])
            return(mark_genes_df)
        return(mark_genes_df)
    return([])

def verify_sig_null_files(sig_genes_df, null_genes_df):
    """
    ensures column headers are correct for sig_genes_df and null_genes_df
    """
    if not null_genes_df.columns.values[0]=='gene':
        raise Exception("Input file for null genes appears to be formatted incorrectly. Please check input.")
    if not ('gene' in sig_genes_df.columns.values):
        raise Exception("Input file for significant genes appears to be formatted incorrectly. Please check input.")

def verify_intron_cluster(sig_genes_df):
    """
    if intron_cluster header in sig_genes_df, return dictionary
    """
    if 'intron_cluster' in sig_genes_df.columns:
        ensembl_intron_dict = pd.Series(sig_genes_df.intron_cluster.values,index=sig_genes_df.gene).to_dict()
        return(ensembl_intron_dict)
    return(False)

def truncate_ensembl_gene_ids(df, col_name):
    """
    if version number in ensembl gene id, remove

    e.g. ENSG00013432.8 -> ENSG00013432
    """
    if '.' in list(df[col_name])[0]:
        df[col_name] = df[col_name].apply(lambda x: x.split('.')[0])
        return(df)
    return(df)       

def extract_genes_expressed(df):
    """
    return list of genes and their lengths
    """
    list_genes = list(set(df['gene']))
    return(list_genes, len(list_genes))

def remove_hla_genes(HLA_remove, HLA_file, gene_values, sig_egenes, null_egenes):
    """
    if user inputs HLA_remove -> return list of genes
    """
    if HLA_remove:
        print("Removing HLA genes from significant & null set...")
        HLA_list = parse_interval_list(HLA_file)
        for interval in HLA_list:
            chr = interval.chr
            range_interval = set(interval.interval_range)
            HLA_genes = [gene for gene in gene_values if gene.chr==chr]
            HLA_genes = [gene for gene in HLA_genes if range_interval.intersection(gene.gene_range)]
            HLA_gene_names = [gene.gene_id for gene in HLA_genes]
        #find gene intersections with significant and null groups
        sig_HLA_genes = [gene for gene in HLA_gene_names if gene in sig_egenes]
        null_HLA_genes = [gene for gene in HLA_gene_names if gene in null_egenes]
        print("HLA genes removed from significant set: "+str(len(sig_HLA_genes)))
        print("HLA genes removed from null set: "+str(len(null_HLA_genes)))
        return(sig_HLA_genes, null_HLA_genes)
    return([], [])

def remove_bed_file_genes(bed_file, sig_egenes, null_egenes):
    """
    removes genes from bed file intersection
    """
    if bed_file:
        print("Removing significant & null genes that intersect intervals from bed_file...")
        interval_list = parse_interval_list(bed_file)
        #flatten intervals by chr, then compare to gencode
        interval_chrs = list(set([interval.chr for interval in interval_list]))
        bed_genes = []
        for chr in interval_chrs:
            chr_intervals = [interval for interval in interval_list if interval.chr == chr]
            chr_region = [interval.interval_range for interval in chr_intervals]
            chr_region = [item for sublist in chr_region for item in sublist]
            chr_region = set(chr_region) #all values in range. Intersect with genes
            chr_genes = [gene for gene in gene_values if gene.chr==chr]
            chr_genes = [gene for gene in chr_genes if chr_region.intersection(gene.gene_range)]
            bed_genes = bed_genes + [gene.gene_id for gene in chr_genes]
        sig_bed_genes = [gene for gene in bed_genes if gene in sig_egenes]
        null_bed_genes = [gene for gene in bed_genes if gene in null_egenes]
        print("bed file genes removed from significant set: "+str(len(sig_bed_genes)))
        print("bed file genes removed from null set: "+str(len(null_bed_genes)))
        return(sig_bed_genes, null_bed_genes)
    return([], [])

def identify_number_unique_genes_significant_set(df, col_name):
    """
    flattens list of genes that are true for a cutoff
    returns unique length of flattened list
    """
    try:
        subset_df = df.loc[df[col_name] == True][['significant_genes', col_name]].copy()
        subset_df['gene_list'] = subset_df['significant_genes'].apply(lambda x: x.split(';'))
        list_sublists = list(subset_df['gene_list'].values)
        flat_list = [item for sublist in list_sublists for item in sublist]
        return(list(set(flat_list)))
    except:
        return('0')

def identify_length_number_unique_genes_significant_set(df, col_name):
    """
    flattens list of genes that are true for a cutoff
    returns unique length of flattened list
    """
    flat_list = identify_number_unique_genes_significant_set(df, col_name)
    if type(flat_list) is list:
        return(str(len(flat_list)))
    return('0')

def populate_rows_gene_centric_table(adv_df, my_genes, p_sig_df, gene_name_dict, sig_genes_df, gsc_df):
    """
    each row corresponds to one gene
    """
    for gene in my_genes: #Make 1 line per gene
        #Get list of gene-sets the gene belongs to
        #ensure genesets are unique
        gss = ';'.join(list(set(adv_df.loc[adv_df.significant_genes.apply(lambda x: gene in x.split(';'))]['gene_set'])))

        if len(gss) > 0:
            # INITIALIZE AS NA
            #Get list of nominal significant gene-sets the gene belongs to
            ngss = extract_genesets_lists(adv_df, gss, 'pass_nominal_significance')

            #Get list of BH significant gene-sets the gene belongs to
            bhgss = extract_genesets_lists(adv_df, gss, 'pass_bh_significance')

            #Get list of FDR significant gene-sets the gene belongs to
            fdrgss = extract_genesets_lists(adv_df, gss, 'pass_fdr_significance')

            #Set most significant geneset gene belongs to
            #TODO: investigate
            bgss = list(set(p_sig_df.loc[p_sig_df.gene_name == gene, 'gene_set']))

            idx = []
            for i in bgss:
                idx.append(adv_df.loc[adv_df.gene_set == i].index)
            try:
                min_idx = min(idx)
                bgss = list(adv_df.iloc[min_idx]['gene_set'])[0]
            except:
                bgss = 'NA'
            #Get the most significant gene-set p_value
            try:
                pgs = "{:.2e}".format(float(list(adv_df.loc[adv_df['gene_set'] == bgss, 'empirical_pval'])[0])) # ARH 29 Jan 2021
            except:
                pgs = 'NA'
            #For sorting, determine genes with most nominal significant genesets
            try:
                # if na, do not count
                num_ngss = len(ngss.split(';'))
            except:
                num_ngss = 0

        else: #If gene does not belong to any gene-set
            gss = 'NA'
            ngss = 'NA'
            bhgss = 'NA'
            fdrgss = 'NA'
            bgss = 'NA'
            pgs = 'NA'
            num_ngss = 0


        #make the row
        if 'variant' in p_sig_df:
            row_to_add = pd.Series({'gene':gene,
                                    'ensembl_gene_id': retrieve_key_from_value(gene_name_dict, gene),
                                    'variant': retrieve_variant_from_gene(gene_name_dict, sig_genes_df, gene),
                                    'number_gene_sets_pass_nominal_significance': num_ngss,
                                    'most_significant_gene_set_membership': bgss,
                                    'most_significant_gene_set_empirical_p_value': pgs, # 8 Feb 2021 ARH changed hypergemeotric to empiricail
                                    'gene_set_membership': gss,
                                    'pass_nominal_significance_gene_sets': ngss,
                                    'pass_bh_significance_gene_sets': bhgss,
                                    'pass_fdr_significance_gene_sets': fdrgss})
        else:
            row_to_add = pd.Series({'gene':gene,
                                    'ensembl_gene_id': retrieve_key_from_value(gene_name_dict, gene),
                                    'number_gene_sets_pass_nominal_significance': num_ngss,
                                    'most_significant_gene_set_membership': bgss,
                                    'most_significant_gene_set_empirical_p_value': pgs, # 8 Feb 2021 ARH changed hypergemeotric to empiricail
                                    'gene_set_membership': gss,
                                    'pass_nominal_significance_gene_sets': ngss,
                                    'pass_bh_significance_gene_sets': bhgss,
                                    'pass_fdr_significance_gene_sets': fdrgss})
        gsc_df = gsc_df.append(row_to_add, ignore_index=True).copy()
    return(gsc_df)

def retrieve_intron_cluster(intron_cluster, ensembl_gene_id):
    """
    retrieves intron cluster from dictionary
    """
    try:
        return(intron_cluster[ensembl_gene_id])
    except:
        return('NA')

def add_intron_cluster_header_gsc_df(gsc_df_populated, intron_cluster):
    """
    adds intron_cluster header name
    """
    gsc_df_populated['intron_cluster'] = gsc_df_populated['ensembl_gene_id'].apply(lambda x: retrieve_intron_cluster(intron_cluster, x))

    gsc_df_populated = gsc_df_populated[['gene', 'ensembl_gene_id', 'intron_cluster','variant', 'number_gene_sets_pass_nominal_significance',
                                         'most_significant_gene_set_membership', 'most_significant_gene_set_empirical_p_value', 
                                         'gene_set_membership', 'pass_nominal_significance_gene_sets', 'pass_bh_significance_gene_sets', 
                                         'pass_fdr_significance_gene_sets', 'resource']]
    return(gsc_df_populated)

def declare_empty_gsc_df(p_sig_df):
    """
    if variant in p_sig_df ensure variant column is present
    """
    if 'variant' in p_sig_df:
        gsc_df = pd.DataFrame({'gene':[],
                               'ensembl_gene_id':[],
                               'variant':[],
                               'number_gene_sets_pass_nominal_significance': [],
                               'most_significant_gene_set_membership': [],
                               'most_significant_gene_set_empirical_p_value': [], # ARH 8 Feb 2021 changed hypergeometric to empirical
                               'gene_set_membership': [],
                               'pass_nominal_significance_gene_sets': [],
                               'pass_bh_significance_gene_sets': [],
                               'pass_fdr_significance_gene_sets': []})
        return(gsc_df)
    gsc_df = pd.DataFrame({'gene':[],
                           'ensembl_gene_id':[],
                           'number_gene_sets_pass_nominal_significance': [],
                           'most_significant_gene_set_membership': [],
                           'most_significant_gene_set_empirical_p_value': [], # ARH 8 Feb 2021 changed hypergeometric to empirical
                           'gene_set_membership': [],
                           'pass_nominal_significance_gene_sets': [],
                           'pass_bh_significance_gene_sets': [],
                           'pass_fdr_significance_gene_sets': []})
    return(gsc_df)

def curate_gene_centric_table(output_df, my_genes, p_sig_df, gene_name_dict, sig_genes_df, resource, intron_cluster):
    """
    creates gene centric table
    """
    adv_df = output_df.copy()
    adv_df['gene_set'] = adv_df.index
    adv_df.reset_index(inplace=True, drop=True)
    #To check for gene membership

    adv_df[['significant_genes']] = adv_df[['significant_genes']].fillna(value='None').copy()
    adv_df['sig_gene_list'] = adv_df.significant_genes.apply(lambda x: x.split(';'))
    gsc_df = pd.DataFrame({'gene':[],
                           'ensembl_gene_id':[],
                           'variant':[],
                           'number_gene_sets_pass_nominal_significance': [],
                           'most_significant_gene_set_membership': [],
                           'most_significant_gene_set_empirical_p_value': [], # ARH 8 Feb 2021 changed hypergeometric to empirical
                           'gene_set_membership': [],
                           'pass_nominal_significance_gene_sets': [],
                           'pass_bh_significance_gene_sets': [],
                           'pass_fdr_significance_gene_sets': []})

    gsc_df_populated = populate_rows_gene_centric_table(adv_df, my_genes, p_sig_df, gene_name_dict, sig_genes_df, gsc_df)

    gsc_df_populated.sort_values(by='number_gene_sets_pass_nominal_significance', ascending=False, inplace=True)
    gsc_df_populated = gsc_df_populated.loc[gsc_df_populated["number_gene_sets_pass_nominal_significance"] > 0]
    gsc_df_populated['number_gene_sets_pass_nominal_significance'] = gsc_df_populated['number_gene_sets_pass_nominal_significance'].astype(int) # ARH 8 February column as integer
    gsc_df_populated['resource'] = resource # ARH 8 Feb 2021 added in resource column

    if type(intron_cluster) is dict:
        gsc_df_populated_added_intron = add_intron_cluster_header_gsc_df(gsc_df_populated, intron_cluster)
        return(gsc_df_populated_added_intron)
    return(gsc_df_populated)

def identify_p_sig_df_columns(variant = True, intron_cluster = False):
    """
    selects relevant columns depending on presence of variant and intron cluster
    """
    if variant and intron_cluster:
        column_list = ['resource', 'gene_set', 'total_number_genes',
                       'ensembl_gene_id', 'EntrezID', 'gene_name',
                       'significant_gene', 'pass_nominal_significance', 'pass_bh_significance', 
                       'pass_fdr_significance', 'intron_cluster', 'variant']
        return(column_list)
    elif variant and not intron_cluster:
        column_list = ['resource', 'gene_set', 'total_number_genes',
                       'ensembl_gene_id', 'EntrezID', 'gene_name', 
                       'significant_gene', 'pass_nominal_significance', 'pass_bh_significance', 
                       'pass_fdr_significance', 'variant']
        return(column_list)
    elif not variant and not intron_cluster:
        column_list = ['resource', 'gene_set', 'total_number_genes',
                       'ensembl_gene_id', 'EntrezID', 'gene_name',
                       'significant_gene', 'pass_nominal_significance',
                       'pass_bh_significance', 'pass_fdr_significance']
        return(column_list)
    else: #not variant and intron_cluster
        column_list = ['resource', 'gene_set', 'total_number_genes',
                       'ensembl_gene_id', 'EntrezID', 'gene_name',
                       'significant_gene', 'pass_nominal_significance', 'pass_bh_significance', 
                       'pass_fdr_significance', 'intron_cluster']
        return(column_list)


def main():
    """
    main section
    """

    # retrieve args
    args = parse_args()

    #Define date
    today = datetime.date.today().strftime('_%d_%b_%Y')

    #Retrieve user-defined parameters
    min_permutations = args.min_permutations
    max_permutations = args.max_permutations
    min_geneset = args.gene_set_size_min
    max_geneset = args.gene_set_size_max
    sig_genes_file = args.significant_genes
    null_genes_file = args.null_genes
    geneset_file = args.gene_set_file
    prefix = args.prefix 
    fast_permute = args.fast_permute
    p_val_cutoff = args.p_val_cutoff
    q_val_cutoff = args.q_val_cutoff
    fdr_cutoff = args.fdr_cutoff
    mark_genes_file = args.genes_to_mark_table
    mark_columns = args.mark_columns
    development = args.development
    output_dir = args.output_dir

    #Arguments for parsing intervals
    HLA_remove = args.HLA_remove
    HLA_file = args.HLA_file
    bed_file = args.bed_remove
    gtf_file = args.gtf_file

    resource = retrieve_resource(geneset_file, args.resource_name) #ARH 29 Jan 2021

    output_dir = format_directory(output_dir) # ARH 5 Mar 2021

    #arguments for gene_expression_levels
    gene_expression_levels = args.gene_expression_levels
    tissue_of_interest = args.tissue_of_interest

    #mark columns
    mark_columns = set_mark_columns(mark_columns) # ARH 5 Mar 2021

    #Parse GENCODE gtf if necessary
    print("Parsing GTF file...")
    gene_dict = parse_gencode_gtf(gtf_file) # ARH 5 Mar 2021
    gene_values = np.array(list(gene_dict.values()))

    #read in various files
    sig_genes_df = pd.read_csv(sig_genes_file, sep='\t', dtype={'gene':str})
    null_genes_df = pd.read_csv(null_genes_file, sep='\t', dtype={'gene':str})
    geneset_df = pd.read_csv(geneset_file, sep='\t')#.iloc[:125000,:]

    # prepare mark_genes_file
    mark_genes_df = mark_genes(mark_genes_file) # ARH 5 Mar 2021

    #verify sig and null files are correct
    verify_sig_null_files(sig_genes_df, null_genes_df) #ARH 5 Mar 2021

    #Verify no '.' in gene ids ARH 5 Mar 2021
    print("Truncating ensembl gene ID version numbers from input as necessary.")
    sig_genes_df = truncate_ensembl_gene_ids(sig_genes_df, 'gene')
    null_genes_df = truncate_ensembl_gene_ids(null_genes_df, 'gene')
    geneset_df = truncate_ensembl_gene_ids(geneset_df, 'ensembl_gene_id')
    full_list = list(geneset_df["ensembl_gene_id"].unique())

    #verify if intron_cluster in sig_genes_df
    intron_cluster = verify_intron_cluster(sig_genes_df)

    #To add extra columns in gene tables
    concat_later = pd.concat([sig_genes_df, null_genes_df], sort = True)
    concat_later = concat_later.reset_index(drop=True)
    concat_later = concat_later.drop_duplicates(["gene"], keep="first")

    #find which genes are expressed in our tissue
    sig_egenes, len_sig_egenes = extract_genes_expressed(sig_genes_df)
    null_egenes, len_null_egenes = extract_genes_expressed(null_genes_df)

    #Remove HLA genes if necessary
    sig_HLA_genes, null_HLA_genes = remove_hla_genes(HLA_remove, HLA_file, gene_values, sig_egenes, null_egenes)

    #Remove bed file intersection genes if necessary
    sig_bed_genes, null_bed_genes = remove_bed_file_genes(bed_file, sig_egenes, null_egenes)

    #Find expression levels if provided
    if gene_expression_levels is not None:
        expr_tissue = tissue_of_interest 
        gene_expr_df = pd.read_csv(gene_expression_levels, sep='\t')
        #get only tissue of interest
        if expr_tissue not in gene_expr_df.columns:
            raise Exception("Indicated tissue of interest is not in expression file. Check your inputs.")
            exit()
        else:
            gene_expr_df = gene_expr_df[['Name', expr_tissue]].copy()
            gene_expr_df.columns = ['gene', 'expression_level']
            #Remove decimal if necessary
            if '.' in list(gene_expr_df['gene'])[0]:
                gene_expr_df['gene'] = gene_expr_df['gene'].apply(lambda x: x.split('.')[0])
            #Find only genes allegedly expressed
            gene_levels = pd.Series(gene_expr_df.expression_level.values, index=gene_expr_df.gene.values) 
            gene_levels.dropna(inplace=True)
            #Check if any genes do not have expression information
            sig_expr_genes = [gene for gene in sig_egenes if gene not in gene_levels.index]
            null_expr_genes = [gene for gene in null_egenes if gene not in gene_levels.index]            
            #Let user know if any were dropped
            if (len(sig_expr_genes)>0) or (len(null_expr_genes)>0):
                print("NOTE: "+str(len(sig_expr_genes))+" genes from significant set and "+
                      str(len(null_expr_genes))+" genes from null set do not have expression "+
                      "values in provided file for tissue of interest. They are not included in analysis.")
    else: 
        sig_expr_genes = []
        null_expr_genes = []
        gene_levels = pd.Series([])

    #Final removal if genes as necessary 
    remove_sig_genes = set(sig_HLA_genes+sig_bed_genes+sig_expr_genes)
    remove_null_genes = set(null_HLA_genes+null_bed_genes+null_expr_genes)
    sig_egenes = set(set(sig_egenes) - set(remove_sig_genes))
    null_egenes = set(set(null_egenes) - set(remove_null_genes))
    if ((len(remove_sig_genes) > 0) or (len(remove_null_genes)>0)):
        print("\nFinal significant gene count removed: "+str(len(remove_sig_genes)))
        print("Final null gene count removed: "+str(len(remove_null_genes)))

    intersect_genes = set(sig_egenes & null_egenes)
    if len(intersect_genes) > 0:
        print("NOTE: "+str(len(intersect_genes))+" genes occur both in your submitted "+
              "significant and null sets. These will be removed from the null set.")
        null_egenes = null_egenes - intersect_genes

    #Check remaining genes
    sig_egenes = list(sig_egenes)
    null_egenes = list(null_egenes)
    genes_expressed = list(set(sig_egenes + null_egenes))

    #Remove genes that aren't part of sig/null from expression series
    gene_levels = gene_levels.loc[gene_levels.index.isin(genes_expressed)].copy()
    geneset_df = geneset_df.loc[geneset_df['ensembl_gene_id'].isin(genes_expressed)].copy()
    

    only_genes_df = geneset_df[['ensembl_gene_id']].copy()
    only_genes_df.columns = ['gene']

    print("Selecting genesets with enough genes expressed in the tissue...")
    #find which genesets to keep
    genesets = geneset_df.groupby('gene_set')
    geneset_df = genesets.filter(lambda x: len(x) >= min_geneset).copy()
    genesets = geneset_df.groupby('gene_set')
    geneset_df = genesets.filter(lambda x: len(x) <= max_geneset).copy()
    geneset_df['size']=geneset_df.groupby('gene_set')['gene_set'].transform('count').copy()
    genesets = geneset_df.groupby('gene_set')

    #Find degree of intersection per geneset, turn into DataFrame
    print("Finding degree of intersection between significant set and genesets...")
    geneset_df['significant'] = geneset_df['ensembl_gene_id'].isin(sig_egenes)
    temp = geneset_df.groupby('gene_set')[['significant']].sum()
    geneset_sig_counts = pd.DataFrame({'gene_set':temp.index,
                                        'number_significant_genes':temp.values[:,0].astype(int)})

    #permute for each geneset
    print("Preparing permutation...")
    keys = genesets.groups.keys() #geneset dfs

    if args.restrict_genes_database:
        sig_overlap = [gene for gene in sig_egenes if gene in full_list]

        #if args.null_set:
        #    genes_expressed_overlap = [gene for gene in null_egenes if gene in full_list]
        #else:
        genes_expressed_overlap = [gene for gene in genes_expressed if gene in full_list]
        sample_size = len(sig_overlap)
        x = len(sig_overlap) #hypergeometric population numerator
        n = len(genes_expressed_overlap) #hypergeometric population denominator
        sig_genes_input = sig_overlap
        genes_expressed_input = genes_expressed_overlap
    else:
        sample_size = len(sig_egenes)
        x = len(sig_egenes) #hypergeometric population numerator
        #if args.null_set:
        #    n = len(null_egenes)
        #    genes_expressed_input = null_egenes
        #else:
        n = len(genes_expressed) #hypergeometric population denominator
        genes_expressed_input = genes_expressed
        sig_genes_input = sig_egenes

    k_df = permutation_manager(keys, sig_genes_input, genesets, genes_expressed_input, sample_size,
                               min_permutations, max_permutations, geneset_sig_counts,
                               fast_permute, gene_levels)

    #Prepare significant df for hypergeometric 
    sig_set_df = geneset_df[['gene_set', 'size']].drop_duplicates().merge(geneset_sig_counts).copy()
    sig_set_df.index = sig_set_df['gene_set']
    sig_set_df.drop('gene_set', axis=1, inplace=True)

    #Prepare null df for hypergeometric
    null_set_df = k_df.merge(geneset_df[['gene_set', 'size']].drop_duplicates(), left_index=True, 
                             right_on='gene_set').copy()
    null_set_df.index = null_set_df.gene_set
    null_set_df.drop('gene_set', axis=1, inplace=True)

    print("Computing hypergeometric p-values for observed values...")
    hyper_p_sig = hypergeometric(x, n, sig_set_df, min_permutations)
    hyper_p_series = pd.Series(hyper_p_sig['number_significant_genes'], index=hyper_p_sig.index)
    hyper_p_null = hypergeometric(x, n, null_set_df, min_permutations)
    #We're only interested in printing this if in development run
    if development:
        print("Computing hypergeometric p-values for permuted values; first "+
              str(min_permutations)+" permutations...")
        hyper_p_null = hypergeometric(x, n, null_set_df, min_permutations)
        hyper_p_null.to_csv(output_dir + "GeneEnrich_"+prefix+today+'_'+str(min_permutations)+'_permut_hyperg_pval.tsv', sep='\t', index=True, index_label='gene_set')
        print("Done computing hypergeometric p-values.")

    #find empirical p-values
    print("Calculating empirical p-values.")
    part_find = partial(find_empirical_p_value, sig_val_df=sig_set_df.number_significant_genes)
    empirical_p_values = k_df.apply(part_find, axis=1)

    #Find fdr
    print("Performing FDR correction.")
    fdr = find_fdr(hyper_p_null, hyper_p_series)
    fdr = pd.Series(fdr, index=fdr.index)
    print("Done. Writing output...")


    #Prepare output
    hyper_p_series.name = 'hypergeometric_pval'
    empirical_p_values.name = 'empirical_pval'
    fdr.name = 'fdr_qval'
    permutation_counts = pd.Series(k_df.count(axis=1))
    permutation_counts.name = 'number_permutations' 

    #geneset_sig_series = pd.Series(geneset_sig_counts.number_significant_genes)
    #geneset_sig_series.index = geneset_sig_counts.gene_set
    #Merge all outputs to create uncleaned output table
    output_df = pd.DataFrame(sig_set_df, index=sig_set_df.index).merge(
                            hyper_p_series, left_index=True, right_index=True).merge(
                            empirical_p_values, left_index=True, right_index=True).merge(
                            fdr, left_index=True, right_index=True).merge(
                            permutation_counts, left_index=True, right_index=True).copy()

    output_df['fraction_significant_genes'] = output_df.number_significant_genes.divide(output_df['size'])
    output_df['bh_adjusted_pval'] = padjust_bh(output_df.empirical_pval.values)

    #accumulate fdr values
    output_df.sort_values(
                   by=['empirical_pval'], ascending=[False], inplace=True)
    output_df['fdr_qval'] = np.minimum(1, np.minimum.accumulate(np.array(output_df.fdr_qval)))
    output_df.sort_values(
                   by=['empirical_pval'], ascending=[True], inplace=True)
    #round fdr where needed
    output_df.loc[output_df.empirical_pval > output_df.fdr_qval, 'fdr_qval'] = (
                  output_df.loc[output_df.empirical_pval > output_df.fdr_qval, 'empirical_pval'])

    #find unique list of significant genes for genesets
    ###NOTE: We probably want everything annotated. ###
    ###These lines are effectively almost useless.  ###
    p_sig = output_df.loc[output_df['empirical_pval'] <= p_val_cutoff].index
    q_sig = output_df.loc[output_df['bh_adjusted_pval'] <= q_val_cutoff].index
    fdr_sig = output_df.loc[output_df['fdr_qval'] <= fdr_cutoff].index
    ###################################################
    all_output_idx = output_df.index

    #output_df.to_csv(output_dir+prefix+'_basic_output.tsv', sep='\t', index=True, index_label='gene_set')
    basic_colnames = list(output_df.columns)

    with open(output_dir+"GeneEnrich_"+prefix+today+'_log.txt', 'w') as f:
        f.write("GeneEnrich-Version 2.0.c\n")
        f.write(str(datetime.datetime.now())+'\n\n')

    output_df, p_sig_df = annotate_sig_genes(geneset_df, p_sig, all_output_idx, only_genes_df, sig_genes_df, 
                                             output_df, gene_dict, gene_colname='significant_genes',
                                             sig_colname='pass_nominal_significance')
    if len(p_sig) > 0:
        p_sig_df['pass_nominal_significance'] = False
        # test this 5 mar 2021
        # try to avoid np where
        p_sig_df["pass_nominal_significance"] = np.where(p_sig_df.gene_set.isin(list(set(p_sig))), True, False)

    p_sig_df['pass_bh_significance'] = False
    p_sig_df['pass_fdr_significance'] = False
 
    output_df, q_sig_df = annotate_sig_genes(geneset_df, q_sig, all_output_idx, only_genes_df, sig_genes_df,
                                             output_df, gene_dict, add_gene_col=False,
                                             sig_colname='pass_bh_significance')

    if len(q_sig) > 0:
        p_sig_df.loc[p_sig_df.gene_set.isin(list(set(q_sig))), 'pass_bh_significance'] = True

    output_df, fdr_sig_df = annotate_sig_genes(geneset_df, fdr_sig, all_output_idx, only_genes_df, sig_genes_df,
                                               output_df, gene_dict, add_gene_col=False,
                                               sig_colname='pass_fdr_significance')

    if len(fdr_sig) > 0:
        p_sig_df.loc[p_sig_df.gene_set.isin(list(set(fdr_sig))), 'pass_fdr_significance'] = True

    if len(mark_genes_df) > 0:
        output_df, GWAS_genes_in_ld = annotate_sig_genes(geneset_df, all_output_idx, only_genes_df, sig_genes_df,
                                                         output_df, gene_dict, gene_colname='marked_genes',
                                                         add_gene_col=True, add_sig_col=False,
                                                         mark_df = mark_genes_df, mark_columns=mark_columns)

    if len(p_sig) > 0:
        #Concat extra columns from inputs
        p_sig_df = p_sig_df.merge(concat_later, how='left', left_on='ensembl_gene_id', right_on='gene').copy()
        p_sig_df = p_sig_df.rename(columns = {'size':'total_number_genes'}) #ARH 29 Jan 2021
        p_sig_df = p_sig_df.drop(['gene'], axis = 1) #ARH 29 Jan 2021
        p_sig_df = p_sig_df.rename(columns = {'database': 'resource'}) # ARH 8 Feb 2021 changed column names
        p_sig_df = p_sig_df.loc[p_sig_df['gene_set'].isin(p_sig)] #ARH 5 Mar 2021; subsetted only to significant genesets


        if type(intron_cluster) is dict:
            p_sig_df['intron_cluster'] = p_sig_df['ensembl_gene_id'].apply(lambda x: retrieve_intron_cluster(intron_cluster, x))
            if 'variant' in p_sig_df.columns:
                p_sig_df_cols = identify_p_sig_df_columns(variant=True, intron_cluster=True)
            else:
                p_sig_df_cols = identify_p_sig_df_columns(variant=False, intron_cluster=True)
        else:
            if 'variant' in p_sig_df.columns:
                p_sig_df_cols = identify_p_sig_df_columns(variant=True, intron_cluster=False)
            else:
                p_sig_df_cols = identify_p_sig_df_columns(variant=False, intron_cluster=False)

        p_sig_df = p_sig_df[p_sig_df_cols]
        p_sig_df.to_csv(f"{output_dir}GeneEnrich_significant_gs_genes_{prefix}{today}.tsv", sep='\t', index=False)
        #Generate gene-centric table
        #In case there's inconsistent names...
        p_sig_df.rename({'ensembl_gene_id':'gene'}, inplace=True)

    gene_name_dict = {key:value.gene_name for key,value in gene_dict.items()}

    # Define genes by nominal or BH significance
    if args.sig_gs_cutoff:
        nominal_sig_genes = identify_number_unique_genes_significant_set(output_df,'pass_bh_significance')
    else:
        nominal_sig_genes = identify_number_unique_genes_significant_set(output_df, 'pass_nominal_significance')

    # gene-centric table
    if nominal_sig_genes != '0':
        print('Generating gene-centric table')
        #my_genes = list(set([gene_name_dict[x] for x in nominal_sig_genes]))
        gsc_df = curate_gene_centric_table(output_df, nominal_sig_genes, p_sig_df, gene_name_dict, sig_genes_df, resource, intron_cluster)

        gsc_df.to_csv(f"{output_dir}GeneEnrich_gene_centric_table_{prefix}{today}.tsv", sep='\t', index=False, na_rep = 'NA')
    else:
        print('No nominally significant genes. Not generating gene-centric table')        

#        mark_genes_df.to_csv(output_dir + "GeneEnrich_"+prefix+today+'_significant_gs_genes.tsv', sep='\t', index=False)
                                 
    #Write basic output ##Ruled out of final version
    #extra_colnames = [x for x in output_df.columns if any(y in x for y in ['significant', 'count'])]
    #basic_columns = basic_colnames+extra_colnames
    #output_df[basic_columns].to_csv(output_dir + prefix+'_basic_output.tsv', sep='\t', index=True, index_label='gene_set')

    #Format to scientific notation
    output_df['hypergeometric_pval'] = output_df['hypergeometric_pval'].astype(float).apply(lambda x: "{:.2e}".format(Decimal(x)) if x > 0 else 0)
    output_df['empirical_pval'] = output_df['empirical_pval'].astype(float).apply(lambda x: "{:.2e}".format(Decimal(x)) if x > 0 else 0)
    output_df['bh_adjusted_pval'] = output_df['bh_adjusted_pval'].astype(float).apply(lambda x: "{:.2e}".format(Decimal(x)) if x > 0 else 0)
    output_df['fdr_qval'] = output_df['fdr_qval'].astype(float).apply(lambda x: "{:.2e}".format(Decimal(x)))  # ARH 29 Jan 2021 replaced :4e with :2e and Decimal format

    #round to 2 decimal places
    output_df['fraction_significant_genes'] = output_df['fraction_significant_genes'].astype(float).apply(lambda x: round(x, 2))
    output_df.name = 'gene_set'
    output_df = output_df.reset_index() # ARH 8 Feb 2021 reset index
    output_df = output_df.rename(columns = {'index': 'gene_set'})
    output_df['resource'] = resource
    output_df = output_df.rename(columns = {'size': 'total_number_genes'}) # ARH 8 Feb 2021 changed name
    output_df[['resource','gene_set','total_number_genes',
               'number_significant_genes', 'fraction_significant_genes',
               'hypergeometric_pval', 'empirical_pval', 
               'bh_adjusted_pval', 'fdr_qval',
               'number_permutations', 'pass_nominal_significance',
               'pass_bh_significance', 'pass_fdr_significance', 'significant_genes']].to_csv(
         output_dir + 'GeneEnrich_gs_results_'+prefix+today+'.tsv', sep='\t', index=None, index_label='gene_set')   

    with open(output_dir+'GeneEnrich_'+prefix+today+'_log.txt', 'a') as f:
        f.write("Number of genesets below empirical_p_value threshold: "+str(len(p_sig))+'\n')
        f.write("Number of genesets below bh_adjusted value threshold: "+str(len(q_sig))+'\n')
        f.write("Number of genesets below fdr threshold: "+str(len(fdr_sig))+'\n\n')

        # number unique genes: nominal signifiant
        p_sig = identify_length_number_unique_genes_significant_set(output_df, 'pass_nominal_significance')
        f.write("Number of unique genes from significant set in genesets nominal significant: "
                +str(p_sig)+'\n')

        # number unique genes: bh signifiant
        q_sig = identify_length_number_unique_genes_significant_set(output_df, 'pass_bh_significance')
        f.write("Number of unique genes from significant set in genesets BH-adjusted significant: "
                    +str(q_sig)+'\n')

        # number unique genes: fdr signifiant
        fdr_sig = identify_length_number_unique_genes_significant_set(output_df, 'pass_fdr_significance')
        f.write("Number of unique genes from significant set in genesets fdr significant: "
                    +str(fdr_sig)+'\n')

        try:
            marked_sig = len(set(GWAS_genes_in_ld['gene']))
            f.write("Number of marked unique genes from significant set in genesets nominal significant: "
                    +str(marked_sig)+'\n\n')
        except:
            marked_sig = None
            f.write("Number of marked unique genes from significant set in genesets nominal significant: 0\n\n")

    if int(p_sig) >= 2:
        print("Building tables & creating heatmaps...")
        #Get subset with significant genesets only
        output_df.index = output_df["gene_set"]
        if args.sig_gs_cutoff:
            p_sig_df = output_df.loc[output_df['bh_adjusted_pval'].astype(float) <= float(q_val_cutoff)].copy()
        else:
            p_sig_df = output_df.loc[output_df['empirical_pval'].astype(float) <= float(p_val_cutoff)].copy()
        p_sig_geneset_df = geneset_df.loc[geneset_df.gene_set.isin(p_sig_df.index.values)].copy()
        #Make heatmaps
        try:  # ARH 8 Feb 2021 added in try/except catch
            gs_gs_plot(sig_genes_input, p_sig_geneset_df, output_dir+prefix)
            gs_gene_plot(sig_genes_input, sig_genes_input, p_sig_geneset_df, output_dir+prefix)
        except:
            print("Not enough significant genes to create heatmaps")

    keys = []
    values = []
    for val in sys.argv[1:]:
        if val.startswith('-'):
            keys.append(val)
        else:
            values.append(val)
    with open(output_dir+"GeneEnrich_"+prefix+today+'_log.txt', 'a') as f:
        for item in zip(keys, values): 
            f.write("\t".join(item)+"\n")

if __name__ == "__main__":
    start_time = time.time()
    print('Starting...')
    main()
    print("--- Done in %s seconds ---" % (time.time() - start_time))


