import argparse
import numpy as np
import pandas as pd
import re


"""
Tools used in GeneEnrich 

Authors: John Rouhana, Andrew Hamel
Affiliation: Massachusetts Eye and Ear, Harvard Medical School
Date: June 2022

"""

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

#BH adjustment
def padjust_bh(p_values):
    """
    Does the same thing as  p.adjust(p, method="BH") from R
    p = np.array of p-values
    """
    n = len(p_values)
    idx = np.arange(n,0,-1) #ranks in descending largeness
    des_ind = np.argsort(p_values)[::-1] #indices in descending largeness
    orig_order = np.argsort(des_ind) #The original order in which the p values are passed
    return np.minimum(1, np.minimum.accumulate(np.float(n)/idx * np.array(p_values)[des_ind]))[orig_order]

#turns a dataframe to a symbol delimited string
def df_to_delim(df, col_delim='|', row_delim=';'):
    """
    Note: drops headers and indices
    df: pandas DataFrame object to turn into string
    col_delim: delimiter for values between columns
    row_delim: delimiter for values between rows
    """ 
    x = df.to_string(header=False, index=False, index_names=False).split('\n')
    x = row_delim.join([re.sub("\s+", col_delim, v.lstrip().strip()) for v in x])
    return(x)

#gets bin values for a pandas series
def get_bins(series, num_quantiles=10.):
    """
    Find the percent of data in each quantile in order to appropriately bin data by bin number
    """
    num_quantiles = float(num_quantiles)
    percent_data_per_quantile = 1./num_quantiles
    return(series.rank(pct = True).apply(lambda x: x/percent_data_per_quantile).apply(np.ceil).astype(np.uint8))
    
