import numpy as np
import pandas as pd
import os
import gzip
import io

#gtf_file = '../data/gencode.v26.annotation.gtf'
#HLA_regions = '../data/HLA_region_file_hg38.tsv'

#Object containing basic information for each gene
class gene:
    def __init__(self, chr, gene_name, gene_id, gene_start, gene_finish):
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.start = int(gene_start)
        self.end = int(gene_finish)
        self.chr = chr
        self.gene_range = range(self.start, self.end)

#Object for intervals
class interval: 
    def __init__(self, chr, interval_start, interval_end):
        self.chr = chr
        self.start = int(interval_start)
        self.end = int(interval_end)
        self.interval_range = range(self.start, self.end)

#Parses gtf file, returns 'gene' objects
def parse_gtf(gtf_file, compression=None):
    """
    gtf_file: path to a GTF file, compressed or not
    """
    gene_dict = {}
    if compression == None:
        gtf = open(gtf_file)
    elif compression == 'gz':
        gtf = gzip.open(gtf_file, 'rb')
    else:
        raise Exception("Compression of gtf file not recognized. Please used "
                        ".gz or no compression.")
    with gtf: 
        for row in gtf:
            if isinstance(row, bytes):
                row = row.decode('utf-8')
            row = row.strip().split('\t')
            if row[0][0]=='#' or row[2]!='gene': continue
            attr = dict([i.split() for i in row[8].replace('"','').split(';') if i!=''])
            gene_dict[attr['gene_id'].split('.')[0]]=(gene(row[0], 
                      attr['gene_name'],attr['gene_id'].split('.')[0], row[3], row[4]))
    return(gene_dict)

#Parses interval file, returns 'interval' objects            
def parse_interval_list(interval_file): 
    """
    interval_file: path to interval file, containing 3 columns, 
    no header: chr, int_start, int_end
    """
    interval_list = []
    with open(interval_file) as file:
        for row in file:
            row = row.strip().split('\t')
            interval_list.append(interval(row[0], row[1], row[2]))
    return(interval_list)

#gene_dict = parse_gtf(gtf_file)

#interval_list = parse_interval_list(HLA_regions)

#gene_values = np.array(list(gene_dict.values()))
#for interval in interval_list:
#    chr = interval.chr
#    range_interval = set(interval.interval_range)
#    genes = [gene for gene in gene_values if gene.chr==chr]
#    print(len(genes))
#    genes = [gene for gene in genes if range_interval.intersection(gene.gene_range)]
#    print(len(genes))
