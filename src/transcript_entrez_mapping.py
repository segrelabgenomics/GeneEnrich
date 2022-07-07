#!/usr/bin/env python3

"""
Author: Andrew Hamel
Affiliation: Massachusetts Eye and Ear, Harvard Medical School
Date: June 2022

This script generates a mapping file between EntrezID, ensembl_gene_id (without decimal version),
and gene symbol (HGNC) for all genes
output file is tab-delimited

Input files can be downloaded from GENCODE website
https://www.gencodegenes.org/human/

"""

import sys
import os
import os.path
import gzip
import re
import argparse
import pandas as pd

from pathlib import Path
from argparse import Namespace, ArgumentParser
from pandas import DataFrame
from typing import List, Dict, Tuple


def parse_args() -> Namespace:
    """
    arguments to run script
    """
    parser = ArgumentParser(description="Curate entrez mapping file.")

    parser.add_argument("--gencode_gtf",
                        type=str, required=True,
                        help="GENCODE gtf file.")
    parser.add_argument("--compression",
                        action="store_true",
                        help="Indicates if --gencode_gtf is compressed.")
    parser.add_argument("--entrez_gene",
                        type=str, required=True,
                        help="GENCODE metadata.EntrezGene file.")
    parser.add_argument("--output_file",
                        type=str, default="entrez_gene_mapping.txt",
                        help="Name of output file.")
    return parser.parse_args()


def die(epitaph: str,
        exit_code: int = 1) -> None:
    """
    In event of error, prints message and exits
    Default exit-code is 1
    """
    print(epitaph)
    sys.exit(exit_code)


def verify_files_exist(args: Namespace) -> bool:
    """
    Ensure specific files exist
    """
    files_exist = [args.gencode_gtf,
                   args.entrez_gene]

    return all(
        [Path(f).is_file() for f in files_exist]
        )


def create_df_lists(gene_list: List[List[str]]) -> DataFrame:
    """
    creates df from list of lists
    """
    df = pd.DataFrame(gene_list,
                      columns = ["ensembl_gene_id_full", "gene_symbol",
                                 "gene_type", "transcript_id_full", "chr", "start",
                                 "stop", "strand"])
    df["ensembl_gene_id"] = df["ensembl_gene_id_full"].apply(lambda x: x.split(".")[0])    
    df["transcript_id"] = df["transcript_id_full"].apply(lambda x: x.split(".")[0])

    df = df[["ensembl_gene_id", "ensembl_gene_id_full",
             "transcript_id", "transcript_id_full","gene_symbol", "gene_type", "chr", "start",
             "stop", "strand"]]

    return df

def parse_gencode_file(file: str,
                       compression: bool = False) -> DataFrame:
    """
    parses file and selects 
    gene name, gene id, start, stop
    """
    gene_list = []

    if compression:
        gtf = gzip.open(file, "rb")
    else:
        gtf = open(file, "r")

    with gtf:
        for row in gtf:
            if isinstance(row, bytes):
                row = row.decode("utf-8")
            row = row.strip().split("\t")

            if row[0][0] == "#" or row[2] != "transcript":
                continue

            chrom, start, stop, strand = row[0], row[3], row[4], row[6]
            # creates dictionary from row[8]
            # attribute 
            attr = dict([i.split() for i in row[8].replace('"','').split(';') if i!=''])
            line = [attr["gene_id"], attr["gene_name"], 
                    attr["gene_type"], attr["transcript_id"], 
                    chrom, start, 
                    stop, strand]        
            gene_list.append(line)

    df = create_df_lists(gene_list)

    return df


def parse_entrez_metadata_file(file: str) -> DataFrame:
    """
    parses gencode.vX.metadata.EntrezGene.gz
    """
    df = pd.read_csv(file, sep="\t", 
                     header=None, names=["transcript_id_full", "EntrezID"])
    df["transcript_id"] = df["transcript_id_full"].apply(lambda x: x.split(".")[0])

    return df


def curate_mapping(gencode_df: DataFrame,
                   entrez_df: DataFrame) -> DataFrame:
    """
    merges dataframe
    extracts unique combinations
    """
    df = pd.merge(gencode_df, entrez_df,
                         on=["transcript_id", "transcript_id_full"])
    df = df[["EntrezID", "ensembl_gene_id", "gene_symbol"]]
    df = df.drop_duplicates(["EntrezID", "ensembl_gene_id", "gene_symbol"])
    df["EntrezID"] = pd.to_numeric(df["EntrezID"], downcast="integer")
    df = df.sort_values(["EntrezID"])

    return df


def main() -> None:
    """
    main section
    """
    args = parse_args()

    if not verify_files_exist(args):
        die(f"One or more required files "
            f"do not exist. Please fix.", 255)

    print("parsing gencode gtf")
    gencode_df = parse_gencode_file(args.gencode_gtf,
                                    args.compression)

    print("parsing entrez metadata file")
    entrez_df = parse_entrez_metadata_file(args.entrez_gene)

    mapping_df = curate_mapping(gencode_df, entrez_df)
    mapping_df.to_csv(args.output_file,
                      sep="\t", header=True, index=None)

if __name__ == "__main__":
    main()






