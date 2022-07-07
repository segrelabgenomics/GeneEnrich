#!/usr/bin/env python3

"""
Author: Andrew Hamel
Affiliation: Massachusetts Eye and Ear, Harvard Medical School
Date: June 2022

This script curates a geneset
to be used when running GeneEnrich

Input files can be downloaded from GENCODE website
example symbols and entrez files downloaded from MSigDB
for REACTOME
https://www.gsea-msigdb.org/gsea/msigdb/

"""

import sys
import os
import os.path
import argparse
import pandas as pd

from pathlib import Path
from typing import List
from pandas import DataFrame
from argparse import Namespace, ArgumentParser

def parse_args() -> Namespace:
    """
    arguments to run script
    """
    parser = ArgumentParser(description="Curate genesets for GeneEnrich")

    parser.add_argument("--symbols_file",
                        type=str, required=True,
                        help="GMT file containing genesets their gene symbols.")
    parser.add_argument("--entrez_file",
                        type=str, required=True,
                        help="GMT file containing genesets their entrezids.")
    parser.add_argument("--entrez_mapping_file",
                        type=str, required=True,
                        help="File containing gene symbol and associated ensembl ids.")
    parser.add_argument("--resource_name",
                        type=str, required=True,
                        help="Name of resource")
    parser.add_argument("--output_dir",
                        type=str, help="Directory to place output files.")

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
    files_exist = [args.symbols_file,
                   args.entrez_file,
                   args.entrez_mapping_file]

    return all(
        [Path(f).is_file() for f in files_exist]
        )


def verify_outdir(args: Namespace) -> str:
    """
    checks if output directory is provided
    ensures it is a valid output directory
    
    if none provided, returns empty string and 
    default will be current working directory
    """
    if args.output_dir:
        if Path(args.output_dir).is_dir():
            return args.output_dir+"/"
        else:
            die(f"{output_dir} does not exist. Please fix.", 1)
    return ""


def parse_gmt(file: str, col_name: str) -> DataFrame:
    """
    opens gmt file
    creates dictionary 
    geneset: [list of genes]
    split[0] is geneset of splitted line
    """
    gene_dict = {}
    with open(file, "r") as f:
        lines = f.readlines()
        for line in lines:
            split = line.split()
            gene_dict[split[0]] = split[2:]  

    # create dataframe from dictionary
    df = pd.DataFrame(gene_dict.items(),
                      columns=["gene_set", col_name])
    df = df.explode(col_name)
    df = df.reset_index(drop=True)

    return df


def parse_entrez_mapping(file: str) -> DataFrame:
    """
    reads file   
    """
    return pd.read_csv(file, sep="\t")


def curate_geneset(symbols_file: str,
                   entrez_file: str,
                   entrez_mapping_file: str,
                   resource: str) -> DataFrame:
    """
    returns dataframe of geneset file for GeneEnrich
    columns are
    1. database, e.g. REACTOME
    2. gene_set
    3. EntrezID
    4. ensembl_gene_id, without decimal version
    """
    symbols_df = parse_gmt(symbols_file, "gene_symbol")
    entrez_df = parse_gmt(entrez_file, "EntrezID")
    entrez_mapping_df = parse_entrez_mapping(entrez_mapping_file)

    # join dataframes 
    full_df = pd.concat([symbols_df, entrez_df], axis=1)
    full_df["EntrezID"] = pd.to_numeric(full_df["EntrezID"], 
                                        downcast="integer")
    full_df.columns = ["gene_set", "gene_symbol", "gene_set2", "EntrezID"]

    # fix geneset
    full_df["gene_set"] = full_df["gene_set"].apply(lambda x: "_".join(x.split("_")[1:]))

    # merge with entrez_mapping df to get ensembl gene ids
    full_df = pd.merge(full_df, entrez_mapping_df, on=["gene_symbol", "EntrezID"])
    full_df["database"] = resource
    full_df = full_df[["database", "gene_set", "EntrezID", "ensembl_gene_id"]]
    full_df = full_df.sort_values(["gene_set"], ascending=True)

    return full_df

def main() -> None:
    """
    main section
    """
    args = parse_args()

    print("verifying input files")
    if not verify_files_exist(args):
        die(f"One or more required files "
            f"do not exist. Please fix.", 255)

    print("verifying output directory")
    output_dir = verify_outdir(args)

    print("curating geneset df")
    geneset_df = curate_geneset(args.symbols_file,
                                args.entrez_file,
                                args.entrez_mapping_file,
                                args.resource_name)

    print(f"writing geneset df to {output_dir}")
    geneset_df.to_csv(f"{output_dir}{args.resource_name}_genesets.txt",
                      header=True, index=None, sep="\t")

if __name__ == "__main__":
    main()




