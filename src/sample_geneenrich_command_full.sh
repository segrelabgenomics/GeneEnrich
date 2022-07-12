#!/bin/bash

# Author: Andrew Hamel
# Affiliation: Massachusetts Eye and Ear, Harvard Medical School
# Date: June 2022 
#
# This is a shell script 
# with a sample run of GeneEnrich
# for all available gene sets provided

# List of significant genes of interest
significant_file="../examples/artery_aorta_cad_sig_genes.txt"
# Background list of genes expressed in given tissue
null_file="../examples/artery_aorta_cad_null_genes.txt"
# Downloaded from GENCODE
gtf_file=""
# hla file
hla_file="../data/HLA_region_file_hg38.tsv"
HLA_remove=True
# prefix
prefix=""
# gene set directory
gene_setdir="../data/gene_sets_June2022/"
# fill in output directory
output_dir=""
mkdir -p $output_dir

# creates directories for each resource
mkdir -p $output_dir/'kegg'
mkdir -p $output_dir/'go_bp'
mkdir -p $output_dir/'go_mf'
mkdir -p $output_dir/'go_cc'
mkdir -p $output_dir/'reactome'
mkdir -p $output_dir/'mgi'
mkdir -p $output_dir/'hallmark'
mkdir -p $output_dir/'all'

# GO MF
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"GO_MF_genesets.txt" \
   --prefix $prefix'_go_mf' \
   --output_dir $output_dir/'go_mf/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'go_mf' \
   --restrict_genes_database

out_file=$output_dir/'go_mf/GeneEnrich_gs_results_'$prefix"*".tsv
if [ -f $out_file ]; then
  cat $out_file1 > $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'go_mf/GeneEnrich_gene_centric_table_'$prefix"*".tsv
if [ -f $out_file ]; then
  cat $out_file > $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

# GO BP
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"GO_BP_genesets.txt" \
   --prefix $prefix'_go_bp' \
   --output_dir $output_dir/'go_bp/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'go_bp' \
   --restrict_genes_database

out_file=$output_dir/'go_bp/GeneEnrich_gs_results_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'go_bp/GeneEnrich_gene_centric_table_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

# GO CC
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"GO_CC_genesets.txt" \
   --prefix $prefix'_go_cc' \
   --output_dir $output_dir/'go_cc/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'go_cc' \
   --restrict_genes_database

out_file=$output_dir/'go_cc/GeneEnrich_gs_results_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'go_cc/GeneEnrich_gene_centric_table_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

# KEGG
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"KEGG_genesets.txt" \
   --prefix $prefix'_kegg' \
   --output_dir $output_dir/'kegg/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'kegg' \
   --restrict_genes_database

out_file=$output_dir/'kegg/GeneEnrich_gs_results_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'kegg/GeneEnrich_gene_centric_table_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

# REACTOME
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"REACTOME_genesets.txt" \
   --prefix $prefix'_reactome' \
   --output_dir $output_dir/'reactome/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'reactome' \
   --restrict_genes_database

out_file=$output_dir/'reactome/GeneEnrich_gs_results_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'reactome/GeneEnrich_gene_centric_table_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

# MGI
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"MGI_genesets.txt" \
   --prefix $prefix'_mgi' \
   --output_dir $output_dir/'mgi/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'mgi' \
   --restrict_genes_database

out_file=$output_dir/'mgi/GeneEnrich_gs_results_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'mgi/GeneEnrich_gene_centric_table_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

# HALLMARK
python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_setdir"HALLMARK_genesets.txt" \
   --prefix $prefix'_hallmark' \
   --output_dir $output_dir/'hallmark/' \
   --min_permutations 1000 --max_permutations 100000 \
   --HLA_remove $hla_remove \
   --gtf_file $gtf_file \
   --HLA_file $hla_file \
   --resource_name 'hallmark' \
   --restrict_genes_database

out_file=$output_dir/'hallmark/GeneEnrich_gs_results_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gs_results_'$prefix'_concatenated.tsv'
fi

out_file=$output_dir/'hallmark/GeneEnrich_gene_centric_table_'$prefix"*".tsv

if [ -f $out_file ]; then
  cat $out_file | tail -n+2 >> $output_dir/'all/GeneEnrich_gene_centric_table_'$prefix'_concatenated.tsv'
fi

sleep 30

