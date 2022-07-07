# GeneEnrich 

_GeneEnrich_ is a tool used to assess enrichment of a selected set of genes in predefined biological pathways or gene sets given a background set of genes expressed in a given tissue. _GeneEnrich_ is capable of controlling for gene expression levels as a potential confounding factor. This program returns an empirical gene set enrichment p-value for each pathway or gene set using the hypergeometric probability distribution, along with a Benjamini-Hochberg FDR or empirical FDR and various summary statistics. The software packages contains gene sets downloaded from MSigDB (Gene ontology, Reactome, KEGG, Hallmark gene sets), and from Mouse Genome Informatics (MGI) (mouse phenotype ontology).

_GeneEnrich_ was in particular designed to test for enrichment of target genes of eQTLs or sQTLs with top ranked GWAS p-values (e.g. P<0.05) in specific pathways, as a follow-up tool for _QTLEnrich_ when enricment of trait associations beyond genome-wide significance is found among a set of eQTLs or sQTLs in a given tissue. However, _GeneEnrich_ can be applied to any list of genes of interest. Gene sets from a number of databases are provided in the software package (Reactome and Gene Ontology from MSigDB and mouse phenotype ontology from Mouse Genome Informatics (MGI)).

Authors: John Rouhana, Andrew Hamel, written in the Lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

Date: June 2022

For questions/comments regarding this module, please contact Andrew Hamel at andrew\_hamel at meei dot harvard dot edu, or Ayellet Segre at ayellet\_segre at meei dot harvard dot edu.

## Dependencies 

_GeneEnrich_ is written in Python (developed on 3.6.4). The module versions used during development are:

  * [pandas 0.24.2](https://pandas.pydata.org/)
  * [numpy 1.16.2](https://www.scipy.org/install.html)
  * [scipy 0.17.1](https://www.scipy.org/install.html)

To export plots, _GeneEnrich_ requires two plotting libraries. The versions used during development are:

  * [matplotlib 3.0.3](https://matplotlib.org/#)
  * [seaborn 0.9.0](https://seaborn.pydata.org/index.html)

## Running Environment and Run Time

We recommend to reserve ~80GB of memory when running _GeneEnrich_.

These are estimated run times for a sample GWAS (UK Biobank Height GWAS) and eQTLs from GTEx tissues using a default of up to 100,000 permutations:

 * ~1,000 gene sets for gene ontology or reactome: 12-15 minutes

## Repository structure

``src``: the directory contains scripts for the software pipeline and for plotting with relevant instructions

``data``: the directory contains input files required to run GeneEnrich for a given set of genes. Some of the input files need to be downloaded by the user (see guidelines below)

``examples``: sample data for running GeneEnrich on a set of genes that are targets of eQTLs in a given tissue with top ranked GWAS p-values for a given trait


# Instructions on how to run GeneEnrich.py

### In order to run _GeneEnrich_, the following steps are required:

A sample command for running _GeneEnrich_ on target genes of eQTLs in GTEx artery aorta with coronary artery disease (CAD) GWAS P-value<0.05 is provided here:

```
../examples/example_geneenrich_command.sh
```

This shell script contains a sample command for running _GeneEnrich_ on multiple gene set databases available in this repo (GO, Reactome, HALLMARK, MGI), and concatenates the gene set enrichment results into a single table:

```
../examples/sample_geneenrich_command_full.sh
```

### A. Data preparation and preproceesing

1. Download gmt files
2. Curate gene sets
3. Prepare set of significant and null (background) genes
4. Prepare expression file

### B. Run _GeneEnrich_

Once input files are prepared, _GeneEnrich_ can be run using GeneEnrich.py. A sample command with input arguments is provided here:

``` 
sample_geneenrich_command.sh 
``` 

### Here is a detailed description of running each step with all possible arguments. 

## A. Data preparation and preproceesing

We provide gene sets in the appropriate format for input into GeneEnrich downloaded from MSigDB, including Gene Ontology, Reactome, KEGG, and HALLMARK gene sets, and mouse phenotype ontology gene sets parsed from Mouse Genome Informatics (MGI).
Two sets of gene sets are provided in the ../data folder, one downloaded on March 2021 and the other on June 2022. The gene sets from March 2021 were used in the manuscript referenced below.
 
If users provide their own gene sets, the following processing is necessary to make gene set files compatible as input into _GeneEnrich_. The gene set file should contain three tab-separated columns, as such:

  * Column 1 with header "database": the name of the database used. Example: KEGG, GO, etc. 
  * Column 2 with header "gene_set": the gene set/pathway to which the gene belongs
  * Column 3 with header "ensembl\_gene\_id": The Ensembl gene ID for each gene in the gene set/pathway. The Ensembl ID does not require its suffix number; if it contains the suffix number, _GeneEnrich_ will truncate it prior to analysis.

If the geneset file contains EntrezIDs instead of Ensembl IDs, **these IDs will need to be converted to Ensembl gene IDs**. Any column not labeled as above in the file will be ignored. 

Here is an example of the top rows of an input gene set (tab-delimited):

```
database gene_set                         EntrezID ensembl_gene_id
GO_BP    MITOCHONDRIAL_GENOME_MAINTENANCE 10000    ENSG00000117020
GO_BP    REPRODUCTION                     100      ENSG00000196839
```

For instructions concerning how to generate appropriate input gene set files from MSigDB, please review the following sections.

### 1. Download gmt files for curating gene sets

For this analysis, gene sets were processed from [MSigDB] (https://www.gsea-msigdb.org/gsea/msigdb/). For a given resource, such as KEGG, both `symbols` and `entrez` files are to be downloaded. 

### 2. Preprocsesing of input gene set file

If gene sets from MSigDB were downloaded, we provide scripts to prepare gene set input files for _GeneEnrich_ in the `src` directory. 

#### 2.1 Entrez-Ensembl gene ID mapping file

First, a mapping file needs to be generated between "EntrezID", "ensembl\_gene\_id", and "gene\_symbol" using the appropriate GENCODE version. 

```
./transcript_entrez_mapping.py --gencode_gtf gencode_gtf_file --entrez_gene gencode_entrez_file
```

##### Inputs to transcript\_entrez\_mapping.py to generate a mapping file

Required arguments:

* ``--gencode_gtf``: GENCODE gtf file. Can be downloaded from [here] (https://www.gencodegenes.org/human/). _Note:_ Please ensure ``--gencode_gtf`` and ``--entrez_gene`` are of the same GENCODE version. The GENCODE gtf version used in GTEx v8 was v26 'gencode.v26.annotation.gtf'.

* ``--entrez_gene``: GENCODE metadata.EntrezGene file that contains two columns: Ensembl ID and Entrez ID. Can be downloaded from [here] (https://www.gencodegenes.org/human/) under Metadata files section. Example file *gencode.v26.metadata.EntrezGene.gz*

Optional arguments: 

* ``--output_file``: Name of output file. Default is entrez\_gene\_mapping.txt

* ``--compression``: Flag to indicate if ``--gencode_gtf`` is compressed. Default is false.

We also provide a shell script to run __transcript\_entrez\_mapping.py__ and the mapping file we generated in concordance with GENCODE v26: ../data/entrez\_ensembl\_mapping\_v26\_16Jun\_2022.txt

```
./run_transcript_entrez_mapping.sh
```

#### 2.2 Preparing gene sets file

Following construction of the mapping file, gene sets can now be compiled using the following command.

```
./compile_genesets.py --symbols_file msigdb_symbols_gmt_file --entrez_file msigdb_entrez_gmt_file --entrez_mapping_file entrez_mapping_file --resource_name name_of_resource
```

##### Inputs to compile\_genesets.py to generate the input gene set file from a given database.

Required arguments:

* ``--symbols_file ``: GMT file containing list of gene sets (one per row) and the gene symbols of all genes that map to each gene set. Can be downloaded from [here] (https://www.gsea-msigdb.org/gsea/msigdb/). 

* ``--entrez_file ``: GMT file containing list of gene sets with the Entrez IDs of all genes that map to each gene set. Can be downloaded from [here] (https://www.gsea-msigdb.org/gsea/msigdb/).

* ``--resource_name``: Name of resource, e.g. Reactome.

Optional arguments: 

* ``--output_dir``: Directory to place compiled gene set. 

We also provide a shell script to run *compile\_genesets.py*.

```
./run_compile\_genesets.sh
```

Note: for the input gene sets file of gene sets curated from the Mouse Genome Informatics, please see the iPython notebook, mouse\_ontology\_parsing.ipynb.

### 3. Prepare set of significant and null genes


#### 3.1 Significant genes file

File containing a list of significant genes of interest for the gene set enrichment analysis. Only one column, "gene" is required. Two additional columns are optional, "variant" and "intron\_cluster" if assessing target genes of significant eQTLs or sQTLs. The "variant" can be the most significant e/sVariant per e/sGene and the "intron\_cluster" is related to the splicing event for the particular sQTL. The values in columns "variant" and "intron\_cluster" can be in any format as long as they are a single string. If file consists of more than one column, the file is to be tab-delimited. The "gene" column is in ENSEMBL ID format. User can include the number after the decimal for the Ensembl gene ID, though _GeneEnrich_ will remove the number after the decimal for matching with gene sets.

The significant gene list generated by _QTLEnrich_ contains target genes of eQTLs or sQTLs in a given tissue of interest with GWAS p-values below a given cutoff (default P<0.05).

Sample file of eGenes with best eVariant outputted from _QTLEnrich_:

```
gene            variant
ENSG00000272512 chr1_995982_G_A_b38
ENSG00000180758 chr1_9131513_G_T_b38
ENSG00000116786 chr1_15666843_A_G_b38

```

For sample variant with "intron\_cluster":

```
gene               variant                intron_cluster
ENSG00000227232    chr1_923955_G_A_b38    chr1:15947:16607:clu_34632:ENSG00000227232.5
ENSG00000164323    chr4_185899038_A_G_b38 chr4:185164179:185175786:clu_29218:ENSG00000164323.12
```

#### 3.2 Null/background gene list

File containing list of null genes for a given analysis. It is recommended to contain a list of all genes expressed in the tissue of interest from which the significant list of genes was taken. The file should have a single column labeled "gene". Genes need to be Ensembl gene IDs; they can contain or not the number suffix. _GeneEnrich_ will remove the number suffix if was added to ensure matching with gene set Ensembl IDs.

If this file is generated by _QTLEnrich_, then the genes in this file will contain all genes expressed in the tissue of interest excluding the significant set of genes. There is the option in _QTLEnrich_ to use only genes expressed in the given tissue with a best e/sVariant per gene that has a GWAS P-value above 0.05. This allows us to test, specifically, for the significance of the e/sQTL-effect toward possibly explaining the GWAS associations through specific pathways.

Sample null file:

```
gene
ENSG00000227232
ENSG00000268903
ENSG00000269981
```

### 4. Prepare expression file (not required if not controlling on expression in the gene set enrichment analysis)

The expression file should have, at minimum, two columns: 

  * Name: This column should contain an Ensembl gene ID for each gene in the Null and Significant gene files. The Ensembl IDs can contain or not a number after the decimal; if it does, _GeneEnrich_ will truncate the ID suffix prior to analysis. Genes not included in this list will be dropped from the Significant/Null lists. 
  * Name of your tissue: This is the name of the tissue being inspected in the analysis. This column should contain median or average expression values per gene (row) and given tissue (column). The file can contain gene expression levels for multiple tissues (columns), as long as it contains the tissue of interest.

All other columns in the file will be ignored. If using GTEx, the gene expression file can be downloaded from the [GTEx portal site](https://gtexportal.org/home/datasets). See filename: GTEx\_Analysis\_2017-06-05\_v8\_RNASeQCv1.1.9\_gene\_median\_tpm.gct.gz. Note: if using file from GTEx site, please remove first two lines and just one header row that begins with 'Name' and contains all tissue names.


### B. Command for running GeneEnrich


Here are the arguments available for running _GeneEnrich_. 

Mandatory arguments:

* --significant\_genes: Path to file with new-line separated genes of the significant set of genes. For more information, see section on significant genes above. 
* --null\_genes: Path to file with new-line separated null/backgroun genes, with optional additional columns with e/sVariant ID and intron cluster for sGenes. One column, "gene" is required. Format for "gene" is the same as for _significant\_genes_ file.
* --geneset\_file: Path to file with gene sets (pathways) to be tested.
* -gtf\_file: Path to gtf file to use for gene start and end definitions. Also used to find gene names. Be certain gtf adheres to desired human genome build version. A GENCODE gtf file from version 26 'gencode.v26.annotation.gtf' was used in GTEx release v8.

Optional arguments:

* --min\_permutations: default=1000. Minimum number of permutations for estimating an empirical gene set enrichment p-value.
* --max\_permutations: default=10000. Maximum number of permutations for estimating an empirical gene set enrichment p-value.
* --geneset\_size\_min: default=10. Minimum number of genes in gene set (pathway) expressed in given tissue (i.e. found in null or significant gene lists). Genesets with fewer genes than this number are not included in the analysis or output files.
* --geneset\_size\_max: default=1000. Maximum number of genes in gene set (pathway) expressed in given tissue (i.e. found in null or significant gene lists). Gene sets with more genes than this number are not included in analysis and removed from the output files.
* --prefix: default='GeneEnrich\_output'. Prefix for file output name. Can include a relative or absolute path.
* --fast\_permute: default=True. If False, each gene set gets its own unique gene permutation indeces. Faster if True.
* --HLA\_remove: default=False. If True, removes genes from the significant and null sets of genes, if they overlap with the HLA region defined in the HLA interval file below. This affects the gene set gene count for geneset\_size\_min and geneset\_size\_max. We tend to remove HLA region when testing a significant set of genes based on top ranked GWAS p-values, not to inflate gene set enrichment p-value due to high LD in the region.
* --HLA\_file: default='../data/HLA_region_file_hg38.tsv'. Path to file indicating HLA region. Default file for HG38 provided. Positions should be 1-based, not 0-based. 
* --bed\_remove: Path to simple BED file with 3 unlabeled columns: chr, start, end. Significant and null genes that intersect with any of these intervals will be removed. Positions should be 1-based, not 0-based. Format should mimic HLA\_file provided.
* --p\_val\_cutoff: default=0.05. Nominal empirical gene set enrichment p-value cutoff at which gene sets are considered significant.
* --q\_val\_cutoff: default=0.1. q-value cutoff at which genesets are considered significant. This is a Benjamini-Hochberg (BH) FDR value. BH FDR is computed per database.
* --fdr\_cutoff: default=0.1. Empirical FDR cutoff at which gene sets are considered significant, computed per database.
* --genes\_to\_mark\_table: Path to table with genes to mark in gene set enrichment results table. At the minimum, one column should be labeled "gene". This column should have truncated Ensembl gene IDs (withouth suffix). 
* --mark\_columns: Comma-separated column names from genes\_to\_mark\_table file to include in main output table.
* --gene\_expression\_levels: Path to a file that contains the median or average gene expression levels per gene in tissue of interest. The file can contain multiple columns for different tissues, as long as it contains the tissue of interest. If this argument is passed, _GeneEnrich_ will use expression level as a confounding factor in the permutation-based gene set enrichment analysis.
* --tissue\_of\_interest: required if -gene\_expression\_levels is passed. Tissue of interest to parse. Tissue name should be identical to header name in -gene\_expression\_levels file. _GeneEnrich_ will ignore other tissue columns in gene expression file.
* --development: Indicates development run. Determines whether to print hypergeometric null permutations table.
* --resource\_name: Name of database resource. Default is None. Note: if not included, resource name is extracted from gene set file. For example, if gene set file is _Reactome\_date.txt_, resource name is Reactome. _GeneEnrich_ extracts the prefix before the first underscore.
* --restrict\_genes\_database: Option to restrict gene-set enrichment analysis considering only genes in the input significant and null list present in the gene set database. Default is false. Include --restrict\_genes\_database in command line to set to true. Note that in cases where the gene set database contains a relatively small number of gene sets and genes (e.g., KEGG or a short custom list of gene sets), restricting to genes only found in database can lead to inflation of FDR.
* --null\_set: If true, the significant genes are not included with the null set of genes in the permutation analysis of the gene set enrichment analysis. Default is false.
* -ouput\_dir: Directory to place outputs. Defaults to working directory.
* --sig\_gs\_cutoff: Default is false. By default an empirical gene set enrichment p-value cutoff as defined by --p\_val\_cutoff is used to determine which significant gene sets and genes will be included in the heatmap and gene-centric table (described below). If true, the BH FDR cutoff as defined by --q\_val\_cutoff is used to determine the significant gene sets for plotting and the gene-centric table.


## Output Tables Format

_GeneEnrich_ outputs the following output files for each database in a separate folder. If _GeneEnrich_ is run on more than one database, all results files across can be concatenated into a single file.
Note FDR is computed per database. A description of the column headers for all of the files can be found in ../geneenrich_output_header.txt.

1. Gene set-level table with hypergeometric and empirical enrichment p-values for each gene set: (GeneEnrich\_gs\_results\_{prefix}\_{resource}\_{date}.tsv)

2. Gene-centric table with a list of significant genes enriched in at least one gene sets (cutoff determined by user) along with list of gene sets in which they are enriched or just belong to (GeneEnrich\_gene\_centric\_table\_{prefix}\_{resource}\_{date}.tsv)

3. Table with all significant gene sets that passed a nominal, Benjamini-Hochberg, or FDR enrichment significance and the significant genes that drove the enrichment (GeneEnrich\_significant\_gs\_genes\_{prefix}\_{date}.tsv)

4. Log file that contains that _GeneEnrich_ command run, and summary statistics of number of significantly enriched gene sets and number of significant genes in at least one significantly enriched gene set given different significance cutoffs (GeneEnrich\_{prefix}\_{resource}\_{date}\_log.txt)

## Output Plots

GeneEnrich plots two heatmaps for each run.

1. Heatmap clustering gene sets by significant genes {prefix}\_{resource}\_gs\_gene\_intersection.pdf
2. Heatmap clustering gene sets by gene sets based on overlap of significant genes driving each gene set's enrichment {prefix}\_{resource}\_gs\_gene_intersection.pdf
The colorbar represents the fraction of significant genes in a given gene set in the row that overlaps with the significant genes in the gene set in the column.
Hierarchical clustering was performed on both plots on rows and columns using the euclidean distance between fraction of overlapping genes.

## Method described in manuscript:

Andrew R. Hamel†, John M. Rouhana†, Alvaro N. Barbeira, Hae Kyung Im, Ayellet V. Segrè. Tissue-dependent expression and splicing QTL and gene-set enrichment analysis of complex trait associations uncovers disease mechanisms Submitted 2022.

