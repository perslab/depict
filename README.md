# Dependencies
* Mac OS X, or UNIX operating system (Microsoft Windows is not supported)
* Java SE 6 (or higher)
  * [Java.com](https://www.java.com/en/download/)
* Python version 2.7 (or higher)
  * [Python.org](https://www.python.org/downloads/)
* PIP (used for install Python libraries)
  * `sudo easy_install pip` 
* Python-bx (on Mac OS X you may be prompted to install XCode)
  * `sudo pip install bx-python`   
* Pandas (version 0.15.2 or higher)
  * `sudo pip install pandas`
* PLINK version 1.9 or higher (only needed if you want to construct loci yourself instead of using the precomputed onces for this example)
  * [PLINK version 1.9](https://www.cog-genomics.org/plink2/) 

# Examples

## Data used in these examples

LDL GWAS [summary statistics](http://csg.sph.umich.edu/abecasis/public/lipids2010/) from [Teslovich Nature 2010](http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html) are used as input in this example. We included all SNPs with P < 5e-8 and manually added chromosome and position columns (hg19/GRCh37).

## Example 1 - Run DEPICT based on the LDL cholesterol example locus file
The following steps outline how to run DEPICT based on a *precomputed LDL cholesterol DEPICT loci file*.  For a particular phenotype, the DEPICT loci file specifices which genes map to given set of associated GWAS SNPs.

1. Clone the DEPICT repository
  * `git clone git@github.com:DEPICTdevelopers/DEPICT.git`
2. Run DEPICT 
  * `python depict.py`
3. Investigate the results which have been written to the following files
  * DEPICT gene prioritization results `ldl_teslovich_nature2010_geneprioritization.txt`
  * DEPICT gene set enrichemtn results `ldl_teslovich_nature2010_genesetenrichment.txt`

## Example 2 - Run DEPICT based LDL cholesterol summary statistics
The following steps outline how to run DEPICT directly on the *LDL cholesterol summary statistics file*. (This file has been precomputed and was used directly in the above example.)

1. Clone the DEPICT repository
  * `git clone git@github.com:DEPICTdevelopers/DEPICT.git`
2. Specify in `depict.py` the path to the PLINK executable on our system
  * `plink_executable = ...` Eg. "/usr/bin/plink"
3. Retrieve the latest precomputed collection of nearest gene and gene to SNP mappings
  * Download  [LD r2 0.5 locus collection (1KG Project Phase 3 data)](http://www.broadinstitute.org/mpg/depict/depict_download/collections/ld0.5_collection_depict_150315.txt.gz)
  * Copy the collection to (do not unzip it) `ld0.5_collection_depict_150315.txt.gz DEPICT/data/collections`
5. Specify in `depict.py` the path to the new collection file
  * `collection_file = "%s/data/ld0.5_collection_depict_150315.txt.gz"%depict_path`
6. Specify in `depict.py` that you would like to clump the LDL cholesterol summary statistics, construct the DEPICT locus file and do the DEPICT analysis (except tissue enrichment analysis for which you need to download additional data as described below)
  * `step_write_plink_output = True`
  * `step_run_plink = True`
  * `step_construct_depict_loci = True`
  * `step_depict_geneprio = True`
  * `step_depict_gsea = True`
  * `step_depict_tissueenrichment = False`
7. Run DEPICT 
  * `python depict.py`
8. Investigate the results which have been written to the following files
  * DEPICT loci `ldl_teslovich_nature2010_loci.txt`
  * DEPICT gene prioritization results `ldl_teslovich_nature2010_geneprioritization.txt`
  * DEPICT gene set enrichemtn results `ldl_teslovich_nature2010_genesetenrichment.txt`

# Analyse your own GWAS summary statistics
The following steps will show you how to run DEPICT on your own GWAS summary statistics. We advice you to run example 2 to make sure that you have all the necessary parts to run a simple example.

## Setup DEPICT

Download all DEPICT files and unzip the zipped archive. Be sure to that you meet all the dependencies described above.  Instead of downloadning all files at once you can use the below steps to setup DEPICT.

1. Clone the DEPICT repository
  * `git clone git@github.com:DEPICTdevelopers/DEPICT.git`
2. Retrieve the latest precomputed collection of nearest gene and gene to SNP mappings
  * Download [LD r2 0.5 locus collection (1KG Project Phase 3 data; 249M)](http://www.broadinstitute.org/mpg/depict/depict_download/collections/ld0.5_collection_depict_150315.txt.gz) to `DEPICT/data/collections` (do not unzip it)
  * Make sure that in `depict.py` the path to the above collection file is correct
    * `collection_file = "%s/data/collections/ld0.5_collection_depict_150315.txt.gz"%depict_path`
3. Retrieve sets of precomputed background loci
  * Download [depict_backgrounds_10-400.tar.gz; 571M](http://www.broadinstitute.org/mpg/depict/depict_download/backgrounds/depict_backgrounds_10-400.tar.gz) to `DEPICT/data/`
  * Extract the zipped archive (say 'yes' to overwrite any existing files in `DEPICT/data/backgrounds/`)
    * `tar xfz depict_backgrounds_10-400.tar.gz`
4. Retrieve the reconstituted gene sets
  * Download the [reconstituted gene sets; 2.4G](http://www.broadinstitute.org/mpg/depict/depict_download/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.tgz) and extract the zipped archive to `DEPICT/data/reconstituted_genesets/`
  * Specify in `depict.py` the path to the reconstituted gene sets (set by default)
    * `reconstituted_genesets_filename = "GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary"`
5. Retrieve the tissue-specific gene expression data
  * Download the [tissue-specific gene expression data; 30M](http://www.broadinstitute.org/mpg/depict/depict_download/tissue_expression/GPL570EnsemblGeneExpressionPerTissue_DEPICT20130820_z.txt.gz) and extract the zipped archive to `DEPICT/data/tissue_expression/`
  * Specify in `depict.py` the path to the reconstituted gene sets (set by default)
    * `tissue_expression_file = "GPL570EnsemblGeneExpressionPerTissue_DEPICT20130820_z.txt"`
6. Tell DEPICT where to find PLINK and genotype data used to clump your summary statistics
  * Specify in `depict.py` the path to the PLINK executable on our system
    * `plink_executable = ...`  Eg. "/usr/bin/plink"
  * Use your own 1000 Genomes Project CEU genotype data (in binary PLINK format) or download and extract our [1000 Genomes phase 3 CEU genotypes files, 349M](http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz) to `DEPICT/data/genotype_data_plink/` ([information on data preprocessing](http://www.broadinstitute.org/mpg/snpsnap/documentation.html))
  * Specify in `depict.py` the path to genotypes. Specify the complete path and filename except the extension). See `depict.py` for an example.
    * `genotype_data_plink_prefix =  ...` 

## Run DEPICT
    
1. In `depict.py` specify the parameters related to your analysis
  * `cutoff =  ...`  E.g. "5e-8" or "1e-5", the GWAS association p value cutoff used in the DEPICT analysis
  * `label = ... `  E.g. "ldl_teslovich_nature2010", the prefix used for all output files
  * `filename = ...` Filename of your GWAS summary statistics file (can be plain text or in gzip format)
  * `pvalue_col_name = 'P'` The p value column in your GWAS summary statistics file (counting starts from 0, ie. first columns is referred to as '0'`)
  * `marker_col_name = None` The SNP identify column in your GWAS summary statistics file. Format: <chr:pos>, ie. '6:2321'.  Should be set to `None` if the below `chr_col_name` and `pos_col_name` are used
  * `chr_col_name = 'Chr'` The chromosome column in your GWAS summary statistics file. Does not need to be set if `marker_col_name` is set
  * `pos_col_name = 'Pos'` The position column in your GWAS summary statistics file. Does not need to be set if `marker_col_name` is set
  * `sep = '\t'` The separator used in your GWAS summary statistics file
2. Specify in `depict.py` that you would like to clump your summary statistics and construct the DEPICT locus file
  * `step_write_plink_output = True`
  * `step_run_plink = True`
  * `step_construct_depict_loci = True`
  * `step_depict_geneprio = True`
  * `step_depict_gsea = True`
  * `step_depict_tissueenrichment = True`
3. Run DEPICT
  * `python depict.py`
