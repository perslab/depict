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
* Pandas
  * `sudo pip install pandas`
* PLINK (only needed if you want to construct loci yourself instead of using the precomputed onces for this example)
  * [PLINK version 1](http://pngu.mgh.harvard.edu/~purcell/plink/) or [PLINK version 2](https://www.cog-genomics.org/plink2/) 

# DEPICT-example

## Data used in this example

LDL GWAS [summary statistics](http://csg.sph.umich.edu/abecasis/public/lipids2010/) from [Teslovich Nature 2010](http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html) are used as input in this example. We included all SNPs with P < 5e-8 and manually added chromosome and position columns (hg19/GRCh37).

## Example 1 - Run DEPICT based on the LDL cholesterol example locus file
The following steps outline how to run DEPICT based on a *precomputed LDL cholesterol DEPICT loci file*.  For a particular phenotype, the DEPICT loci file specifices which genes map to given set of associated GWAS SNPs.

1. Run DEPICT 
  * `python depict.py`
2. Investigate the results which have been written to the following files
  * DEPICT gene prioritization results `ldl_teslovich_nature2010_geneprioritization.txt`
  * DEPICT gene set enrichemtn results `ldl_teslovich_nature2010_genesetenrichment.txt`

## Example 2 - Run DEPICT based LDL cholesterol summary statistics
The following steps outline how to run DEPICT directly on the *LDL cholesterol summary statistics file*. (This file has been precomputed and was used directly in the above example.)

1. Specify in `depict.py` the path to the PLINK executable on our system
  * `plink_binary = ...` # Eg. "/usr/bin/plink"
2. Download the latest precomputed collection of nearest gene and gene to SNP mappings
  * [LD r2 0.5 locus collection (1KG Project Phase 1 data)](http://www.broadinstitute.org/mpg/depict/depict_download/collections/ld0.5_collection_depict_150302.txt.gz)
3. Copy the collection to (do not unzip it)
  * `cp ld0.5_collection_depict_150302.txt.gz DEPICT/data/`
3. Specify in `depict.py` the path to the new collection file
  * `collection_file = "%s/data/ld0.5_collection_depict_150302.txt.gz"%depict_path`
4. Specify in `depict.py` that you would like to clump the LDL cholesterol summary statistics and construct the DEPICT locus file
  * `step_write_plink_output = True`
  * `step_run_plink = True`
  * `step_construct_depict_loci = True`
5. Run DEPICT 
  * `python depict.py`
6. Investigate the results which have been written to the following files
  * DEPICT loci `ldl_teslovich_nature2010_loci.txt`
  * DEPICT gene prioritization results `ldl_teslovich_nature2010_geneprioritization.txt`
  * DEPICT gene set enrichemtn results `ldl_teslovich_nature2010_genesetenrichment.txt`

# Getting started on your own GWAS summary statistics
The following steps will show you how to run DEPICT on your own GWAS summary statistics.  If Example 2 worked, the following steps should be straightforward.

1. Retrieve the latest precomputed collection of nearest gene and gene to SNP mappings
  * Download [LD r2 0.5 locus collection (1KG Project Phase 1 data; 249M)](http://www.broadinstitute.org/mpg/depict/depict_download/collections/ld0.5_collection_depict_150302.txt.gz) to `DEPICT/data/` (do not unzip it)
  * Make sure that in `depict.py` the path to the above collection file is correct
    * `collection_file = "%s/data/ld0.5_collection_depict_150302.txt.gz"%depict_path`
2. Retrieve sets of precomputed background loci
  * Download [depict_backgrounds_10-400.tar.gz; 571M](http://www.broadinstitute.org/mpg/depict/depict_download/backgrounds/depict_backgrounds_10-400.tar.gz) to `DEPICT/data/`
  * Extract the zipped archive (say 'yes' to overwrite any existing files in `DEPICT/data/backgrounds/`)
    * `tar xfz depict_backgrounds_10-400.tar.gz`
3. Retrieve reconstituted gene sets
  * Download [;M]() and extract the zipped archive to `DEPICT/data/`
  * Specify in `depict.py` the path to the reconstituted gene sets
    * `reconstituted_genesets_file = "%s/data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary"%depict_path` 
4. Tell DEPICT where to find tools/data for clumping our summary statistics
  * Specify in `depict.py` the path to the PLINK executable on our system
    * `plink_binary = ...` # Eg. "/usr/bin/plink"
  * Download our [1000 Genomes Project PLINK binary-formatted genotypes (CEU super population), 405M](http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz) to a directory on your system ([information on data preprocessing](http://www.broadinstitute.org/mpg/snpsnap/documentation.html))
  * Extract the zipped archive
  * Specify in `depict.py` the path to genotypes. Specific the complete path and filename (except the extension, see `depict_example.py` for an example.
    * `genotype_data_plink_prefix =  ...` 
5. In depict.py modify the following parameters
  * `cutoff =  ...` # E.g. "5e-8"
  * `label = ... ` # E.g. "ldl_teslovich_nature2010"
  * `filename_extension = ...` # E.g. ".txt"
  * `pvalue_col = 3` # NB: Counting starts from 0, ie. first columns is referred to as '0'`
  * `marker_col = None` # Format: <chr:pos>, ie. '6:2321'.  Should be set to `None` if the below `chr_col` and `pos_col` are used
  * `chr_col = 1` # Does not need to be set if `marker_col` is set
  * `pos_col = 2` # Does not need to be set if `marker_col` is set
  * `sep = '\t'`
6. Specify in `depict.py` that you would like to clump your summary statistics and construct the DEPICT locus file
  * `step_write_plink_output = True`
  * `step_run_plink = True`
  * `step_construct_depict_loci = True`
