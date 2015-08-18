# Dependencies
* Mac OS X, or UNIX operating system (Microsoft Windows is not supported)
* Java SE 6 (or higher)
  * [Java.com](https://www.java.com/en/download/)
* Python version 2.7 (Python version 3 or higher is not supported)
  * [Python.org](https://www.python.org/downloads/)
* PIP (used for install Python libraries)
  * `sudo easy_install pip` 
* Python intervaltree library
  * `sudo pip install intervaltree`   
* Pandas (version 0.15.2 or higher)
  * `sudo pip install pandas`
* PLINK version 1.9 August 1 relaease (or newer)
  * [PLINK version 1.9](https://www.cog-genomics.org/plink2/) 


# DEPICT
The following description explains how to download DEPICT, test run it on example files and how to run it on your GWAS summary statistics.


## Download DEPICT
Download the compressed [DEPICT version 1 rel152](http://www.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_rel152.tar.gz) files and unzip the archive to where you would like the DEPICT tool to live on your system. Note that you when using DEPICT can write your analysis files to a different folder.  Be sure to that you meet all the dependencies described above.


## Test run DEPICT
The following steps outline how to test run DEPICT on LDL cholesterol GWAS summary statistics from [Teslovich, Nature 2010](http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html). This example is available in both the 1000 Genomes Project pilot phase DEPICT version and the 1000 Genomes Project phase 3 DEPICT version.

1. Edit `DEPICT/example/ldl_teslovich_nature2010.cfg`
  * Point `plink_executable` to where PLINK executable (version 1.9 or higher) is on our system (e.g. `/usr/bin/plink`)
2. Run DEPICT on the LDL summary statistics
  * E.g. `./src/python/depict.py example/ldl_teslovich_nature2010.cfg`
3. Investigate the results (see the [Wiki](https://github.com/perslab/DEPICT/wiki) for a description of the output format).
  * DEPICT loci `ldl_teslovich_nature2010_loci.txt`
  * Gene prioritization results `ldl_teslovich_nature2010_geneprioritization.txt`
  * Gene set enrichment results `ldl_teslovich_nature2010_genesetenrichment.txt`
  * Tissue enrichment results `ldl_teslovich_nature2010_tissueenrichment.txt`


## Run DEPICT based on your GWAS
The following steps allow you to run DEPICT on your GWAS summary statistics. We advice you to run the above LDL cholesterol example before this point to make sure that you meet all the necessary dependencies to run DEPICT.

1. Make sure that you use hg19 genomic SNP positions
2. Make an 'analysis folder' in which your trait-specific DEPICT analysis will be stored
3. Copy the template config file from `src/python/template.cfg` to your analysis folder and give the config file a more meaningful name
4. Edit your config file
  * Point `analysis_path` to your analysis folder.  This is the directory to which output files will be written
  * Point `gwas_summary_statistics_file` to your GWAS summary statistics file.  This file can be either in plain text or gzip format (i.e. having the .gz extension)
  * Specify the GWAS association p value cutoff (`association_pvalue_cutoff`). We recommend using `5e-8` or `1e-5`
  * Specify the label, which DEPICT uses to name all output files (`label_for_output_files`)
  * Specify the name of the association p value column in your GWAS summary statistics file (`pvalue_col_name`)
  * Specify the name of the marker column (`marker_col_name`). Format: <chr:pos>, ie. '6:2321'.  If this column does not exist chr_col and pos_col will be used, then leave if empty
  * Specify the name of the chromosome column (`chr_col_name`).  Leave empty if the above `marker_col_name` is set
  * Specify the name of the position column (`pos_col_name`).  Leave empty if the above `marker_col_name` is set. Please make sure that your SNP positions used human genome build GRCh37
  * Specify the separator used in the GWAS summary statistics file (`separator`). Options are
    * `tab`
    * `comma`
    * `semicolon`
    * `space`
  * Point `plink_executable` to where PLINK 1.9 executable is on yur system (e.g. `/usr/bin/plink`)
  * If you are using other genotype data than the data part of DEPICT then point `genotype_data_plink_prefix` to where your PLINK binary format 1000 Genomes Project genotype files are on your system. Specify the entire path of the filenames except the extension
5. Run DEPICT
  * `<path to DEPICT>/src/python/depict.py <path to your config file>`
6. Investigate the results which have been written to your analysis folder. See the [Wiki](https://github.com/perslab/DEPICT/wiki) for details on the output format
  * Associated loci in file ending with `_loci.txt`
  * Gene prioritization results  in file ending with `_geneprioritization.txt`
  * Gene set enrichment results  in file ending with `_genesetenrichment.txt`
  * Tissue enrichment results in file ending with `_tissueenrichment.txt`


# Troubleshooting
Please send the log file (ending with `_log.txt`) with a brief description of the problem to Tune H Pers (tunepers@broadinstitute.org).

The overall version of DEPICT follows the DEPICT publications. The current version is `v1` from [Pers, Nature Communications, 2015](http://www.nature.com/ncomms/2015/150119/ncomms6890/full/ncomms6890.html) and the release follows the number of commits of the DEPICT git repository (`git log --pretty=format:'' | wc -l`).  The latest 1000 Genomes Project pilot phase DEPICT version is `rel138`, the latest 1000 Genomes Project phase 3 version is `rel137`.

# How to cite

[Pers, Nature Communications 2015](http://www.ncbi.nlm.nih.gov/pubmed/25597830)

[1000 Genomes Project](http://www.ncbi.nlm.nih.gov/pubmed/20981092), because DEPICT makes extensively use of their data.


# Data used in these examples

LDL GWAS [summary statistics](http://csg.sph.umich.edu/abecasis/public/lipids2010/) from [Teslovich, Nature 2010](http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html) are used as input in this example. We included all SNPs with P < 5e-8 and manually added chromosome and position columns (hg19/GRCh37). 

1000 Genomes Consortium pilot release and phase 3 release data are used in DEPICT.  Please remember to cite [their paper](http://www.nature.com/nature/journal/v467/n7319/full/nature09534.html) in case you use our tool.
