# DEPICT-example

## Dependencies
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

## Data used in this example

LDL GWAS [summary statistics](http://csg.sph.umich.edu/abecasis/public/lipids2010/) from [Teslovich Nature 2010](http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html) are used as input in this example. We included all SNPs with P < 5e-8 and manually added chromosome and position columns (hg19/GRCh37).

## Example 1 - Run DEPICT based on the LDL cholesterol example locus file
The following steps outline how to run DEPICT based on a precomputed LDL locus file.  The locus file specifices which genes map to a given LDL associated SNP.

1. Run DEPICT 
  * `python depict.py`
2. Investigate the results which have been written to the following files
  * DEPICT loci `ldl_teslovich_nature2010_loci.txt`
  * DEPICT gene prioritization results `ldl_teslovich_nature2010_geneprioritization.txt`
  * DEPICT gene set enrichemtn results `ldl_teslovich_nature2010_genesetenrichment.txt`

## Example 2 - Run DEPICT based LDL cholesterol summary statistics
The following steps outline how to run DEPICT directly on the LDL cholesterol summary statistics file. (This file has been precomputed and was used directly in the above example.)

1. Point the `plink_binary` to the PLINK executable on your system
2. Set `step_write_plink_output`, `step_run_plink` and `step_construct_depict_loci` equal to `True`
3. Run DEPICT 
  * `python depict.py`
4. Inspect the DEPICT output, that is the files listed in the above example step #2.
