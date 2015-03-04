# DEPICT-example

## Dependencies
* Java version 8 (or higher)
  * [Java.com](https://www.java.com/en/download/)
* Python version 2.7 (or higher)
  * [Python.org](https://www.python.org/downloads/)
* Python-bx
  * `sudo pip install bx-python`   
* Pandas
  * `sudo pip install pandas`

## Data used in this example

LDL GWAS [summary statistics](http://csg.sph.umich.edu/abecasis/public/lipids2010/) from [Teslovich Nature 2010](http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html) are used as input in this example. We included all SNPs with P < 5e-8 and manually added chromosome and position columns (hg19/GRCh37).

## Quick start
* Modify depict.py
  * Change the `analysis_path` to the path of the `DEPICT-example` folder
* Run DEPICT
  * `python depict.py`
