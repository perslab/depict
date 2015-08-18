#!/usr/bin/env python2.7

import subprocess
import pandas as pd

def run_plink_clumping(plink_binary, plink_genotype_data_plink_prefix, plink_clumping_pvalue, plink_clumping_distance, plink_clumping_r2, path, gwas_filename, output_filename):
	plink_clumping_snp_column = "SNP" 
	plink_clumping_pvalue_column = "P"
	"""
	Function to run PLINK
	"""
	#"--exclude", "/tmp/mylist.txt",
	cmd = [plink_binary, 
		"--bfile", plink_genotype_data_plink_prefix,
		"--clump-p1", str(plink_clumping_pvalue),
		"--clump-kb", str(plink_clumping_distance),
		"--clump-r2", str(plink_clumping_r2),
		"--clump-snp-field", plink_clumping_snp_column,
		"--clump-field", plink_clumping_pvalue_column,
		"--clump", path + "/" + gwas_filename,
		"--out", path + "/" + output_filename
	]
	return subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
