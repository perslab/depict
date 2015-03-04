#!/usr/bin/python

import pdb
import pandas as pd

# Downloaded LDL Teslovich et al. Nature 2010 data
# Extracted genome-wide significant SNPs
# 	zcat /home/data/gwas/teslovich/LDL_ONE_Eur.tbl.sorted.gz | gawk '{if($6<5e-8)print $0}' > /tmp/ldl_teslovich_nature2010.tmp

# Helper function to convert column to int
def convert(x):
        try:
                return x.astype(int)
        except:
                return x

# Add chr and pos columns
snps_df = pd.read_csv("/tmp/ldl_teslovich_nature2010.tmp", header=0, delimiter="\t")
bim = pd.read_csv("/home/data/1000G/data/phase1/bed_CEU_GBR_TSI_unrelated/CEU_GBR_TSI_unrelated.phase1_release_v3.20101123.snps_indels_svs.genotypes.bim", header=None, delimiter="\t")
bim.columns = ['Chr','snp_id','strand','Pos','allele1','allele2']
merged_df = pd.merge(snps_df, bim, left_on='SNP', right_on='snp_id', how='inner')
merged_df.Chr = merged_df.Chr.apply(convert)
merged_df.Pos = merged_df.Pos.apply(convert)
merged_df.to_csv("/home/projects/depict/DEPICT-example/ldl_teslovich_nature2010.txt", index=False, quoting=0, doublequote = False, sep='\t',columns=["SNP","Chr","Pos","P"])
