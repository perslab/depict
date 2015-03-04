#! /usr/bin/python

import pandas as pd
import pdb

snps_df = pd.read_csv("/home/projects/depict/DEPICT-example/ldl_teslovich_nature2010.txt", header=0, delimiter="\t")
snps_df['marker'] = snps_df.apply(lambda x: "%s:%s"%(x.Chr,x.Pos),axis=1)
collection_df = pd.read_csv("../data/ld0.5_collection_depict_150302.txt.gz", index_col=0, header=0, delimiter="\t", compression = 'gzip')
merged_df = pd.merge(snps_df,collection_df, left_on='SNP', right_on='snp_id', how='inner')
merged_df.set_index(merged_df.marker,inplace=True)
merged_df.to_csv("/home/projects/depict/DEPICT-example/data/ld0.5_collection_depict_150302_ldl_teslovich_nature2010.txt", index=True, quoting=0, doublequote = False, sep='\t',columns=collection_df.columns)
