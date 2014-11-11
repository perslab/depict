#!/usr/bin/python

import os,glob

depict_jar = "/Users/tp/Google\ Drive/Projects/DEPICT/Code/DEPICT/dist/Depict.jar"
path_to_depict_folder = "/Users/tp/Google\ Drive/Projects/DEPICT/Code/DEPICT_example"
study_label = "ldl_teslovich_nature2010"
locus_file = "%s/%s.txt"%(path_to_depict_folder,study_label)
path_results = path_to_depict_folder 
ncores = 2
gene_annotation = "ensembl_v75_GRCh37_depictaffymetrix.tab"
depict_genelist_file = "ensembl_GRCh37_depictaffymetrix_genes.tab"
reconstituted_genesets_file = "%s/data/reconstituted_genesets_example.txt"%path_to_depict_folder

# Gene prioritization and reconstituted gene set enrichment analysis
os.system("java -Xms512m -Xmx16000m -XX:+UseParallelGC -XX:ParallelGCThreads=3 -jar %s %s %s %s 1 1 0 %s %s %s %s %s "%(depict_jar,"%s/data/"%path_to_depict_folder,locus_file,study_label,ncores,path_results,reconstituted_genesets_file,gene_annotation,depict_genelist_file))
