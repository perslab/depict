#! /usr/bin/python

import pandas as pd
import gzip,sys
from scipy.stats import ranksums
#from scipy.stats import ttest_ind
#import pdb
#import numpy

outfile_prefix = sys.argv[1]
start_gs = int(sys.argv[2])
last_gs = int(sys.argv[3])

# The below files can be downloaded from the DEPICT website
gene_sets_discrete = "/cvar/jhlab/tp/depict_download/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355_term_gene_assignments.txt.gz"
reconstituted_gene_sets_df = pd.read_csv("/cvar/jhlab/tp/depict_download/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz",compression="gzip",index_col=0,sep="\t")

def compute_ranksum_p(start_gs,last_gs):
	res = {}
	with gzip.open(gene_sets_discrete,"r") as infile:
	    gs = infile.readlines()
	    for line in gs[start_gs:min(last_gs,len(gs))]:
	        words = line.strip().split("\t")
	        gs_genes = words[2].split("|")
	
	        # Stratify 
	        gs_gene_scores = []
	        other_genes_scores = []
	        for g in reconstituted_gene_sets_df.index:
	            if g in gs_genes:
	                gs_gene_scores.append(reconstituted_gene_sets_df.iloc[reconstituted_gene_sets_df.index.get_loc(g),reconstituted_gene_sets_df.columns.get_loc(words[0])])
	            else:
	                other_genes_scores.append(reconstituted_gene_sets_df.iloc[reconstituted_gene_sets_df.index.get_loc(g),reconstituted_gene_sets_df.columns.get_loc(words[0])])    
	
	        # Test
	        z, p1 = ranksums(gs_gene_scores, other_genes_scores)
	        t, p2 = ttest_ind(gs_gene_scores, other_genes_scores, equal_var=False)
	        print "{}: gs_median={}, other_median={}, p_utest={} p_ttest={})".format(words[0],numpy.median(gs_gene_scores),numpy.median(other_genes_scores),p1,p2)
	        res[words[0]] = p1
	    
	    # Write to file
	    with open("{}_{}_{}.tab".format(outfile_prefix,start_gs,last_gs), "w") as f:
	    	for gs in res:
	       		f.write("{}\t{}\n".format(gs,res[gs]))


compute_ranksum_p(start_gs,last_gs)
