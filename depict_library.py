#!/usr/bin/env python2.7

#use Python-2.7

import pdb,math
import pandas as pd
from bx.intervals.cluster import ClusterTree
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalTree


# Read a file with a single column
def read_single_column_file(infilename):
	with open (infilename,'r') as infile:
		return [line.strip() for line in infile.readlines()]


# Function to merge loci and return updated dataframe
def merge_loci(df):

	# Cluster the loci
	tree = ClusterTree(0,0)
	for i,marker in enumerate(df.index):
		tree.insert(df.ix[marker,'locus_gene_boundaries'][0], df.ix[marker,'locus_gene_boundaries'][1],i)

	# Create new dataframe with overlapping loci
	df_addon = pd.DataFrame()
	rows2drop = []
	for i, (start, end, overlapping_loci) in enumerate(tree.getregions()):
		if len(overlapping_loci) > 1:
			marker = ";".join(df.index[overlapping_loci])
			df_addon.ix[marker,'snp_id'] = ";".join(df.ix[overlapping_loci,'snp_id'])
			df_addon.ix[marker,'locus_start'] = start
			df_addon.ix[marker,'locus_end'] = end
			df_addon.ix[marker,'genes_in_locus'] = ";".join(set((";".join([ (str(x) if not isinstance(x, float) else '') for x in df.ix[overlapping_loci,'genes_in_locus']])).split(";")))
			df_addon.ix[marker,'nearest_gene'] = ";".join(df.ix[overlapping_loci,'nearest_gene'])
			df_addon.ix[marker,'chr'] = df.ix[overlapping_loci[0],'chr']
			rows2drop.extend(overlapping_loci)
    
	# Add merged locus and drop ovlerapping loci
	if not df.empty:
		df = df.drop(df.index[rows2drop])
		df = df.append(df_addon)
    	return df


# Function to identify and record missing SNPs
def record_missing_snps(df,df_log):
	missing_index = df.locus_start.isnull()
	if len(df.index[missing_index]) > 0:
		df_log.set_value('snps_not_found', 'Markers', ';'.join(df.index[missing_index]))
	return missing_index


# Function to identify and record HLA SNPs
def record_hla_snps(df,df_log,hla_start,hla_stop):
	hla_index = df.apply(lambda x : True if ( (x.chr == '6') and ( ( x.locus_stop > hla_start and x.locus_stop < hla_stop ) or ( x.locus_start > hla_start and x.locus_start < hla_stop ) ) ) else False, axis = 1)
	if len(df.index[hla_index]) > 0:
		df_log.set_value('snps_hla', 'Markers', ';'.join(df.index[hla_index]))
	return hla_index	


# Function to identify and record non-autosomal SNPs
def record_sex_chr_snps(df,df_log):
	non_autosomal_index = df.apply(lambda x : True if x.chr in [str(y) for y in range(1,23,1)] else False, axis = 1)
	if len(df.index[non_autosomal_index]) > 0:
		df_log.set_value('snps_non-autosomal', 'Markers', ';'.join(df.index[non_autosomal_index]))
	return non_autosomal_index 


# Function to construct gene interval tree
def get_nearest_gene_intervall_tree(depict_gene_information_file, depict_genes):
	trees = {}
	for i in range(1, 23, 1):
		trees[str(i)] = IntervalTree()
	with open (depict_gene_information_file,'r') as infile:
		for line in infile.readlines()[1:]:
			words = line.strip().split('\t')
			if words[0] in depict_genes and words[4] in [str(x) for x in range(1,23,1)]:
				tss = int(words[5]) if words[7] == '1' else int(words[6])
				trees[words[4]].insert_interval(Interval(tss, tss, value=words[0])) if words[0] in depict_genes and words[4] in [str(x) for x in range(1,23,1)] else None
	return trees


# Function to extend boundaries in case of partly overlapping genes
def get_locus_gene_boundaries(row,df_gene_boundaries):
	locus_min = row.locus_start 
	locus_max = row.locus_end
	locus_genes = row.nearest_gene.split(';')
	locus_genes.extend(row.genes_in_locus.split(';')) if not isinstance(row.genes_in_locus, float)  else None # If it empty then then it is NaN which is a float
	for gene in set(locus_genes):
		if gene in df_gene_boundaries.index:
			gene_start = df_gene_boundaries.ix[gene,'Gene Start (bp)']
			gene_end = df_gene_boundaries.ix[gene,'Gene End (bp)']
			locus_min = gene_start if gene_start < locus_min else locus_min
			locus_max = gene_end if gene_end > locus_max else locus_max
		else:
			sys.exit("%s not found in gene information file used for gene boundary computations") 
	return locus_min, locus_max


# Function to populate with DEPICT genes' TSSs (See http://pydoc.net/Python/bx-python/0.6.0/bx.intervals.intersection_tests/)
def get_nearest(chrom, pos, trees):
	chr_1_bps = 248956422
	my_max_dist = chr_1_bps / 10
	gene_up = trees[chrom].before(pos,num_intervals=1,max_dist=my_max_dist)
	gene_down = trees[chrom].after(pos,num_intervals=1,max_dist=my_max_dist)

	# Assuming there will always be a gene_down if there is no gene_up (i.e. we are at the start of the chr)
	if not gene_up:
		return gene_down[0].value #, abs(gene_down[0].start - pos)
    
	# Assuming there will always be a gene_up if there is no gene_down (i.e. we are at the end of the chr)
	elif not gene_down:
		return gene_up[0].value #, abs(gene_up[0].start - pos)

	# Test whether upstream or downstream gene
	return gene_up[0].value if abs(gene_up[0].start - pos) < abs(gene_down[0].start - pos) else gene_down[0].value #, abs(gene_up[0].start - pos) if abs(gene_up[0].start - pos) < abs(gene_down[0].start - pos) else abs(gene_down[0].start - pos)


# Helper function to convert column to int
def convert(x):
	try:
		return x.astype(int)
	except:
		return x


# Read SNPsnap collection and input SNPs from user
def write_loci(collectionfile,depict_gene_file,depict_gene_information_file,snpfile,outfile,hla_start,hla_stop):

	# Read users SNPs
	with open (snpfile,'r') as infile:
	        usersnps = [line.strip() for line in infile.readlines() ]

	# DEPICT gene universe
	depict_genes = read_single_column_file(depict_gene_file)
	
	# Read SNPsnap collection
	collection = pd.read_csv(collectionfile, index_col=0, header=0, delimiter="\t") #,compression = 'bz2')

	# Extract user SNPs
	collection_usersnps = collection.loc[usersnps,:]
	
	# Identify and remove SNPs missing	
	df_log = pd.DataFrame()
	df_index_remove = pd.DataFrame()
	index_snps_missing = record_missing_snps(collection_usersnps,df_log)
	index_snps_hla = record_hla_snps(collection_usersnps,df_log,hla_start,hla_stop)
	index_snps_non_autosomal = record_sex_chr_snps(collection_usersnps,df_log)
	collection_usersnps.drop(collection_usersnps.index[ ( index_snps_missing | index_snps_hla ) | index_snps_non_autosomal] , inplace=True)
	
	# Redefine locus boundaries based on gene boundaries
	df_gene_boundaries = pd.read_csv(depict_gene_information_file, header = 0, sep='\t', index_col = 0, usecols = [0,5,6])
	collection_usersnps['locus_gene_boundaries'] = collection_usersnps.apply(lambda x : get_locus_gene_boundaries(x, df_gene_boundaries), axis = 1)
	
	# Merge loci
	df_final = pd.DataFrame()
	for chrom in set(collection_usersnps.chr):
	    df_final = df_final.append(merge_loci(collection_usersnps[collection_usersnps.chr == chrom]))

	# Write loci and log to file
	df_final.chr = df_final.chr.apply(convert)
	df_final.locus_start = df_final.locus_start.apply(convert)
	df_final.locus_end = df_final.locus_end.apply(convert)
	df_final.to_csv(outfile, index=False, quoting=0, doublequote = False, sep='\t',columns=["snp_id","chr","locus_start","locus_end","nearest_gene","genes_in_locus"])
	df_log.to_csv(outfile.replace('.txt','.log'), index=True, quoting=0, doublequote = False, sep='\t')

	return 1
