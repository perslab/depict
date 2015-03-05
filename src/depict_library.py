#!/usr/bin/env python2.7

import pdb,math,sys,gzip,os
import pandas as pd
from glob import glob
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
def record_missing_snps(df,log_dict):
	missing_index = df.locus_start.isnull()
	if len(df.index[missing_index]) > 0:
		log_dict['snps_not_found'] = ';'.join(df.index[missing_index])
	return missing_index,log_dict


# Function to identify and record HLA SNPs
def record_hla_snps(df,log_dict,hla_start,hla_stop):
	hla_index = df.apply(lambda x : True if ( (x.chr == '6') and ( ( x.locus_stop > hla_start and x.locus_stop < hla_stop ) or ( x.locus_start > hla_start and x.locus_start < hla_stop ) ) ) else False, axis = 1)
	if len(df.index[hla_index]) > 0:
		log_dict['snps_hla'] = ';'.join(df.index[hla_index])
	return hla_index, log_dict


# Function to identify and record non-autosomal SNPs
def record_sex_chr_snps(df,log_dict):
	non_autosomal_index = df.apply(lambda x : True if x.chr in [str(y) for y in range(1,23,1)] else False, axis = 1)
	if len(df.index[non_autosomal_index]) > 0:
		log_dict['snps_non-autosomal'] = ';'.join(df.index[non_autosomal_index])
	return non_autosomal_index, log_dict


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
def construct_depict_loci(analysis_path,label,cutoff,collectionfile,depict_gene_file,depict_gene_information_file,outfile,hla_start,hla_stop):

	print "\nWriting DEPICT input loci"

	# Dictionary with log information
	log_dict = {}

	# Read users SNPs
	user_snps = get_plink_index_snps(analysis_path,label,cutoff)

	# DEPICT gene universe
	depict_genes = read_single_column_file(depict_gene_file)
	
	# Read SNPsnap collection
	collection = pd.read_csv(collectionfile, index_col=0, header=0, delimiter="\t", compression = 'gzip')

	# Extract user SNPs
	collection_usersnps = collection.loc[user_snps,:]
	
	# Identify and remove SNPs missing	
	df_index_remove = pd.DataFrame()
	index_snps_missing, log_dict = record_missing_snps(collection_usersnps,log_dict)
	index_snps_hla, log_dict = record_hla_snps(collection_usersnps,log_dict,hla_start,hla_stop)
	index_snps_non_autosomal, log_dict = record_sex_chr_snps(collection_usersnps,log_dict)
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

	return log_dict


# Helper function to save SNPs that are in my data
def write_plink_input(path,filename,label,marker_col,p_col,chr_col,pos_col,sep,genotype_data_plink_prefix):

	print "\nReading user input and writing PLINK output"

	# Read mapping (Faster than usgin pandas dataframe, because constructing a chr:pos index takes a hell of time
	infile = open("%s.bim"%genotype_data_plink_prefix,'r')
	lines = infile.readlines()
	infile.close()
	mapping = {}
	for line in lines:
		words = line.strip().split()
		mapping["%s:%s"%(words[0],words[3])] = words[1] 

	# Read users SNP file
	with ( gzip.open("%s/%s"%(path,filename),'r') if '.gz' in filename else open("%s/%s"%(path,filename),'r') ) as infile, open("%s/%s.tab"%(path,label),'w') as outfile:
		outfile.write("SNP_chr_pos\tSNP\tChr\tPos\tP\n")
		for line in infile.readlines()[1:]:
			words = line.strip().split(sep)
			if marker_col is not None:
				chrom = int(words[marker_col].split(":")[0])
				pos = int(words[marker_col].split(":")[1])
				marker_id = words[marker_col] 
			elif chr_col is not None and pos_col is not None:
				chrom = int(words[chr_col])
				pos = int(words[pos_col])
				marker_id = "%s:%s"%(chrom,pos)
			else:
				sys.exit('Please specify either a column for the marker (format <chr:pos>) or the chromosome and position columns.')	

			if marker_id in mapping:
				outfile.write("%s\t%s\t%s\t%s\t%s\n"%(marker_id,mapping[marker_id],chrom,pos,words[p_col]))
			else:
				outfile.write("%s\t-\t%s\t%s\t%s\n"%(marker_id,chrom,pos,words[p_col]))
	return 0
	

# Helper function to run PLINK
def run_plink(path, label, genotypes_1kg, plink_binary, plink_extra_params, cutoff, distance, r2): 

	print "\nRunning PLINK"

	plink_prefix = "%s --bfile %s %s"%(plink_binary,genotypes_1kg,plink_extra_params)
	cmd = "%s --clump-p1 %s --clump-kb %s --clump-r2 %s --clump %s/%s.tab --out %s/%s"%(plink_prefix,cutoff,distance,r2,path,label,path,label) 
	return os.popen(cmd).readlines()


# Helper function to retrieve PLINK Index SNPs
def get_plink_index_snps(path,label,cutoff):

	# Read file with correct SNP indices
	id_df = pd.read_csv("%s/%s.tab"%(path,label),index_col=0,header=0,sep="\t")

	# Read PLINK results
	index_snp_col = 7
	index_snps = []
	with open("%s/%s.clumped"%(path,label),'r') as infile:
		lines = infile.readlines()[1:]
		for line in lines:
			words = line.strip().split(' ') 
			if len(words) >= index_snp_col-1: 
				index_snps.append(words[index_snp_col]) if words[index_snp_col] != '' else None

	# Re-map to chr:pos
	return id_df.index[id_df.SNP.isin(index_snps)]

