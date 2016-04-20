#/usr/bin/env python2.7

import warnings
warnings.filterwarnings("ignore")
import pdb,math,sys,gzip,os,re,subprocess
import pandas as pd
import glob
from time import time
from intervaltree import Interval, IntervalTree

# Read a file with a single column
def read_single_column_file(infilename):
	with open (infilename,'r') as infile:
		return [line.strip() for line in infile.readlines()]

# Function to translate separator
def get_separator(sep_symbol):
	sep = None
	if sep_symbol == 'tab':
		sep = '\t'
	elif sep_symbol == 'space':
		sep = ' '
	elif sep_symbol == 'comma':
		sep = ','
	elif sep_symbol == 'semicolon':
		sep = ';'
	else:
		raise Exception("Please make sure your GWAS summary statistics file separator is specified correctly")
	return sep
	

# Function to build a nonoverlapping interval tree
def insert_interval(tree, begin, end, value):
	overlap = tree[begin:end]
	if len(overlap): # If overlap, remove overlapping region(s) and add merged locus
		merged_intervall_begin = min(begin, min([a for (a,b,v) in overlap]))
		merged_intervall_end = max(end, max([b for (a,b,v) in overlap]))
		merged_intervall_tupple = [item for sublist in [ v for (a,b,v) in overlap] for item in sublist]
		merged_intervall_tupple.append(value) 
		tree.remove_overlap(begin, end)
		tree[merged_intervall_begin:merged_intervall_end] = merged_intervall_tupple
	else:	# Add interval
		tree[begin:end] = [value]
	return 1


# Function to merge loci and return updated dataframe
def merge_loci(df):

	# Cluster the loci
	tree = IntervalTree()
	for i,marker in enumerate(df.index):
		#if marker == "1:109727284":
		#if df.ix[marker,'locus_gene_boundaries'][0] > 109656301-1 and df.ix[marker,'locus_gene_boundaries'][1] < 110060710+1:
		#	print "{}\t{}\t{}\t{}".format(i,df.ix[marker,'locus_gene_boundaries'][0],df.ix[marker,'locus_gene_boundaries'][1],df.ix[marker,'snp_id'])
		#if i ==22:
		insert_interval(tree, df.ix[marker,'locus_gene_boundaries'][0]-1, df.ix[marker,'locus_gene_boundaries'][1]+1,i)

	# Create new dataframe with overlapping loci
	df_addon = pd.DataFrame()
	rows2drop = []
	for i, (start, end, overlapping_loci) in enumerate(tree):
		if len(overlapping_loci) > 1:
			marker = ";".join(df.index[overlapping_loci])
			df_addon.ix[marker,'snp_id'] = ";".join(df.ix[overlapping_loci,'snp_id'])
			df_addon.ix[marker,'locus_start'] = start
			df_addon.ix[marker,'locus_end'] = end
			df_addon.ix[marker,'gwas_pvalue'] = min(df.ix[overlapping_loci,'gwas_pvalue'])
			genes_in_locus_set = set((";".join([ (str(x) if not isinstance(x, float) else '') for x in df.ix[overlapping_loci,'genes_in_locus']])).split(";"))
			genes_in_locus_set.discard('')
			df_addon.ix[marker,'genes_in_locus'] = ";".join(genes_in_locus_set)
			df_addon.ix[marker,'nearest_gene'] = ";".join(df.ix[overlapping_loci,'nearest_gene'])
			df_addon.ix[marker,'chr'] = df.ix[overlapping_loci[0],'chr']
			rows2drop.extend(overlapping_loci)
    
	# Add merged locus and drop ovlerapping loci
	if not df.empty:
		df = df.drop(df.index[rows2drop])
		df = df.append(df_addon)
    	return df


# Function to identify and record HLA SNPs
def record_mhc_snps(df,log_dict,mhc_start,mhc_end):
	mhc_index = df.apply(lambda x : True if ( (x.chr == 6) and ( ( x.locus_end > mhc_start and x.locus_end < mhc_end ) or ( x.locus_start > mhc_start and x.locus_start < mhc_end ) ) ) else False, axis = 1)
	if len(df.index[mhc_index]) > 0:
		log_dict['snps_mhc'] = ';'.join(df.index[mhc_index])
	return mhc_index, log_dict


# Function to identify and record non-autosomal SNPs
def record_sex_chr_snps(df,log_dict):
	non_autosomal_index = df.apply(lambda x : True if x.chr in [str(y) for y in range(1,23,1)] else False, axis = 1)
	if len(df.index[non_autosomal_index]) > 0:
		log_dict['snps_non-autosomal'] = ';'.join(df.index[non_autosomal_index])
	return non_autosomal_index, log_dict


# Function to extend boundaries in case of partly overlapping genes
def get_locus_gene_boundaries(row,df_gene_boundaries):
	locus_min = row.locus_start 
	locus_max = row.locus_end
	locus_genes = row.nearest_gene.split(';')
	locus_genes.extend(row.genes_in_locus.split(';')) if not isinstance(row.genes_in_locus, float)  else None # If it is empty then then it is NaN which is a float
	for gene in set(locus_genes):
		if gene in df_gene_boundaries.index:
			gene_start = df_gene_boundaries.ix[gene,'ensembl_bp_start']
			gene_end = df_gene_boundaries.ix[gene,'ensembl_bp_end']
			locus_min = gene_start if gene_start < locus_min else locus_min
			locus_max = gene_end if gene_end > locus_max else locus_max
		else:
			sys.exit("Exiting: {} not found in gene information file used for gene boundary computations".format(gene)) 
	return locus_min, locus_max


# Helper function to convert column to int
def convert(x):
	try:
		return x.astype(int)
	except:
		return x


# Helper function to run PLINK clumping
def run_plink_clumping(plink_binary, plink_genotype_data_plink_prefix, association_pvalue_cutoff, plink_clumping_distance, plink_clumping_r2, input_filename, output_filename, plink_clumping_snp_column_header, plink_clumping_pvalue_column_header):
	"""
	Function to run PLINK
	"""

	# Make sure the version is correct
	version = subprocess.Popen([plink_binary,"--version"],stdout=subprocess.PIPE).communicate()[0].split()
	month = version[len(version)-2]
	year = int(version[-1].replace(")",""))
	if year < 2015 or ( year == 2015 and month in ["Jan","Feb","Mar","Apr","Jun","Jul"]):
		sys.exit("\nExiting.. Please install a PLINK version 1.9 August relase (or higher)\n")

	#"--exclude", "/tmp/mylist.txt",
	cmd = [plink_binary, 
		"--bfile", plink_genotype_data_plink_prefix,
		"--clump-p1", str(association_pvalue_cutoff),
		"--clump-kb", str(plink_clumping_distance),
		"--clump-r2", str(plink_clumping_r2),
		"--clump-snp-field", plink_clumping_snp_column_header,
		"--clump-field", plink_clumping_pvalue_column_header,
		"--clump", input_filename,
		"--out", output_filename
	]
	return subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()


# Read SNPsnap collection and input SNPs from user
def write_depict_loci_helper(infile,collection,df_gene_boundaries,mhc_start,mhc_end):

	# Dictionary with log information
	log_dict = {}

	# Read users SNPs
	t0 = time()
	try:
		clumps_df = pd.read_csv(infile,delimiter=r"\s+",skip_blank_lines=True)
	except IOError:
		print('\nError when trying to open {}.  Please check the PLINK log.'.format(infile))
		sys.exit()
	clumps_df.index = ["{}:{}".format(a,b) for a,b in zip(clumps_df.CHR.astype('int').tolist(),clumps_df.BP.astype('int').tolist())]
	t1 = time()
	#print 'Time 1 %f' %(t1-t0)

	# Extract user SNPs
	collection_usersnps = collection.loc[collection.index.isin(clumps_df.index)]
	if not collection_usersnps.shape[0]:
		print('\nNo SNPs found after clumping and discarding SNPs not in collection. Please check the PLINK log and your genome build.')
		sys.exit()
	t2 = time()
	#print 'Time 2 %f' %(t2-t1)

	# Identify and remove SNPs missing	
	df_index_remove = pd.DataFrame()
	index_snps_mhc, log_dict = record_mhc_snps(collection_usersnps,log_dict,mhc_start,mhc_end)
	index_snps_non_autosomal, log_dict = record_sex_chr_snps(collection_usersnps,log_dict)
	usersnps_df = collection_usersnps.drop(collection_usersnps.index[ index_snps_mhc | index_snps_non_autosomal])
	t3 = time()
	#print 'Time 3 %f' %(t3-t2)
	
	# Redefine locus boundaries based on gene boundaries
	usersnps_df['locus_gene_boundaries'] = usersnps_df.apply(lambda x : get_locus_gene_boundaries(x, df_gene_boundaries), axis = 1)
	usersnps_df = usersnps_df.join(clumps_df['P'],how="inner")
	usersnps_df.rename(columns={'P': 'gwas_pvalue'}, inplace=True)
	t4 = time()
	#print 'Time 4 %f' %(t4-t3)
	
	# Merge loci
	depictloci_df = pd.DataFrame()
	for chrom in set(usersnps_df.chr):
	    depictloci_df = depictloci_df.append(merge_loci(usersnps_df[usersnps_df.chr == chrom]))
	t5 = time()
	#print 'Time 5 %f' %(t5-t4)

	# Write loci and log to file
	depictloci_df.chr = depictloci_df.chr.apply(convert)
	depictloci_df.locus_start = depictloci_df.locus_start.apply(convert)
	depictloci_df.locus_end = depictloci_df.locus_end.apply(convert)
	depictloci_df.drop(['locus_gene_boundaries','pos'], axis=1, inplace=True)
	t6 = time()
	#print 'Time 6 %f' %(t6-t5)

	return depictloci_df, log_dict


# Construct DEPICT loci
def write_depict_loci(analysis_path,label,association_pvalue_cutoff,collection_file,depict_gene_annotation_file,locus_file,mhc_start,mhc_end,plink_executable,genotype_data_plink_prefix,plink_clumping_distance,plink_clumping_r2,plink_input_file,number_random_runs,background_plink_clumping_pvalue,plink_clumping_snp_column_header,plink_clumping_pvalue_column_header,null_gwas_prefix,depict_contact_email,req_fraction_of_background_files,background_loci_dir_suffix,background_data_path):

	print "\nReading precomputed 1KG SNP collection file"

	# Read SNPsnap collection and gene information
	t0 = time()
	collection = pd.read_csv(collection_file, index_col=0, header=0, delimiter="\t", compression = 'gzip')

	# Read gene information
	df_gene_boundaries = pd.read_csv(depict_gene_annotation_file, header = 0, sep='\t', index_col = 0, usecols = [0,2,3])
	t1 = time()
	#print '%f sec' %(t1-t0)

	print "\nWriting DEPICT loci"
	log_a_out, log_a_err = run_plink_clumping(plink_executable, genotype_data_plink_prefix, association_pvalue_cutoff, plink_clumping_distance, plink_clumping_r2, "{}/{}".format(analysis_path,plink_input_file), "{}/{}".format(analysis_path,label), plink_clumping_snp_column_header, plink_clumping_pvalue_column_header )
	t1 = time()
	depictloci_df, log_dict = write_depict_loci_helper("{}/{}.clumped".format(analysis_path,label),collection,df_gene_boundaries,mhc_start,mhc_end)
	t2 = time()
	#print '%f sec' %(t2-t1)

	# Saving observed loci to file
	depictloci_df.to_csv(locus_file, index=False, quoting=0, doublequote = False, sep='\t',columns=["snp_id","chr","locus_start","locus_end","gwas_pvalue","nearest_gene","genes_in_locus"],float_format="%10.2e")

	print "\nRetrieving background loci"

	# Construct background loci
	def write_background_loci(loci_requested):

		background_loci_dir = "{}/nloci{}_{}".format(background_data_path,loci_requested,background_loci_dir_suffix)
		if not os.path.exists(background_loci_dir):
			os.makedirs(background_loci_dir)
	
			# Draw progress bar
			end_val = number_random_runs
			bar_length=50
	
			# Construct bagground loci sets
			for i in range(0,number_random_runs):
				t3 = time()
				log_back_out, log_back_err = run_plink_clumping(plink_executable,\
					genotype_data_plink_prefix,\
					background_plink_clumping_pvalue,\
					plink_clumping_distance,\
					plink_clumping_r2,\
					"{}_{}.tab.gz".format(null_gwas_prefix,i+1),\
					"{}/{}".format(background_loci_dir,i+1),\
					plink_clumping_snp_column_header,\
					plink_clumping_pvalue_column_header)
				clumped_file = "{}/{}.clumped".format(background_loci_dir,i+1)
				if not os.path.exists(clumped_file):
					sys.stdout.write("PLINK problems with null GWAS {} ignoring\n".format(i))
					continue
				depictloci_background_df, log_back_dict = write_depict_loci_helper(clumped_file,\
					collection, df_gene_boundaries, mhc_start,mhc_end)
				depictloci_background_df.sort('gwas_pvalue',inplace=True)
				background_plink_clumping_pvalue_relaxed = background_plink_clumping_pvalue
				while len(depictloci_background_df) < loci_requested:
					print('Not enough background loci in iteration {} (need: n={}, found: {}).\
						Relaxing clumping p value, consider increasing the parameter "background_plink_clumping_pvalue" \
						in your config file.'.format(loci_requested,len(depictloci_background_df),i))	
					background_plink_clumping_pvalue_relaxed *= 5
					log_back_out, log_back_err = run_plink_clumping(plink_executable,\
						genotype_data_plink_prefix,\
						background_plink_clumping_pvalue_relaxed,\
						plink_clumping_distance,\
						plink_clumping_r2,\
						"{}_{}.tab.gz".format(null_gwas_prefix,i+1),\
						"{}/{}".format(background_loci_dir,i+1),\
						 plink_clumping_snp_column_header,\
						 plink_clumping_pvalue_column_header)
					depictloci_background_df, log_back_dict = write_depict_loci_helper("{}/{}.clumped".format(background_loci_dir,i+1),\
						collection, df_gene_boundaries, mhc_start,mhc_end)
					depictloci_background_df.sort('gwas_pvalue',inplace=True)
				output_file = "{}/permutation{}.txt".format(background_loci_dir,i+1)
				depictloci_background_df.iloc[0:loci_requested,:].to_csv(\
					output_file,\
					index=False,\
					quoting=0,\
					doublequote = False,\
					sep='\t',
					columns=["snp_id","chr","locus_start","locus_end","gwas_pvalue","nearest_gene","genes_in_locus"],float_format="%10.2e")
				os.system("gzip {} ".format(output_file))
				t4 = time()
				#print '\tBackground loci iteration {}: {} sec'.format(i+1,t4-t3)
			        percent = float(i+1) / end_val
			        hashes = '#' * int(round(percent * bar_length))
	       			spaces = ' ' * (bar_length - len(hashes))
			        sys.stdout.write("\rPercent: [{0}] {1}%".format(hashes + spaces, int(round(percent * 100))))
			        sys.stdout.flush()
	
			sys.stdout.write("\n")
		else:
			# Let's make sure that at least 90% of the requested background files are available
			background_files_count = len([name for name in glob.glob('{}/*.gz'.format(background_loci_dir)) if os.path.isfile(name)])
			if background_files_count < (number_random_runs * req_fraction_of_background_files):
				sys.exit("Exiting.. To few background files in {}. Please remove the folder, rerun DEPICT and contact {} if the error prevails.".format(background_loci_dir,depict_contact_email)) 

		return 1 

	write_background_loci(len(depictloci_df))
	#for i in range(450,501):
	#	write_background_loci(i)
	
	return log_a_out + "\n".join(["{}: {}".format(key,log_dict[key]) for key in log_dict]) #, log_a_err 
	

# Helper function to save SNPs that are in my data
def write_plink_input(path, filename, label, marker_col_name, p_col_name, chr_col_name, pos_col_name, sep, genotype_data_plink_prefix, association_pvalue_cutoff):

	print "\nReading GWAS and mapping by chromosome and position to genotype data"

	# Read mapping (Faster than using pandas dataframe, because constructing a chr:pos index takes a hell of time
	with open("%s.bim"%genotype_data_plink_prefix,'r') as infile:
		mapping = {}
		chr_col = 0
		pos_col = 3
		rs_id_col = 1
		for line in infile.readlines():
			words = line.strip().split()
			mapping["%s:%s"%(words[chr_col],words[pos_col])] = words[rs_id_col] 

	# Write PLINK input file
	with ( gzip.open(filename,'r') if '.gz' in filename else open(filename,'r') ) as infile, open("%s/%s_depict.tab"%(path,label),'w') as outfile:
		outfile.write("SNP_chr_pos\tSNP\tChr\tPos\tP\n")

		header = 1
		missing_snps = []
		for line in infile.readlines():
			words = line.strip().split(get_separator(sep))
			if header:
				p_col = words.index(p_col_name) if p_col_name is not None else None
				marker_col = words.index(marker_col_name) if marker_col_name is not None else None
				chr_col = words.index(chr_col_name) if chr_col_name is not None else None
				pos_col = words.index(pos_col_name) if pos_col_name is not None else None
				header = 0
				continue
			if marker_col_name is not None:
				chrom = 23 if words[marker_col].split(":")[0] == "X" else int(words[marker_col].split(":")[0])
				pos = int(words[marker_col].split(":")[1])
				marker_id = words[marker_col] 
			elif chr_col_name is not None and pos_col_name is not None:
				if words[chr_col] == "X":
					chrom = 23
				else:
					chrom = int(words[chr_col])
				pos = int(words[pos_col])
				marker_id = "%s:%s"%(chrom,pos)
			else:
				sys.exit('Please specify either a column for the marker (format <chr:pos>) or the chromosome and position columns.')	

			if marker_id in mapping:
				outfile.write("%s\t%s\t%s\t%s\t%s\n"%(marker_id,mapping[marker_id],chrom,pos,words[p_col]))
			else:
				outfile.write("%s\tNA\t%s\t%s\t%s\n"%(marker_id,chrom,pos,words[p_col]))

				# Only report missing SNPs if below cutoff
				if float(words[p_col]) < association_pvalue_cutoff:
					missing_snps.append(marker_id)

	return {"{} SNPs met your association p value cutoff, but were in the HLA region, on a sex chromosome, or not found in the 1000 Genomes Project data:".format(len(missing_snps)): "{}".format(";".join(missing_snps))} if len(missing_snps)>0 else {}


# Helper function to retrieve PLINK Index SNPs
def get_plink_index_snps(path,label,cutoff,index_snp_col):

	# Read file with correct SNP indices
	id_df = pd.read_csv("%s/%s_depict.tab"%(path,label),index_col=0,header=0,sep="\t")

	# Read PLINK results
	index_snps = []
	with open("%s/%s.clumped"%(path,label),'r') as infile:
		lines = infile.readlines()[1:]
		for line in lines:
			if line.strip() != '':
				words = re.sub(' +','\t',line.strip()).split('\t') 
				index_snps.append(words[index_snp_col]) if words[index_snp_col] != '' else None

	# Re-map to chr:pos
	return id_df.index[id_df.SNP.isin(index_snps)]

# Function to run DEPICT
def run_depict(java_executable, depict_jar, background_data_path, locus_file, label, do_geneprio, do_gsea, do_tissue, ncores, analysis_path, reconstituted_genesets_file, depict_gene_annotation_file, tissue_expression_file, max_top_genes_for_gene_set, nr_repititions, nr_permutations, mhc_start_bp, mhc_end_bp, go_mapping_file, mgi_mapping_file, inweb_mapping_file, tissue_mapping_file, eqtl_mapping_file, eqtl_file, heap_size_in_mb, prioritize_genes_outside_input_loci, leave_out_chr):

	def get_cmd(geneprio_flag, gsea_flag, tissue_flag):
		cmd = [java_executable,	"-Xms512M", "-Xmx{}M".format(heap_size_in_mb),"-XX:+UseParallelGC", '-XX:ParallelGCThreads=3', "-jar",
			depict_jar,			# Below are the arguments to the DEPIT java file
			background_data_path,		# 0  String dataDirectory
	       		locus_file, 			# 1  String filenameLociDefinedBySignificantSNPssAndLDInformation
	       	 	label, 				# 2  String outputFileLabel
       		 	str(int(geneprio_flag)), 	# 3  boolean conductNetworkAnalysis
       			str(int(gsea_flag)),		# 4  boolean conductPathwayAnalysis
       			str(int(tissue_flag)),		# 5  boolean conductTissueAnalysis
       			str(ncores), 			# 6  int nrCores
       			analysis_path,  		# 7  String resultsDirectory
      	 		reconstituted_genesets_file, 	# 8  String cofuncMatrixFile
	       		depict_gene_annotation_file, 	# 9  String filenameGeneAnnotation 
		        tissue_expression_file, 	# 10 String tissueMatrixFile
			str(max_top_genes_for_gene_set),# 11 int maxTopGenesPerGeneSet
			str(nr_repititions),	        # 12 int nrReps
        		str(nr_permutations),		# 13 int nrPerms
        		str(mhc_start_bp),		# 14 int HLAstart
			str(mhc_end_bp), 		# 15 int HLAend
		        go_mapping_file,		# 16 String goMappingFile
		        mgi_mapping_file,		# 17 String mgiMappingFile
		        inweb_mapping_file,		# 18 String inwebMappingFile
		        tissue_mapping_file,		# 19 String tissueMappingFile
			eqtl_mapping_file,		# 20 String filenameGenericIlluminaProbeIDToEnsembl
			eqtl_file,			# 21 String filenameGenericIlluminaProbeIDEQTLs
			str(int(prioritize_genes_outside_input_loci)), 	# 22 boolean calculateGenePrioritizationPValueForGenesOutsideLoci
			leave_out_chr			# 23 String chrToBeLeftOut
		]
		return cmd

	log_out = log_err = None
	if do_gsea and do_tissue:
		print("\nRunning DEPICT{}{}{}".format(" gene prioritization" if do_geneprio else "", " and" if do_geneprio and do_gsea else ""," gene set enrichment analysis" if do_gsea else ""))
		log_a_out,log_a_err = subprocess.Popen(get_cmd(do_geneprio, True, False), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		print("\nRunning DEPICT tissue enrichment analysis".format())
		log_b_out,log_b_err = subprocess.Popen(get_cmd(False, False, True), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		log_out = log_a_out + log_b_out
		log_err = log_a_err + log_b_err
	else:
		print("\nRunning DEPICT{}{}{}".format(" gene prioritization" if do_geneprio else "", " and" if do_geneprio and do_gsea else ""," gene set enrichment analysis" if do_gsea else " tissue enrichment analysis" if do_tissue else ""))
		log_out,log_err = subprocess.Popen(get_cmd(do_geneprio, do_gsea, do_tissue), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	print ""

	return log_out, log_err
