#!/usr/bin/env python2.7

import warnings
warnings.filterwarnings("ignore")
import os,sys,pdb,logging,subprocess,ConfigParser,argparse
from glob import glob
import pandas as pd

# Read path to config file
parser = argparse.ArgumentParser(description='DEPICT')
parser.add_argument('cfg_file', metavar='DEPICT configuration file', type=str, help='DEPICT configuration file, all user inputs are specified in this file')
args = parser.parse_args()


# Paths
path_to_depict = os.path.realpath(__file__).split('/')
depict_path = "/".join(path_to_depict[0:(len(path_to_depict)-3)])
data_path = "data/"
sys.path.append("%s/src/python"%depict_path)
from depict_library import write_depict_loci,write_plink_input,run_depict


# Parse the config file
cfg = ConfigParser.ConfigParser()
cfg_file_found = cfg.read(args.cfg_file)
analysis_path = cfg.get("PATHS",'analysis_path')
step_write_depict_loci = cfg.getboolean("DEPICT SETTINGS", "step_construct_depict_loci")
step_depict_geneprio = cfg.getboolean("DEPICT SETTINGS", "step_depict_geneprio")
step_depict_gsea = cfg.getboolean("DEPICT SETTINGS", "step_depict_gsea")
step_depict_tissueenrichment = cfg.getboolean("DEPICT SETTINGS", "step_depict_tissueenrichment")


# GWAS summary statistics input file parameters (only autosomal SNPs are included)
association_pvalue_cutoff =  cfg.getfloat("GWAS FILE SETTINGS",'association_pvalue_cutoff')
filename = cfg.get("GWAS FILE SETTINGS",'gwas_summary_statistics_file')
label = cfg.get("GWAS FILE SETTINGS",'label_for_output_files')
pvalue_col_name = cfg.get("GWAS FILE SETTINGS",'pvalue_col_name')
marker_col_name = cfg.get("GWAS FILE SETTINGS",'marker_col_name') if cfg.get("GWAS FILE SETTINGS",'marker_col_name') != '' else None
chr_col_name = cfg.get("GWAS FILE SETTINGS",'chr_col_name') if cfg.get("GWAS FILE SETTINGS",'chr_col_name') != '' else None
pos_col_name = cfg.get("GWAS FILE SETTINGS",'pos_col_name') if cfg.get("GWAS FILE SETTINGS",'pos_col_name') != '' else None
sep = cfg.get("GWAS FILE SETTINGS",'separator')


# PLINK and genotype files
plink_executable = cfg.get("PLINK SETTINGS",'plink_executable')
genotype_data_plink_prefix = "{}/{}".format(depict_path,cfg.get("PLINK SETTINGS",'genotype_data_plink_prefix')) if cfg.get("PLINK SETTINGS",'genotype_data_plink_prefix').startswith("data/genotype_data_plink") else cfg.get("PLINK SETTINGS",'genotype_data_plink_prefix')
plink_clumping_snp_column_header = "SNP"
association_pvalue_cutoff_column_header = "P"
plink_clumping_distance = 500 
plink_clumping_r2 = 0.1

# Locus construction paramenters
collection_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS",'collection_file'))
locus_file = "%s/%s_loci.txt"%(analysis_path,label)
number_random_runs = cfg.getint("MISC SETTINGS",'number_random_runs')
background_plink_clumping_pvalue = cfg.getfloat("MISC SETTINGS",'background_plink_clumping_pvalue')
background_data_path = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","background_data_path"))
null_gwas_prefix = "{}/null_gwas/dgi_1kg_hg19_null".format(background_data_path)
req_fraction_of_background_files = 0.9


# DEPICT parameters
java_executable = "java"
depict_jar = "%s/src/java/dist/Depict.jar"%depict_path
ncores = cfg.getint("MISC SETTINGS",'number_of_threads')
heap_size_in_mb = cfg.getint("MISC SETTINGS",'heap_size_in_mb')
max_top_genes_for_gene_set = cfg.getint("MISC SETTINGS",'max_top_genes_for_gene_set')
nr_repititions = cfg.getint("MISC SETTINGS",'nr_repititions')
nr_permutations = cfg.getint("MISC SETTINGS",'nr_permutations')
mhc_start_bp = cfg.getint("MISC SETTINGS",'mhc_start_bp')
mhc_end_bp = cfg.getint("MISC SETTINGS",'mhc_end_bp')
reconstituted_genesets_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","reconstituted_genesets_file"))
tissue_expression_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","tissue_expression_file"))
depict_gene_annotation_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","depict_gene_annotation_file"))
go_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","go_mapping_file"))
mgi_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","mgi_mapping_file"))
inweb_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","inweb_mapping_file"))
tissue_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","tissue_mapping_file"))
eqtl_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","eqtl_mapping_file"))
eqtl_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","eqtl_file"))
prioritize_genes_outside_input_loci= cfg.getboolean("MISC SETTINGS",'prioritize_genes_outside_input_loci')
leave_out_chr = cfg.get("MISC SETTINGS",'leave_out_chr')


# Construct directory if it does not exist
if not os.path.exists(analysis_path):
	print "\nCreating directory were results will be saved ({})".format(analysis_path)
	os.makedirs(analysis_path)
else:
	print "\nWill store result files to {}".format(analysis_path)


# Logging and contact
log_file = '%s/%s_log.txt'%(analysis_path,label)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO)
depict_contact_email = "tunepers@broadinstitute.org"


# Background loci folder
background_loci_dir_suffix = "nperm{numruns}_kb{dis}_rsq{rsq}_mhc{mhcsta}-{mhcend}_col{collection}/".format(\
	numruns=number_random_runs,\
	dis=plink_clumping_distance,\
	rsq=plink_clumping_r2,\
	mhcsta=mhc_start_bp,\
	mhcend=mhc_end_bp,\
	collection=collection_file.split("/")[-1].replace(".txt.gz","").replace('_','-'))


# Read PLINK index SNPs and write DEPICT locus file
if step_write_depict_loci:
	map_log = write_plink_input(analysis_path, filename, label, marker_col_name, pvalue_col_name, chr_col_name, pos_col_name, sep, genotype_data_plink_prefix, association_pvalue_cutoff)
	logging.info(map_log)

	loci_log = write_depict_loci(analysis_path,label,association_pvalue_cutoff,collection_file,depict_gene_annotation_file,locus_file,mhc_start_bp,mhc_end_bp,plink_executable,genotype_data_plink_prefix,plink_clumping_distance, plink_clumping_r2,"%s_depict.tab"%label,number_random_runs,background_plink_clumping_pvalue,plink_clumping_snp_column_header,association_pvalue_cutoff_column_header,null_gwas_prefix,depict_contact_email,req_fraction_of_background_files,background_loci_dir_suffix,background_data_path)
	logging.info(loci_log)


# Run DEPICT
if step_depict_geneprio or step_depict_gsea or step_depict_tissueenrichment:

	# Check that locus file exists and get the number of loci
	if os.path.isfile(locus_file):
		background_loci_dir = "{}/nloci{}_{}".format(background_data_path,pd.read_csv(locus_file).shape[0],background_loci_dir_suffix)
	else:
		sys.exit("Exiting: Locus file ({}) not found. Please make sure DEPICT loci have been constructed.".format(locus_file)) 

	depict_log = run_depict(java_executable, depict_jar, background_loci_dir, locus_file, label, step_depict_geneprio, step_depict_gsea, step_depict_tissueenrichment, ncores, analysis_path, reconstituted_genesets_file, depict_gene_annotation_file, tissue_expression_file, max_top_genes_for_gene_set, nr_repititions, nr_permutations, mhc_start_bp, mhc_end_bp, go_mapping_file, mgi_mapping_file, inweb_mapping_file, tissue_mapping_file, eqtl_mapping_file, eqtl_file, heap_size_in_mb, prioritize_genes_outside_input_loci, leave_out_chr)
	logging.info(depict_log)
