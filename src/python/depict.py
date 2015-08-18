#!/usr/bin/env python2.7

import os,sys,pdb,logging,subprocess,ConfigParser,argparse
from glob import glob

# Read path to config file
parser = argparse.ArgumentParser(description='DEPICT')
parser.add_argument('cfg_file', metavar='DEPICT configuration file', type=str, help='DEPICT configuration file, all user inputs are specified in this file')
args = parser.parse_args()


# Paths
path_to_depict = os.path.realpath(__file__).split('/')
depict_path = "/".join(path_to_depict[0:(len(path_to_depict)-2)])
data_path = "data/"
sys.path.append("%s/src/"%depict_path)
from depict_library import construct_depict_loci,write_plink_input,run_depict
from plink_library import run_plink_clumping


# Parse the config file
cfg = ConfigParser.ConfigParser()
cfg_file_found = cfg.read(args.cfg_file)
analysis_path = cfg.get("PATHS",'analysis_path')
step_write_plink_input = cfg.getboolean("DEPICT SETTINGS", "step_write_plink_input")
step_run_plink = cfg.getboolean("DEPICT SETTINGS", "step_run_plink")
step_construct_depict_loci = cfg.getboolean("DEPICT SETTINGS", "step_construct_depict_loci")
step_depict_geneprio = cfg.getboolean("DEPICT SETTINGS", "step_depict_geneprio")
step_depict_gsea = cfg.getboolean("DEPICT SETTINGS", "step_depict_gsea")
step_depict_tissueenrichment = cfg.getboolean("DEPICT SETTINGS", "step_depict_tissueenrichment")


# GWAS summary statistics input file parameters (only autosomal SNPs are included)
cutoff =  cfg.getfloat("GWAS FILE SETTINGS",'association_pvalue_cutoff')
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
plink_clumping_distance = 500 
plink_clumping_r2 = 0.1


# Locus construction paramenters
collection_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS",'collection_file'))
locus_file = "%s/%s_loci.txt"%(analysis_path,label)
hla_start = 25000000
hla_end = 35000000


# DEPICT parameters
java_executable = "java"
depict_jar = "%s/dist/Depict.jar"%depict_path
ncores = cfg.getint("MISC SETTINGS",'number_of_threads')
heap_size_in_mb = cfg.getint("MISC SETTINGS",'heap_size_in_mb')
max_top_genes_for_gene_set = cfg.getint("MISC SETTINGS",'max_top_genes_for_gene_set')
nr_repititions = cfg.getint("MISC SETTINGS",'nr_repititions')
nr_permutations = cfg.getint("MISC SETTINGS",'nr_permutations')
hla_start_bp = cfg.getint("MISC SETTINGS",'hla_start_bp')
hla_end_bp = cfg.getint("MISC SETTINGS",'hla_end_bp')
reconstituted_genesets_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","reconstituted_genesets_file"))
tissue_expression_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","tissue_expression_file"))
depict_gene_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","depict_gene_file"))
depict_gene_annotation_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","depict_gene_annotation_file"))
depict_gene_information_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","depict_gene_information_file"))
go_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","go_mapping_file"))
mgi_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","mgi_mapping_file"))
inweb_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","inweb_mapping_file"))
tissue_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","tissue_mapping_file"))
eqtl_mapping_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","eqtl_mapping_file"))
eqtl_file = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","eqtl_file"))
background_data_path = "{}/{}".format(depict_path,cfg.get("MISC SETTINGS","background_data_path"))
prioritize_genes_outside_input_loci= cfg.getboolean("MISC SETTINGS",'prioritize_genes_outside_input_loci')
leave_out_chr = cfg.get("MISC SETTINGS",'leave_out_chr')

# Logging
log_file = '%s/%s_log.txt'%(analysis_path,label)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO)


# Parse and discards SNPs not in my 1KG data
if step_write_plink_input:
	map_log = write_plink_input(analysis_path, filename, label, marker_col_name, pvalue_col_name, chr_col_name, pos_col_name, sep, genotype_data_plink_prefix, cutoff)
	logging.info(map_log)


# Run PLINK to get index SNPs
if step_run_plink:
	plink_log = run_plink_clumping(plink_executable, genotype_data_plink_prefix, cutoff, plink_clumping_distance, plink_clumping_r2, analysis_path, "%s_depict.tab"%label, label)
	logging.info(plink_log)


# Read PLINK index SNPs and construct DEPICT locus file
if step_construct_depict_loci:
	loci_log = construct_depict_loci(analysis_path,label,cutoff,collection_file,depict_gene_information_file,locus_file,hla_start,hla_end)
	logging.info(loci_log)


# Run DEPICT
if step_depict_geneprio or step_depict_gsea or step_depict_tissueenrichment:
	loci_log = run_depict(java_executable, depict_jar, background_data_path, locus_file, label, step_depict_geneprio, step_depict_gsea, step_depict_tissueenrichment, ncores, analysis_path, reconstituted_genesets_file, depict_gene_annotation_file, depict_gene_file, tissue_expression_file, max_top_genes_for_gene_set, nr_repititions, nr_permutations, hla_start_bp, hla_end_bp, go_mapping_file, mgi_mapping_file, inweb_mapping_file, tissue_mapping_file, eqtl_mapping_file, eqtl_file, heap_size_in_mb, prioritize_genes_outside_input_loci, leave_out_chr)
	logging.info(loci_log)

