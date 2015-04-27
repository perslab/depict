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
genotype_data_plink_prefix = cfg.get("PLINK SETTINGS",'genotype_data_plink_prefix')
plink_extra_params = ""
plink_clumping_distance = 500 
plink_clumping_r2 = 0.1


# Locus construction paramenters
collection_file = "%s/data/collections/%s"%(depict_path,cfg.get("MISC SETTINGS",'collection_file'))
locus_file = "%s/%s_loci.txt"%(analysis_path,label)
hla_start = 25000000
hla_end = 35000000


# DEPICT parameters
java_executable = "java"
depict_jar = "%s/dist/Depict.jar"%depict_path
ncores = cfg.getint("MISC SETTINGS",'number_of_threads')
max_top_genes_for_gene_set = cfg.getint("MISC SETTINGS",'max_top_genes_for_gene_set')
nr_repititions = cfg.getint("MISC SETTINGS",'nr_repititions')
nr_permutations = cfg.getint("MISC SETTINGS",'nr_permutations')
hla_start_bp = cfg.getint("MISC SETTINGS",'hla_start_bp')
hla_end_bp = cfg.getint("MISC SETTINGS",'hla_end_bp')
depict_gene_annotation_file = "%s/data/mapping_and_annotation_files/GPL570ProbeENSGInfo+HGNC_reformatted.txt"%depict_path
depict_genelist_file = "%s/data/mapping_and_annotation_files/GPL570ProbeENSGInfo+HGNC_reformatted.ens"%depict_path
reconstituted_genesets_file = "%s/data/reconstituted_genesets/%s"%(depict_path,cfg.get("MISC SETTINGS","reconstituted_genesets_file"))
depict_gene_file = "%s/data/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt"%depict_path
depict_gene_information_file = "%s/data/mapping_and_annotation_files/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV65.txt"%depict_path
tissue_expression_file = "%s/data/tissue_expression/GPL570EnsemblGeneExpressionPerTissue_DEPICT20130820_z.txt"%depict_path


# Logging
log_file = '%s/%s_log.txt'%(analysis_path,label)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO)


# Parse and discards SNPs not in my 1KG data
if step_write_plink_input:
	map_log = write_plink_input(analysis_path, filename, label, marker_col_name, pvalue_col_name, chr_col_name, pos_col_name, sep, genotype_data_plink_prefix)
	logging.info(map_log)


# Run PLINK to get index SNPs
if step_run_plink:
	plink_log = run_plink_clumping(plink_executable, genotype_data_plink_prefix,  plink_extra_params, cutoff, plink_clumping_distance, plink_clumping_r2, analysis_path, "%s_depict.tab"%label, label)
	logging.info(plink_log)


# Read PLINK index SNPs and construct DEPICT locus file
if step_construct_depict_loci:
	loci_log = construct_depict_loci(analysis_path,label,cutoff,collection_file,depict_gene_file,depict_gene_information_file,locus_file,hla_start,hla_end)
	logging.info(loci_log)


# Run DEPICT
if step_depict_geneprio or step_depict_gsea or step_depict_tissueenrichment:
	loci_log = run_depict(java_executable, depict_jar, "%s/data/"%depict_path, locus_file, label, step_depict_geneprio, step_depict_gsea, step_depict_tissueenrichment, ncores, analysis_path, reconstituted_genesets_file, depict_gene_annotation_file, depict_genelist_file, tissue_expression_file, max_top_genes_for_gene_set, nr_repititions, nr_permutations, hla_start_bp, hla_end_bp)
	logging.info(loci_log)

