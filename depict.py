#!/usr/bin/env python2.7

import os,sys,pdb,logging,subprocess
from glob import glob


# Paths
depict_path = os.getcwd() 
analysis_path = depict_path # Path to directory where input file lives and all output files are written
sys.path.append("%s/src/"%depict_path)
from depict_library import construct_depict_loci,write_plink_input,run_depict
from plink_library import run_plink_clumping


# PLEASE SPECIFY: Steps that shall be run
step_write_plink_output = False
step_run_plink = False
step_construct_depict_loci = False # If you want to run this and the preceeding step please specificy path to PLINK (see below)
step_depict_geneprio = True
step_depict_gsea = True
step_depict_tissueenrichment = False


# PLEASE SPECIFY: GWAS summary statistics input file parameters (only autosomal SNPs are included)
cutoff =  "5e-8"
label = "ldl_teslovich_nature2010"
filename = "%s.txt"%(label) 
pvalue_col = 3 # NB: Counting starts from 0, ie. first columns is referred to as '0'
marker_col = None # Format: <chr:pos>, ie. '6:2321'.  If this column does not exist chr_col and pos_col will be used and this column should be set to 'None'
chr_col = 1 # Set ot 'None' if this column does not exist
pos_col = 2 # Set ot 'None'if this column does not exist
sep = '\t'


# PLEASE SPECIFY: PLINK and genotype files
plink_executable = "/home/tools/plink/plink_v1-90_stable_beta_3f_2-Mar/plink" # Change to your PLINK executable
genotype_data_plink_prefix = "%s/data/genotype_data_plink/ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes"%depict_path # Point this to the PLINK formated genotype data (point to filename but leave out the filename extension)
plink_extra_params = ""
plink_clumping_distance = 500 # (ADVICE: keep default settings for the below entries)
plink_clumping_r2 = 0.1
plink_clumping_snp_column = "SNP" 
plink_clumping_pvalue_column = "P"


# Locus construction paramenters
collection_file = "%s/data/collections/ld0.5_collection_depict_150315.txt.gz"%depict_path
locus_file = "%s/%s_loci.txt"%(analysis_path,label)
hla_start = 25000000
hla_end = 35000000


# DEPICT parameterds
java_executable = "java"
depict_jar = "%s/dist/Depict.jar"%depict_path
ncores = 2
gene_annotation = "GPL570ProbeENSGInfo+HGNC_reformatted.txt"
depict_genelist_file = "GPL570ProbeENSGInfo+HGNC_reformatted.ens"
reconstituted_genesets_filename = "GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary"
reconstituted_genesets_file = "%s/data/reconstituted_genesets/%s"%(depict_path,reconstituted_genesets_filename) 
depict_gene_file = "%s/data/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt"%depict_path
depict_gene_information_file = "%s/data/mapping_and_annotation_files/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV65.txt"%depict_path


# Logging
log_filename = '%s/%s_log.txt'%(analysis_path,label)
logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)


# Parse and discards SNPs not in my 1KG data
if step_write_plink_output:
	write_plink_input(analysis_path,filename,label,marker_col,pvalue_col,chr_col,pos_col,sep,genotype_data_plink_prefix)


# Run PLINK to get index SNPs
if step_run_plink:
	plink_log = run_plink_clumping(plink_executable, genotype_data_plink_prefix,  plink_extra_params, cutoff, plink_clumping_distance, plink_clumping_r2, plink_clumping_snp_column, plink_clumping_pvalue_column, analysis_path, "%s.tab"%label, label)
	logging.info(plink_log)


# Read PLINK index SNPs and construct DEPICT locus file
if step_construct_depict_loci:
	loci_log = construct_depict_loci(analysis_path,label,cutoff,collection_file,depict_gene_file,depict_gene_information_file,locus_file,hla_start,hla_end)
	logging.info(loci_log)


# Gene prioritization and reconstituted gene set enrichment analysis
if step_depict_geneprio or step_depict_gsea or step_depict_tissueenrichment:
	loci_log = run_depict(java_executable, depict_jar, "%s/data/"%depict_path, locus_file, label, step_depict_geneprio, step_depict_gsea, step_depict_tissueenrichment, "%s"%ncores, analysis_path, reconstituted_genesets_file, gene_annotation, depict_genelist_file )
	logging.info(loci_log)

