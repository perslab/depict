#!/usr/bin/env python2.7

# Install
#	Python-2.7
#	sudo pip install bx-python
#	pandas
#		Latest versiion (sudo pip install --upgrade pandas)


# pre-req
#	plnk
#	plink binary 1kg files

# TODO
# GIT
# 	new repo
#	version number in git


import os,sys,pdb,logging
from glob import glob
from depict_library import construct_depict_loci,write_plink_input,run_plink


# PLEASE SPECIFY: Paths
analysis_path = "/home/projects/depict/DEPICT-example" # Path to directory where input file lives and all output files are written
depict_path = "/home/projects/depict/DEPICT-example" # Path to DEPICT-example (clone it from GibHub)
depict_code_path = "/home/projects/depict/DEPICT"  # Path to DEPICT jar file (clone it from GitHub)


# PLEASE SPECIFY: Steps that shall be run
step_write_plink_output = False
step_run_plink = False
step_construct_depict_loci = True
step_run_depict = True


# PLEASE SPECIFY: GWAS summary statistics input file parameters (only autosomal SNPs are included)
cutoff =  "5e-8" # "1e-5" 
label = "ldl_teslovich_nature2010" # "P5e-4.N100k.EduYears.STDERR.METAANALYSIS1.TBL.oldSTR" # n=58 SNPs
filename = "%s.txt"%(label) 
pvalue_col = 3 # NB: Counting starts from 0, ie. first columns is referred to as '0'
marker_col = None # Format: <chr:pos>, ie. '6:2321'.  If this column does not exist chr_col and pos_col will be used and this column should be set to 'None'
chr_col = 1 # Does not need to be set if marker_col is set
pos_col = 2 # Does not need to be set if marker_col is set
sep = '\t'


# PLEASE SPECIFY: PLINK and genotype files
plink_binary = "/home/tools/plink/plink_v1-90_stable_beta_3f_2-Mar/plink"
plink_extra_params = "--exclude /home/data/1000G/data/phase1/bed_CEU_GBR_TSI_unrelated/duplicate_markers.rsID"
genotype_data_plink_prefix = "/home/data/1000G/data/phase1/bed_CEU_GBR_TSI_unrelated/CEU_GBR_TSI_unrelated.phase1_release_v3.20101123.snps_indels_svs.genotypes" 


# Locus construction paramenters  (ADVICE: keep default settings)
distance = 1000 # 1000 #distance = 500 
r2 = 0.5 # 0.5 
collection_file = "%s/data/ld0.5_collection_depict_150302_ldl_teslovich_nature2010.txt.gz"%depict_path
locus_file = "%s/%s_loci.txt"%(analysis_path,label)
hla_start = 25000000
hla_stop = 35000000


# DEPICT parameterds (ADVICE: keep default settings)
depict_jar = "%s/dist/Depict.jar"%depict_code_path
ncores = 2
gene_annotation = "GPL570ProbeENSGInfo+HGNC_reformatted.txt"
depict_genelist_file = "GPL570ProbeENSGInfo+HGNC_reformatted.ens"
reconstituted_genesets_file = "%s/data/reconstituted_genesets_example.txt"%depict_path
depict_gene_file = "%s/data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt"%depict_path
depict_gene_information_file = "%s/data/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV65.txt"%depict_path


# Logging (ADVICE: keep default settings)
log_filename = '%s/%s_log.txt'%(analysis_path,label)
logging.basicConfig(filename=log_filename, level=logging.INFO)


# Parse and discards SNPs not in my 1KG data
if step_write_plink_output:
	write_plink_input(analysis_path,filename,label,marker_col,pvalue_col,chr_col,pos_col,sep,genotype_data_plink_prefix)


# Run PLINK to get index SNPs
if step_run_plink:
	plink_log = run_plink(analysis_path, label, genotype_data_plink_prefix, plink_binary, plink_extra_params, cutoff, distance, r2) 
	logging.info(plink_log)


# Read PLINK index SNPs and construct DEPICT locus file
if step_construct_depict_loci:
	loci_log = construct_depict_loci(analysis_path,label,cutoff,collection_file,depict_gene_file,depict_gene_information_file,locus_file,hla_start,hla_stop) #if not os.path.isfile(locus_file) else None
	logging.info(loci_log)


# Gene prioritization and reconstituted gene set enrichment analysis
if step_run_depict:
	depict_log = os.popen("java -Xms512m -Xmx16000m -XX:+UseParallelGC -XX:ParallelGCThreads=3 -jar %s %s %s %s 1 1 0 %s %s %s %s %s "% \
		(depict_jar,"%s/data/"%depict_path,locus_file,label,ncores,analysis_path,reconstituted_genesets_file,gene_annotation,depict_genelist_file)).readlines()
	logging.info(depict_log)
