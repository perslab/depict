#!/usr/bin/env python2.7

import os,sys,pdb,logging
from glob import glob


# PLEASE SPECIFY: Paths
depict_path = os.getcwd() # Path to DEPICT-example (clone it from GibHub)
analysis_path = depict_path # Path to directory where input file lives and all output files are written
sys.path.append("%s/src/"%depict_path)
from depict_library import construct_depict_loci,write_plink_input,run_plink


# PLEASE SPECIFY: Steps that shall be run
step_write_plink_output = True
step_run_plink = True
step_construct_depict_loci = True # If you want to run this and the preceeding step please specificy path to PLINK (see below)
step_run_depict = True


# PLEASE SPECIFY: GWAS summary statistics input file parameters (only autosomal SNPs are included)
cutoff =  "5e-8"
label = "ldl_teslovich_nature2010"
filename = "%s.txt"%(label) 
pvalue_col = 3 # NB: Counting starts from 0, ie. first columns is referred to as '0'
marker_col = None # Format: <chr:pos>, ie. '6:2321'.  If this column does not exist chr_col and pos_col will be used and this column should be set to 'None'
chr_col = 1 # Does not need to be set if marker_col is set
pos_col = 2 # Does not need to be set if marker_col is set
sep = '\t'


# PLEASE SPECIFY: PLINK and genotype files
plink_binary = "/home/tools/plink/plink_v1-90_stable_beta_3f_2-Mar/plink"
plink_extra_params = ""
genotype_data_plink_prefix = "%s/data/CEU_GBR_TSI_unrelated.phase1_release_v3.20101123.snps_indels_svs.genotypes_ldl_teslovich_nature2010"%depict_path


# Locus construction paramenters  (ADVICE: keep default settings)
distance = 1000 
r2 = 0.5 
collection_file = "%s/data/ld0.5_collection_depict_150302.txt.gz"%depict_path
locus_file = "%s/%s_loci.txt"%(analysis_path,label)
hla_start = 25000000
hla_stop = 35000000


# DEPICT parameterds (ADVICE: keep default settings)
depict_jar = "%s/dist/Depict.jar"%depict_path
ncores = 2
gene_annotation = "GPL570ProbeENSGInfo+HGNC_reformatted.txt"
depict_genelist_file = "GPL570ProbeENSGInfo+HGNC_reformatted.ens"
reconstituted_genesets_file = "%s/data/reconstituted_genesets_example.txt"%depict_path
depict_gene_file = "%s/data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt"%depict_path
depict_gene_information_file = "%s/data/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV65.txt"%depict_path


# Logging (ADVICE: keep default settings)
log_filename = '%s/%s_log.txt'%(analysis_path,label)
logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)


# Parse and discards SNPs not in my 1KG data
if step_write_plink_output:
	write_plink_input(analysis_path,filename,label,marker_col,pvalue_col,chr_col,pos_col,sep,genotype_data_plink_prefix)


# Run PLINK to get index SNPs
if step_run_plink:
	plink_log = run_plink(analysis_path, label, genotype_data_plink_prefix, plink_binary, plink_extra_params, cutoff, distance, r2) 
	logging.info(plink_log)


# Read PLINK index SNPs and construct DEPICT locus file
if step_construct_depict_loci:
	loci_log = construct_depict_loci(analysis_path,label,cutoff,collection_file,depict_gene_file,depict_gene_information_file,locus_file,hla_start,hla_stop)
	logging.info(loci_log)


# Gene prioritization and reconstituted gene set enrichment analysis
if step_run_depict:
	# Arguments to Java binary
        #  0  String dataDirectory
        #  1  String filenameLociDefinedBySignificantSNPssAndLDInformation
        #  2  String outputFileLabel
        #  3  boolean conductNetworkAnalysis (Integer coding)
        #  4  boolean conductPathwayAnalysis
        #  5  boolean conductTissueAnalysis
        #  6  int nrCores
        #  7  String resultsDirectory
        #  8  String cofuncMatrixPath
        #  9  String filenameGeneAnnotation = dataDirectory + "/" + args[9]
        #  10 String confineAnalysisToSubsetOfEnsemblGenes = dataDirectory + "/" + args[10]
	depict_log = os.popen("java -Xms512m -Xmx16000m -XX:+UseParallelGC -XX:ParallelGCThreads=3 -jar %s %s %s %s 1 1 0 %s %s %s %s %s "% \
		(depict_jar,"%s/data/"%depict_path,locus_file,label,ncores,analysis_path,reconstituted_genesets_file,gene_annotation,depict_genelist_file)).readlines()
	logging.info(depict_log)
