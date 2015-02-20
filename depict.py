#!/usr/bin/env python2.7
#use Python-2.7

import os
from depict_library import get_loci

depict_jar = "/cvar/jhlab/tp/DEPICT/dist/Depict.jar" # Path to DEPICT jar file (clone it from GitHub)
depict_folder = os.getcwd() # Path to DEPICT-example (clone it from GibHub)
collection_file = "/cvar/jhlab/snpsnap/data/step3/1KG_snpsnap_production_v1_uncompressed_incl_rsID/ld0.5/ld0.5_collection.tab" # Make compressed version
depict_gene_file = "/cvar/jhlab/tp/depict/data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt" # Copy file into Depict/data
depict_gene_information_file = '/cvar/jhlab/tp/depict/data/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV65.txt'

study_label = "ldl_teslovich_nature2010"
marker_file = "ldl_teslovich_nature2010.chrpos"
locus_file = "%s/%s_loci.txt"%(depict_folder,study_label)
results_folder = depict_folder 
ncores = 2
gene_annotation = "GPL570ProbeENSGInfo+HGNC_reformatted.txt"
depict_genelist_file = "GPL570ProbeENSGInfo+HGNC_reformatted.ens"
reconstituted_genesets_file = "%s/data/reconstituted_genesets_example.txt"%depict_folder
hla_start = 25000000
hla_stop = 35000000

# Construct DEPICT locus file
get_loci(collection_file,depict_gene_file,depict_gene_information_file,marker_file,locus_file,hla_start,hla_stop) if not os.path.isfile(locus_file) else None

# Gene prioritization and reconstituted gene set enrichment analysis
#os.system("java -Xms512m -Xmx16000m -XX:+UseParallelGC -XX:ParallelGCThreads=3 -jar %s %s %s %s 1 1 0 %s %s %s %s %s "%(depict_jar,"%s/data/"%depict_folder,locus_file,study_label,ncores,results_folder,reconstituted_genesets_file,gene_annotation,depict_genelist_file))
