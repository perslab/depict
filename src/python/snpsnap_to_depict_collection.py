#!/usr/bin/python

import pdb,math
import pandas as pd
from bx.intervals.cluster import ClusterTree
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalTree
from datetime import date


depict_data_path = "/home/projects/depict/git/DEPICT/data/"
depictgenefile = "%s/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt"%depict_data_path
depict_gene_information_file = "%s/mapping_and_annotation_files/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV65.txt"%depict_data_path


# SNPsnap collections for 1000 Genomes Project phase 3 can be downloaded from http://www.broadinstitute.org/mpg/snpsnap/database_download.html
geno_flag =  "1000genomespilot"
#geno_flag = "1000genomesphase3" 
collection_id = "ld0.5_collection"
collection_file = "/home/data/snpsnap/data/{}/{}.tab.gz".format(geno_flag,collection_id)


# Function construct tree with gene intervals
def get_nearest_gene_intervall_tree(depict_gene_information_file, depictgenes):
    trees = {}
    for i in range(1, 23, 1):
            trees[str(i)] = IntervalTree()
    with open (depict_gene_information_file,'r') as infile:
        for line in infile.readlines()[1:]:
            words = line.strip().split('\t')
            if words[0] in depictgenes and words[4] in [str(x) for x in range(1,23,1)]:
                tss = int(words[5]) if words[7] == '1' else int(words[6])
                trees[words[4]].insert_interval(Interval(tss, tss, value=words[0])) if words[0] in depictgenes and words[4] in [str(x) for x in range(1,23,1)] else None
    return trees


# Function to get nearest gene
def get_nearest(chrom, pos):
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


# Function to filter out genes not in DEPICT
def limit_to_depict_genes(row):
    genes = []
    if not isinstance(row['ID_genes_in_matched_locus'], float):
        for gene in row['ID_genes_in_matched_locus'].split(';'):
            genes.append(gene) if gene in depictgenes else None
    return ";".join(set(genes))


# Read SNPsnap collection
collection = pd.read_csv(collection_file, index_col=0, header=0, delimiter="\t",compression="gzip")


# Add chromosome and position columns
collection['chr'] = pd.Series([x.split(":")[0] for x in collection.index], index=collection.index)
collection['pos'] = pd.Series([x.split(":")[1] for x in collection.index], index=collection.index)


# Exclude non-autosomal SNPs
non_autosomal_index = collection.apply(lambda x : True if x.chr in [str(y) for y in range(23,25,1)] else False, axis = 1)
collection.drop(collection.index[non_autosomal_index], inplace=True)


# Read genes covered by DEPICT
with open (depictgenefile,'r') as infile:
    depictgenes = [line.strip() for line in infile.readlines()]


# Construct gene interval tree and identify nearest genes for all SNPs in the collection
trees = get_nearest_gene_intervall_tree(depict_gene_information_file,depictgenes)
chr_1_bps = 248956422
my_max_dist = chr_1_bps / 10
collection['nearest_gene'] =  collection.apply(lambda x : get_nearest(x.chr, int(x.pos)),axis=1)


# Limit to genes in DEPICT
collection['genes_in_locus'] = collection.apply(limit_to_depict_genes,axis=1) # 2h 48min per loop"


# Rename columns and save selected columns
collection.rename(columns={'rsID': 'snp_id', 'loci_upstream': 'locus_start','loci_downstream': 'locus_end'}, inplace=True)
collection.rename(columns={'nearest_genes': 'nearest_gene'}, inplace=True)
collection.to_csv("%s/collections/%s_%s_depict_%s.txt"%(depict_data_path,collection_id,geno_flag,date.today().strftime("%y%m%d")), index=True, quoting=0, doublequote = False, sep='\t',columns=["snp_id","chr","pos","locus_start","locus_end","nearest_gene","genes_in_locus"])
