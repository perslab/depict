#!/usr/bin/env python2.7

import os
import sys
import time
import argparse
import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation ### installation: pip install scikit-learn

import subprocess

import math
import re
import glob

import pdb

import ConfigParser

###################################### EXTRA LOADING ######################################

### To turn of "SettingWithCopyWarning" warnings:
#pd.options.mode.chained_assignment = None  # default='warn'

###################################### DEPENDENCIES ######################################

# 1) Cytoscape 3.2.1 or greater. Cytoscape 3.1.1 does not work.
# 2) Java 1.7, 1.8 or greater
# 3) Python packages
	# numpy
	# pandas
	# scikit-learn (pip install scikit-learn)

### Check your Java version ### 
# OSX: /usr/libexec/java_home -V
# Linux/UNIX: java -version

### Download Java 1.8
# http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html

###################################### USAGE ######################################

### Example Usage
#python network_plot.py <config_file.cfg>

###################################### Known bugs ######################################
# 2015-09-19: "Missing edges" with Pearson_correlation close to the cut-off
	# For some reason, Cytoscape does not draw edges between nodes with Pearson_correlation close to the cut-off.
	# This is a harmless bug in Cytoscape, but might confuse you if you try to recontruct the graph from the network table.




#######################################################################################
###################################### CONSTANTS ######################################
#######################################################################################

#NOTE: These are internal variables to the program and should only be updated if the DEPICT output is updated from DEPICT version 1.0 rel155
COLS2READ_GENESETENRICHMENT = ["Original gene set ID", "Original gene set description", "Nominal P value", "False discovery rate"]

RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL = "-" # update this symbol if the header symbol of the gene names changes

EXPECTED_DEPICT_FDR_CUTOFFS = ["< 0.01","<= 0.01","< 0.05","<= 0.05","< 0.20","<= 0.20",">= 0.20","> 0.20"] # expected DEPICT cut-offs
EXPECTED_DEPICT_FDR_CUTOFFS += ["<0.01","<=0.01","<0.05","<=0.05","<0.20","<=0.20",">=0.20",">0.20"] # PT added 2018-08-25 for compatibility with RPgwas

EXPECTED_CUTOFF_TYPES = ["fdr", "pvalue"] # expected cutoff_type


#######################################################################################
###################################### FUNCTIONS ###################################### 
#######################################################################################

def subset_genesetenrichment_data_frame(df, cutoff_type):
	"""
	Subsets the df_genesetenrichment.
	"""
	# OBS: list_of_fdr_cutoffs is a global variable. Not so nice, but should be fine...
	
	df_genesetenrichment_raw = df # *IS DEEP COPY NEEDED?* --> only if you change columns value in the data frame
		# copy (legacy reasons - I am too lazy to change variable names...)
	

	if cutoff_type == "fdr":
		print "cutoff_type={} | Using the following FDR values as significant gene sets: [{}]".format( cutoff_type, ",".join(list_of_fdr_cutoffs) )

		### SUBSET data based on FDR
		df_genesetenrichment_subset = df_genesetenrichment_raw[df_genesetenrichment_raw['False discovery rate'].isin(list_of_fdr_cutoffs)]
			# ^^ REF: http://stackoverflow.com/questions/12096252/use-a-list-of-values-to-select-rows-from-a-pandas-dataframe
	elif cutoff_type == "pvalue":
		print "cutoff_type={} | Using the following P-value as significant gene sets: [p-value <= {}]".format( cutoff_type, pvalue_cutoff )
		
		df_genesetenrichment_subset = df_genesetenrichment_raw[df_genesetenrichment_raw['Nominal P value'] <= pvalue_cutoff]

	else:
		raise Exception("INTERNAL ERROR in source code: 'if statements' related to determining cutoff_type is not updated/synchronized with EXPECTED_CUTOFF_TYPES.")

	
	n_gene_set_significant = df_genesetenrichment_subset.shape[0]
	if n_gene_set_significant == 0:
		print "Found no gene sets with the specified cutoff_type and cutoff value in file_genesetenrichment. The program cannot continue the analysis."
		print "To solve this issue, you may try to run network_plot.py using P-values as 'cutoff_type' with a non-stringent cut-off value (see the config file for further description)."
		print "Will exit the programme."
		sys.exit(0)
	else:	
		print "Extracted {} number of significant gene sets. These gene sets will be used for further analysis".format(n_gene_set_significant)

	return df_genesetenrichment_subset



def add_cluster_results_to_data_frame(df):
	""" 
	Function will add columns with clustering results to data frame
	The data frame must contain the column 'Original gene set ID' (data frame for the DEPICT enrichments)
	"""

	df.set_index('Original gene set ID', drop=False, inplace=True) # do not delete columns after.
	df.index.name = 'idx' # just to be nice... not important

	### Initialyze columns with APPROPRIATE default values
	df['Cluster ID'] = np.nan # will become integer later
	df['Cluster center (boolean)'] = False
	df['Cluster minimum P value (boolean)'] = False
	df['Cluster gene set with minimum P value'] = np.nan
	df['Cluster minimum P value'] = np.nan


	################## Assigning cluster labels to data frame ##################
	### using Pandas match
	#http://pandas.pydata.org/pandas-docs/stable/comparison_with_r.html#match
	#http://stackoverflow.com/a/15866201
	#pd.Series(pd.match(s,[2,4],np.nan)) # <-- *CONSIDER USING THIS*

	### Using index
	# *by indexing df according to "df_reconstituted_genesets.columns" we order to data frame correctly when assignning the labels*
	# IMPORTANT --> REQUIRES index to be 'Original gene set ID'
	#df.ix[['GO:0003712', 'GO:0003713'],:]
	df.ix[df_reconstituted_genesets.columns,'Cluster ID'] = labels
		# ^*OBS: these labels will be reassigned later | New September 8th 2015

	################## Assigning Cluster center (boolean) and Cluster minimum P value ##################

	for k in range(n_clusters): # looping over clusters
		### REMEMBER ### 
		# We are mapping *FROM* column index/names in df_reconstituted_genesets, that was used for similarity calculations, *TO* df/df_genesetenrichment to get descriptions (and p-values) of Gene Sets.
		# We cannot just used the cluster_center_indices as index, because these are *NOT* in sync with the df/df_genesetenrichment
		
		### saving cluster center (boolean)s
		cluster_center_ID = df_reconstituted_genesets.columns[cluster_centers_indices[k]] # cluster_center_ID --> a string
		df.ix[cluster_center_ID, 'Cluster center (boolean)'] = True
		
		class_members_ID = df_reconstituted_genesets.columns[labels == k] # array of strings/genesetsIDs. selecting gene sets names that are members of the current cluster
		
		### saving class member with the lowest p-value (including the cluster center (boolean) observation)
		idx_min_pval = df.ix[class_members_ID, "Nominal P value"].idxmin()
		df.ix[idx_min_pval, 'Cluster minimum P value (boolean)'] = True
		

		### saving more...
		df.ix[class_members_ID, 'Cluster minimum P value'] = df.ix[idx_min_pval, "Nominal P value"] # scalar
		df.ix[class_members_ID, 'Cluster gene set with minimum P value'] = df.ix[idx_min_pval, "Original gene set ID"] # scalar


	## safety check
	assert( sum(df["Cluster center (boolean)"])==n_clusters )
	assert( sum(df["Cluster minimum P value (boolean)"])==n_clusters )

	################## [SEMI HACK] *REASSIGNING CLUSTER LABELS/IDs*  ##################
	# Here we assign the lowest cluster ID/label to the cluster with the smallest "Cluster minimum P value"
	tmp_current_clusterID = None
	tmp_clusterID_rank = 0
	df = df.sort_values(['Cluster minimum P value', 'Cluster ID']) # SORT!
	for row_index, row in df.iterrows():
		# NB: we could just loop over the index instead!
		if tmp_current_clusterID is None:
			tmp_current_clusterID = row['Cluster ID']

		if tmp_current_clusterID != row['Cluster ID']:
			tmp_clusterID_rank += 1
			tmp_current_clusterID = row['Cluster ID']

		df.ix[row_index, 'Cluster ID'] = tmp_clusterID_rank


	### OLD CODE
	# tmp_cluster_rank = 0
	# for k in range(n_clusters): # looping over clusters | we could also loop over unique values in df['Cluster ID']
	# 	df[df['Cluster ID']==k]

		# http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.rank.html
		# Series.rank(method='average', na_option='keep', ascending=True, pct=False)
		# rank supports different tie-breaking methods, specified with the method parameter:
			# average : average rank of tied group
			# min : lowest rank in the group
			# max : highest rank in the group
			# first : ranks assigned in the order they appear in the array


	return df


def set_cytoscape_node_label_text(string):
	"""
	Function to set a PROPER name for the cytoscape node labels
	
	ISSUE. The main issue here is to deal with is: KEGG, REACTOME and "PPI" (PPI subnetwork) gene sets. Examples:
		1) REACTOME: 	REACTOME_DEPOLARIZATION_OF_THE_PRESYNAPTIC_TERMINAL_TRIGGERS_THE_OPENING_OF_CALCIUM_CHANNELS --> Depolarization of the presynaptic terminal triggers the opening of calcium channels
		2) KEGG: 		KEGG_ERBB_SIGNALING_PATHWAY --> erbb signaling pathway
		3) PPI:			KLF1 PPI subnetwork --> KLF1 protein complex
	
	SOLUTION.
		KEGG/REACTOME:
			a) strip underscores
			b) strip leading KEGG/REACTOME
			c) convert string to sentence case (all lowercase, except first letter)
		PPI:
			a)


	ISSUES THAT ARE UNFIXABLE (hard to fix):
	1) Capitalized gene names in KEGG/REACTOME will be *lost*. E.g. the label will be "erbb signaling pathway" instead of "ERBB signaling pathway"
	"""
	pattern_REACTOME = re.compile(r"^REACTOME_", re.IGNORECASE)
	pattern_KEGG = re.compile(r"^KEGG_", re.IGNORECASE)
	pattern_PPI = re.compile(r"PPI subnetwork$", re.IGNORECASE)

	if pattern_REACTOME.search(string):
		#print "REACTOME"
		tmp = pattern_REACTOME.sub("", string) # strip pattern
		tmp = tmp.replace ("_", " ") # replace underscore with whitespace
		tmp = tmp.capitalize() # return a copy of word with only its first character capitalized. [str.lower()/str.upper()]
	elif pattern_KEGG.search(string):
		#print "KEGG"
		tmp = pattern_KEGG.sub("", string) # strip pattern
		tmp = tmp.replace ("_", " ") # replace underscore with whitespace
		tmp = tmp.capitalize() # return a copy of word with only its first character capitalized. [str.lower()/str.upper()]
	elif pattern_PPI.search(string):
		#print "PPI"
		#tmp = pattern_PPI.sub("protein complex", string) # strip pattern
		tmp = pattern_PPI.sub("PPI", string) # strip pattern | *2015-12 TUNE SUGGESTION*
	else: # WE HAVE A NICELY FORMATTED string - DO NOT DO ANYTHING
		tmp = string

	node_label_text = tmp # hmm, a bit redundant...

	return node_label_text

def set_cytoscape_node_height_and_width(string):
    """ Function to compute node height and width """
    string_length = len(string)
    
    CHAR_WIDTH = 6
    NODE_MAX_WIDTH = 120
    DEFAULT_HEIGHT = 30
    ADD_HEIGHT_PER_LINE = 20
    BUFFER = 20

    # Width
    width = min(string_length * CHAR_WIDTH, NODE_MAX_WIDTH) 
    node_width = width + BUFFER
    label_width = node_width

    # Height
    height_factor = int(max(1, math.ceil((string_length*CHAR_WIDTH)/float(NODE_MAX_WIDTH))))  # ceil returns float, therefor we use int(). ceil: "rounds up to nearest integer"
    node_height = DEFAULT_HEIGHT + (height_factor - 1) * ADD_HEIGHT_PER_LINE

    return (node_width, node_height, label_width)

def write_network_table_file(df,file_out):
	""" 
	Function to write cytoscape network file (.network_table)
	*IMPORTANT* the INDEX of the input df MUST contain the index of the genesets to draw edges between.
	
	TODO IMPROVEMENT: 
		you could also re-calculate the correlation matrix for the input data frame. 
		You would then subset the "df_reconstituted_genesets" based in the index of the input df.

	"""
	print "Writing network table file..."
	print "Got data frame with {} rows".format(len(df))


	f_network_table = open(file_out, 'w')
	#header_network_table = ["Source","Target","S_name","T_name","Pearson_correlation","Pearson_correlation_discrete"]
	header_network_table = ["Source","Target","Pearson_correlation","Pearson_correlation_discrete"]
	f_network_table.write( "\t".join(header_network_table)+"\n" )

	### Looping over *UPPER* triangle of the SYMMETRIC matrix/dataframe 
	n_rows, n_cols = matrix_corr.shape[0], matrix_corr.shape[1]
	assert(n_rows == n_cols) # SYMMETRIC
	n = n_rows # abitrary

	### ***** FIX THIS **** : you need to loop over selected variables in data frame ***####

	for i_row in range(n-1):
		flag_singleton = True # "raise" this flag for each geneset.
		gsID_source = matrix_corr.index[i_row] # same as matrix_corr.columns[i_row] <-- because it is symmetric
		if not df.index.isin([gsID_source]).any(): continue # *OBS* important

		for j_col in range(i_row+1,n):
			gsID_target = matrix_corr.columns[j_col] # same as matrix_corr.index[j_col] <-- because it is symmetric
			if not df.index.isin([gsID_target]).any(): continue # *OBS* important
				### ---> both source and target MUST be in the input data frame.

			corr = matrix_corr.ix[i_row, j_col] # look-up correlation in the correlation matrix
			if corr >= network_correlation_cutoff:
				flag_singleton = False # withdraw the flag: we found a correlation partner. 
				line = "{}\t{}\t{}\t{}\n".format(gsID_source, gsID_target, corr, int(round(corr,1)*10) )
				f_network_table.write(line)
		if flag_singleton:
			line = "{}\t{}\t{}\t{}\n".format(gsID_source, "", "", "" )
			f_network_table.write(line)

	f_network_table.close()
	print "Done writing network table file"

	################## Other ways to working with networks in Python ##################
	### NetworkX - good
	#http://stackoverflow.com/questions/21207872/construct-networkx-graph-from-pandas-dataframe
	#graph = nx.from_numpy_matrix(df2.values)

	### iGraph
	# http://igraph.org/python/doc/tutorial/tutorial.html

	### pynetconv
	# http://pynetconv.sourceforge.net/api/
	
	### THINGs to try ###
	### 1)
	#SEE http://pandas.pydata.org/pandas-docs/stable/advanced.html#creating-a-multiindex-hierarchical-index-object
	# iterables = [df.index, df.index]
	# idx = pd.MultiIndex.from_product(iterables, names=['first', 'second'])
	### 2) 
	#TRY --> numpy.triu (triangle-upper)
	### 3) 
	# matrix_corr[matrix_corr > 0.3]


def write_genesetenrichment_cluster_result_file(df):
	### Indexing Panda objects can return two fundamentally different objects: a view or a copy.
	### REF: http://stackoverflow.com/questions/17995328/changing-values-in-pandas-dataframe-doenst-work
	### OBS: any *modification* to the data frame in "view", modifies the "original" data frame
	df_write = df.copy(deep=True) # deep "copy" - *OBS: INEFFICIENT MEMORY USAGE*

	### Columns
	# Original gene set ID
	# Original gene set description
	# Nominal P value
	# False discovery rate
	# Cluster ID
	# Cluster center (boolean)
	# Cluster minimum P value (boolean)
	# Cluster gene set with minimum P value
	# Cluster minimum P value
	print "Writing geneset enrichment cluster result file..."
	df_write["Cluster ID"] = df_write["Cluster ID"] + 1 # looks better in the output
	df_write.sort_values(['Cluster ID','Nominal P value'], inplace=True)
	df_write.to_csv(file_genesetenrichment_cluster_result, sep="\t", index=False)
	print "Done"


########### COPY FROM IPython NOTEBOOK - *NOT TESTED* ###########
# def discretize_pval_minuslogten(pval, scale, return_value_type):
#     # pval is *FLOORED* "minus_log10" values. 
#         # EXAMPLES
#         # pval=0 --> real_pval=10^-0=1
#         # pval=1 --> real_pval=10^-1=0.1
#         # pval=5 --> real_pval=10^-5=0.00001
#     # *Remember*: larger "pval" values corrresponds to smaller (MORE significant) p-values
# 	# "scale" is a numpy array generated by linspace - increasing order.
# 	#*OBS*: we do not use the first [lowest] element of "scale"
# 	boundary_lowest_possible_value = 0 # this is the lowest possible "minus_log10" pvalue. Remember that 10^-0=1
# 	boundary_low = scale[1]
# 	boundary_medium = scale[2]
# 	boundary_high = scale[3]

# 	value_discretized_numeric = -99 # TODO: find better initial value!? It MUST be numeric.
# 	value_discretized_string = "dummy" 
#         # ^ *REMEMBER* that this string representation has to "invert" the boundaries because "pval" is minus_log10 values.
# 	if boundary_lowest_possible_value <= pval <= boundary_low: # notice that this boundary is a bit different than the others intervals. But it is correct :-)
# 		value_discretized_numeric = 1
# 		value_discretized_string = "10^-{} < pval <= 10^-{}".format(boundary_low, boundary_lowest_possible_value)
# 	elif boundary_low < pval <= boundary_medium:
# 		value_discretized_numeric = 2
# 		value_discretized_string = "10^-{} < pval <= 10^-{}".format(boundary_medium, boundary_low)
# 	elif boundary_medium < pval <= boundary_high:
# 		value_discretized_numeric = 3
# 		value_discretized_string = "10^-{} < pval <= 10^-{}".format(boundary_high, boundary_medium)

# 	if return_value_type == "numeric_type":
# 		return_value = value_discretized_numeric
# 	elif return_value_type == "string_type":
# 		return_value = value_discretized_string
# 	else:
# 		raise Exception("Internal error: got unexpected argument for 'return_value_type': {}.".format(return_value_type))

# 	return return_value


# TODO: use "bisect" library for the intervals. E.g: bisect.bisect_left([100,300,500,800,1000], p)
def discretize_pval(pval, scale, return_value_type):
	# "scale" is a numpy array generated by linspace - increasing order.
	boundary_lowest_possible_value = scale[0] # this is the lowest possible pvalue.
	boundary_low = scale[1]
	boundary_medium = scale[2]
	boundary_high = scale[3]

	value_discretized_numeric = -99 # TODO: find better initial value!? It MUST be numeric.
	value_discretized_string = "dummy" 
        # ^ *REMEMBER* that this string representation has to "invert" the boundaries because "pval" is minus_log10 values.
	if boundary_lowest_possible_value <= pval <= boundary_low: # notice that this boundary is a bit different than the others intervals. But it is correct :-)
		value_discretized_numeric = 3
		value_discretized_string = "{:.2e} < pval <= {:.2e}".format(boundary_lowest_possible_value, boundary_low)
	elif boundary_low < pval <= boundary_medium:
		value_discretized_numeric = 2
		value_discretized_string = "{:.2e} < pval <= {:.2e}".format(boundary_low, boundary_medium)
	elif boundary_medium < pval <= boundary_high:
		value_discretized_numeric = 1
		value_discretized_string = "{:.2e} < pval <= {:.2e}".format(boundary_medium, boundary_high)

	if return_value_type == "numeric_type":
		return_value = value_discretized_numeric
	elif return_value_type == "string_type":
		return_value = value_discretized_string
	else:
		raise Exception("Internal error: got unexpected argument for 'return_value_type': {}.".format(return_value_type))

	return return_value


def generate_node_attributes(df):
	""" Generate the df_node_attributes """

	df_node_attributes = df # copy
	### Set Cytoscape label
	df_node_attributes['node_label_text'] = df_node_attributes['Original gene set description'].map(set_cytoscape_node_label_text)
	### set node and label dimension. TODO: consider using CYTOSCAPE for positioning labels: "Label Force-Directed Layout: The Label Force-Directed Layout is an automatic label repositioning" 
	#df_node_attributes['node_width'], df_node_attributes['node_height'], df_node_attributes['label_width'] = zip(*df_node_attributes['Original gene set description'].map(set_cytoscape_node_height_and_width))
	df_node_attributes['node_width'], df_node_attributes['node_height'], df_node_attributes['label_width'] = zip(*df_node_attributes['node_label_text'].map(set_cytoscape_node_height_and_width))
	#df_node_attributes['minuslogten_pval_discrete_floor'] = -np.log10(df_node_attributes['Nominal P value']).round(0) # returns integer Series | PASCAL ORIGINAL
	df_node_attributes['within_cluster_min_pval_minuslogten'] = -np.log10(df_node_attributes['Cluster minimum P value'])
	df_node_attributes['minuslogten_pval_discrete'] = -np.log10(df_node_attributes['Nominal P value'])
	df_node_attributes['minuslogten_pval_discrete_floor'] = np.floor(-np.log10(df_node_attributes['Nominal P value'])) # returns integer Series
		# NOTE: by "flooring" "minus_log10" values, we are making the p-values *less* significant. Hence the following in-equality is true: real_pvalue <= [pvalue represented in "minuslogten_pval_discrete_floor"]

	###### BEFORE 2016-01-11 #########
	# ### ***Discretizing the "Cluster minimum P value"***
	# #df_node_attributes['within_cluster_min_pval_minuslogten_floor'] = -np.log10(df_node_attributes['Cluster minimum P value']).round(0) # PASCAL ORIGINAL
	# df_node_attributes['within_cluster_min_pval_minuslogten_floor'] = np.floor(-np.log10(df_node_attributes['Cluster minimum P value'])) # *2015-12 TUNE SUGGESTION*
	# tmp_scale_largest_pval = df_node_attributes['within_cluster_min_pval_minuslogten_floor'].min() # "upper" boundary: "largest values, least significant"
	# 	# Result of using "floor" for the upper/lower boundary: ***results in a LARGER "upper bundary". Hence the p-values can squeezed *below* this value
	# tmp_scale_smallest_pval = df_node_attributes['within_cluster_min_pval_minuslogten_floor'].max() # "lower" boundary: "smallest p-values, most significant"
	# 	# Result of using "floor" for the upper/lower boundary: ***results in a LARGER "upper bundary". Hence the p-values can squeezed *below* this value
	# tmp_scale = np.floor(np.linspace(tmp_scale_largest_pval, tmp_scale_smallest_pval, num=4, endpoint=True)) # *2015-12 TUNE SUGGESTION*
	# 	#Returns num evenly spaced samples, calculated over the interval [start, stop ]. The endpoint of the interval is by default INCLUDED
	# 	# *OBS IMPORTANT* "num=X", where X is an integer, must match how the function "discretize_pval" processes the values.

	# ### Mapping via the function "discretize_pval()"
	# df_node_attributes['within_cluster_min_pval_discrete'] = df_node_attributes['within_cluster_min_pval_minuslogten_floor'].apply(discretize_pval, scale=tmp_scale, return_value_type="numeric_type")
	# 	# column "within_cluster_min_pval_discrete" contains ONLY the values (currently 1,2,3) from discretize_pval()
	# 	#SEE http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.apply.html#pandas.Series.apply
	# df_node_attributes['within_cluster_min_pval_discrete_string'] = df_node_attributes['within_cluster_min_pval_minuslogten_floor'].apply(discretize_pval, scale=tmp_scale, return_value_type="string_type")


	### ***Discretizing the "Cluster minimum P value"***
	tmp_scale_smallest_pval = df_node_attributes['Cluster minimum P value'].min() 
	tmp_scale_largest_pval = df_node_attributes['Cluster minimum P value'].max() 
	tmp_scale = np.linspace(tmp_scale_smallest_pval,tmp_scale_largest_pval, num=4, endpoint=True)
		# REMEMBER THAT "start" and "stop" points are BY DEFAULT included in the sequence. Hence tmp_scale_smallest_pval and tmp_scale_largest_pval will be INCLUDED!

	### Mapping via the function "discretize_pval()"
	df_node_attributes['within_cluster_min_pval_discrete'] = df_node_attributes['Cluster minimum P value'].apply(discretize_pval, scale=tmp_scale, return_value_type="numeric_type")
		# column "within_cluster_min_pval_discrete" contains ONLY the values (currently 1,2,3) from discretize_pval()
		#SEE http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.apply.html#pandas.Series.apply
	df_node_attributes['within_cluster_min_pval_discrete_string'] = df_node_attributes['Cluster minimum P value'].apply(discretize_pval, scale=tmp_scale, return_value_type="string_type")



	# =============================== *TODO* find a better place for this code block ============================
	# It is not nice to have these print statements and global variables like "file_summary" in this block

	### Write summary file
	print "Writing summary file containing information about the 'Cytoscape Node Fill Color Mapping.'"
	print "Please check out the summary mapping file, to be able to construct legends for the network. Path to file: {}".format(file_summary)
	with open(file_summary, 'w') as f:
		f.write("### Cytoscape Node Fill Color Mapping ###\n")
		f.write("attribute mapping variable: within_cluster_min_pval_discrete. This is the within cluster minium pval, transformed into a discrete variable in the interval [1-3].")
		f.write("Largest P-value [min(within_cluster_min_pval_discrete)]: {}\n".format(tmp_scale_largest_pval))
		f.write("Smallest P-value [max(within_cluster_min_pval_discrete)]: {}\n".format(tmp_scale_smallest_pval))
		f.write("Node color intervals: group 1=[{}-{}]; group 2=]{}-{}]; group 3=]{}-{}]\n".format(tmp_scale[0],tmp_scale[1], tmp_scale[1],tmp_scale[2], tmp_scale[2],tmp_scale[3]))
		f.write("Here is how you may setup your lengend:\n")
		f.write("Color 1: P-val <= {} \n".format(tmp_scale[1]))
		f.write("Color 2: P-val <= {} \n".format(tmp_scale[2]))
		f.write("Color 3: P-val <= {} \n".format(tmp_scale[3]))
	# ===========================================================================================================

	return df_node_attributes

def write_node_attribute_file(df):
	### Indexing Panda objects can return two fundamentally different objects: a view or a copy.
	### REF: http://stackoverflow.com/questions/17995328/changing-values-in-pandas-dataframe-doenst-work
	### OBS: any *modification* to the data frame in "view", modifies the "original" data frame
	df_write = df.copy(deep=True) # deep "copy" - *OBS: INEFFICIENT MEMORY USAGE*

	### Columns
	# Original gene set ID
	# Original gene set description
	# Nominal P value
	# False discovery rate
	# Cluster ID
	# Cluster center (boolean)
	# Cluster minimum P value (boolean)
	# Cluster gene set with minimum P value
	# Cluster minimum P value
	# node_label_text
	# node_width
	# node_height
	# label_width
	# minuslogten_pval_discrete_floor
	# within_cluster_min_pval_minuslogten_floor
	# within_cluster_min_pval_discrete
	print "Writing node attribute file..."
	df_write["Cluster ID"] = df_write["Cluster ID"] + 1 # looks better in the output
	df_write.sort_values(['Cluster ID','Nominal P value'], inplace=True)
	df_write.to_csv(file_node_attribute, sep="\t", index=False)
	print "Done"
	
def write_cytoscape_script():
	"""
	Function to write out cytoscape script 

	*IMPORTANT*: Cytoscape does *NOT accept single quotes* - only double quotes.
	"""
	
	print "Writing cytoscape script..."

	f = open(file_cytoscape_script, 'w')

	### Import network
	f.write("""network import file file="{file}" firstRowAsColumnNames=true startLoadRow=1 indexColumnSourceInteraction=1 indexColumnTargetInteraction=2""".format(file=file_network_table)+"\n")

	### Import node attributes
	# REF: http://wiki.cytoscape.org/Cytoscape_3/UserManual/Attributes
	# *OBS*: tab seperated
	f.write(r"""table import file file="{file}" firstRowAsColumnNames=true startLoadRow=1 keyColumnIndex=1 delimiters="\t" """.format(file=file_node_attribute)+"\n")
		### ^ use raw string to inhibit python from setting a "	" character

	### Set layout
	f.write("""layout kamada-kawai EdgeAttribute=Pearson_correlation unweighted=false"""+"\n")
		# Kamada and Kawai (1988): same as "Spring-Embedded Layout"
		# Network nodes are treated like physical objects that repel each other, such as electrons. 
		# The layout algorithm sets the positions of the nodes in a way that minimizes the sum of forces in the network
	# 1) [yWorks - not available through command line] yFiles Organic Layout: organic layout algorithm is a kind of *spring-embedded algorithm*
	# 2) layout force-directed EdgeAttribute=Pearson_correlation unweighted=false
	# ---> NOTE: the applying a WEIGHTED layout algorithm may fail if the network is small or has few edges. This is tricky...

	### Load and apply visual style
	# *OBS*: KEEP STYLE NAME UPDATED - must match name in .xml file
	f.write("""vizmap load file file="{file}" """.format(file=cytoscape_style)+"\n")
	f.write("""vizmap apply styles="DEPICT-style-v1" """+"\n")
	
	### Set view to fit display
	f.write("""view fit content"""+"\n")

	### Export
	f.write("""view export OutputFile="{file_out}" options=PDF""".format(file_out=file_cytoscape_graphics)+"\n")
	f.write("""view export OutputFile="{file_out}" options=PNG""".format(file_out=file_cytoscape_graphics)+"\n")

	f.write("""session save as file="{file_out}" """.format(file_out=file_cytoscape_session)+"\n")

	

	### OPTIONAL: exit script when done.
	if not flag_interactive_cytoscape_session:
		f.write("""command quit"""+"\n")
	

	f.close()


def run_cytoscape_script():
	""" Function to run Cytoscape with the input script generated by this program """

	print "The program is now about to launch cytoscape ({}) and run the script".format(cytoscape_executable)

	list_existing_graphics = glob.glob(file_cytoscape_graphics+"*")
	if list_existing_graphics: # previous results files exists
		print "Any existing graphics will be overwritten/deleted..."
		for elem in list_existing_graphics:
			os.remove(elem)
			print "WARNING: found and deleted existing graphics file {}".format(elem)



	### SHELL VERSION
	#cmd = "{executable} -S {script}".format(executable=cytoscape_executable, script=file_cytoscape_script)
	### EXECUTABLE VERSION
	cmd = [cytoscape_executable, "-S", file_cytoscape_script] # <-- remember the structure of the list an arguments.
		# shlex.split('foo -a -b --bar baz') --> ['foo', '-a', '-b', '--bar', 'baz']
		# subprocess.list2cmdline(cmd) <-- undocumented, not sure it works perfectly

	print "Running command: {}".format(cmd)
	with open(os.devnull) as fnull:
		### SHELL VERSION
		#p = subprocess.Popen(cmd, stdout=fnull, stderr=None, shell=True) # stderr=None --> STDERR goes to the terminal output
		### EXECUTABLE VERSION
		p = subprocess.Popen(cmd, stdout=fnull, stderr=None) # stderr=None --> STDERR goes to the terminal output
	print "Launched Cytoscape session. PID (Process ID): {}".format(p.pid)
	# print "Waiting for Cytoscape to finish..."
	# p.wait()
	# print "Cytoscape done. Graphics are now created."

#############################################################################################
###################################### GET CONFIG FILE ######################################
#############################################################################################

# Get config file delivered as command line argument
parser = argparse.ArgumentParser(description='DEPICT network plot')
parser.add_argument('cfg_file', metavar='DEPICT network plot configuration file', type=str, help='DEPICT network plot configuration file, all user inputs are specified in this file')
args = parser.parse_args()

# Read path to config file
cfg = ConfigParser.ConfigParser() 
	# DEFAULT is "allow_no_value=False", hence *not allowing* empty values in the config file.
	# consider using ConfigParser.SafeConfigParser: SafeConfigParser does not make any attempt to understand the option type
cfg_succes = cfg.read(args.cfg_file) # Attempt to read and parse a *list* of filenames, returning a list of filenames which were successfully parsed. if not files are found, the list will be empty 
if not cfg_succes: # check if list is empty
	print "Could not read config file: {}".format(FILE_CONFIG)
	print "Please make sure that the config file is located in the same directory as the network_plot.py script."
	print "Will exit..."
	sys.exit(0)

file_reconstituted_genesets_matrix = os.path.abspath(cfg.get("RECONSTITUTED GENE SETS","file_reconstituted_genesets_matrix"))
cytoscape_executable = os.path.abspath(cfg.get("CYTOSCAPE",'cytoscape_executable'))
cytoscape_style = os.path.abspath(cfg.get("CYTOSCAPE",'cytoscape_style'))
file_genesetenrichment = os.path.abspath(cfg.get("INPUT PARAMETERS","file_genesetenrichment"))
cutoff_type = cfg.get("INPUT PARAMETERS","cutoff_type")
fdr_cutoff = cfg.get("INPUT PARAMETERS", "fdr_cutoff")
pvalue_cutoff = cfg.getfloat("INPUT PARAMETERS", "pvalue_cutoff") # may be empty string
out = cfg.get("OUTPUT PARAMETERS","output_label")
node_selection = cfg.get("OUTPUT PARAMETERS","node_selection")
network_correlation_cutoff = cfg.getfloat("OUTPUT PARAMETERS",'network_correlation_cutoff') # coerces option to a floating point number
genesetID_network = cfg.get("OUTPUT PARAMETERS","genesetID_network") # may be empty string
flag_interactive_cytoscape_session = cfg.getboolean("OUTPUT PARAMETERS","flag_interactive_cytoscape_session")
flag_network_plot = cfg.getboolean("OUTPUT PARAMETERS","flag_network_plot")



################## PROCESS/VALIDATE config arguments ##################

### Check cytoscape settings
if flag_network_plot: # if true --> RUNNING in "FULL MODE"
	if not os.path.exists(cytoscape_executable):
		raise Exception("ERROR IN CONFIG: Cytoscape executable {} does not exists.".format(cytoscape_executable))
	if not os.path.exists(cytoscape_style):
		raise Exception("ERROR IN CONFIG: Cytoscape style {} does not exists.".format(cytoscape_style))
else: # cluster only mode
	print "Running in 'cluster only' mode: will not generate network plot | flag_network_plot=[{}]".format(flag_network_plot)

### Check cutoff_type | *must be done BEFORE checking fdr or pvalue values*
if not cutoff_type in EXPECTED_CUTOFF_TYPES:
	raise Exception("ERROR IN CONFIG: Got unexpected cutoff_type value '{}'".format(cutoff_type))

### Check cutoff values
if cutoff_type == "fdr":
	list_of_fdr_cutoffs = [x.strip() for x in fdr_cutoff.split(",")] # split on comma and remove whitespace
		# ^^ this will never fail? --> no need to try/catch?
	for fdr in list_of_fdr_cutoffs:
		if not fdr in EXPECTED_DEPICT_FDR_CUTOFFS:
			raise Exception("ERROR IN CONFIG: Got unexpected FDR cutoff '{}'".format(fdr))
elif cutoff_type == "pvalue":
	#if not isinstance(pvalue_cutoff, float): raise Exception("ERROR IN CONFIG: Expected pvalue_cutoff to be of type float. Check the config file.")
		# ^^ TOO MUCH
	pass
else:
	raise Exception("INTERNAL ERROR in source code: 'if statements' related to determining cutoff_type is not updated/synchronized with EXPECTED_CUTOFF_TYPES.")


time_script_start = time.time() # *START TIME*

### *LEGACY from parsing arguments via argparse (command line arguments)* ### 
### The check is not needed when using configparser, because allow_no_value=False is used by Configparser.
# if out is None:
# 	out = os.path.join(os.getcwd(), 'network_plot', 'network_plot')
# 	out = os.path.abspath(out) # convert to absolute path.
# 	# e.g. /Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src/network_plot/network_plot
# elif out == "": # empty string
# 	raise Exception("No output_label label given. Please correct the config file")
# else:
# 	out = os.path.abspath(out) # convert to absolute path. # E.g.:
# 		# os.path.abspath(".") --> '/Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src'
# 		# os.path.abspath("./SOMEDIR") --> '/Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src/SOMEDIR'
# 		# os.path.abspath("SOMEDIR/ANOTHERDIR") --> '/Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src/SOMEDIR/ANOTHERDIR'
# 		# os.path.abspath("/Users/pascaltimshel/") --> '/Users/pascaltimshel'

### Create output directory
out_dir = os.path.dirname(out)
if not os.path.exists(out_dir):
	print "WARNING: out directory does not exist. Will create the directory: {}".format(out_dir)
	os.makedirs(out_dir)


print "Output file-prefix: {}".format(out)


if genesetID_network != "":
	if not isinstance(genesetID_network, str): # just to be sure...
		raise Exception("Value of genesetID_network argument is not of type string.")
	print "genesetID_network argument given: {}".format(genesetID_network)


################## Setting output filenames ##################
if genesetID_network:
	gsID_clean = re.sub(r'[^\w]', '', genesetID_network)  # regex: stripping all symbols from string. keeping only alpha-numeric characters.
	file_network_table = "{out}_genesetID-{gsID}_network_table.{ext}".format(out=out, gsID=gsID_clean, ext="txt")
	file_cytoscape_graphics = "{out}_genesetID-{gsID}_network_diagram".format(out=out, gsID=gsID_clean) # OBS: no extension - cytoscape will do this
	file_cytoscape_session = "{out}_genesetID-{gsID}_cytoscape_session".format(out=out, gsID=gsID_clean) # OBS: no extension - cytoscape will do this
else:
	file_network_table= out + "_network_table.txt"
	file_cytoscape_graphics = out + "_network_diagram" # OBS: no extension - cytoscape will do this
	file_cytoscape_session = out + "_cytoscape_session" # OBS: no extension - cytoscape will do this

file_cytoscape_script = out + "_tmp_cytoscape_script.txt"


file_node_attribute = out + "_nodeattributes.txt" # *OBS* the file CANNOT be loaded by Cytoscape if the extension is ".attrs"
file_genesetenrichment_cluster_result = out + "_cluster_results.txt"

file_summary = out + "_summary.txt" # file contains the discretizing range mapping and MORE




### *KEEP ME UPDATED*
list_of_out_files = [file_network_table,
					file_node_attribute,
					file_genesetenrichment_cluster_result,
					file_summary,
					file_cytoscape_script,
					file_cytoscape_graphics,
					file_cytoscape_session] # USED for checking for existing files

flag_existing_out_files = False
for elem in list_of_out_files:
	if os.path.exists(elem):
		flag_existing_out_files = True

### Make sure that the genotype prefix is correct ###
if flag_existing_out_files:
	ans = ""
	print "WARNING: detected that one or more output file already existed in the directory: {}".format(out_dir)
	print "These files will be overwritten."
	print "You can accept and continue by typing 'yes'."
	print "You can stop and kill the program by typing 'no' (or hit Ctrl-c)."
	while ans != 'yes':
	 	ans = raw_input("Continue ('yes'/'no'): ")
	 	if ans == "no":
	 		sys.exit(0)
	print "Ok let's start..."


#######################################################################################
###################################### READ DATA ######################################
#######################################################################################


################## Geneset Enrichment file ##################
tmp_header_cols = pd.read_csv(file_genesetenrichment, sep="\t", nrows=1).columns[0:len(COLS2READ_GENESETENRICHMENT)] # read only header and first entry | pandas Index
if not (tmp_header_cols == COLS2READ_GENESETENRICHMENT).all():
	raise Exception( "The DEPICT geneset enrichment file did not fit the correct format for the header. The following header fields is expected: [{}]".format("|".join(COLS2READ_GENESETENRICHMENT)) )

df_genesetenrichment_raw = pd.read_csv(file_genesetenrichment, sep="\t", usecols=COLS2READ_GENESETENRICHMENT) # usecols: either column names or position numbers
print "Gene enrichment file read: {}".format(file_genesetenrichment)

### **SUBSET data based on cutoff_type**
df_genesetenrichment = subset_genesetenrichment_data_frame(df_genesetenrichment_raw, cutoff_type)


################## Reconstituted geneset matrix ##################
## Obs: compressed via .gz
## Dimensions: [genes x genesets] --> [19987 x 14463]
## Mem required for full matrix --> 2.3 GB
## *ASSUMPTIONS*: the first column with genes names should have the "cols2read_reconstituted_genesets_matrix_with_rowname_symbol" in the header!

### Columns to read
cols2read_reconstituted_genesets_matrix = df_genesetenrichment['Original gene set ID'] # cols to read gs from enrichment file. 
cols2read_reconstituted_genesets_matrix_with_rowname_symbol = pd.Series([RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL]).append(cols2read_reconstituted_genesets_matrix) # pushing "rownames" symbol header to the cols2read - OTHERWISE IT WILL NOT BE READ.


### Check for existence of all gene sets in the file_reconstituted_genesets_matrix.
# This will *AVOID* the exception "ValueError: 'GeneSet_XYZ' is not in list"
# Reading only the header to check: 
	# 1) presence of all cols2read
	# 2) "check correct file format".
tmp_header_cols = pd.read_csv(file_reconstituted_genesets_matrix, sep="\t", nrows=1, index_col=False, compression='gzip').columns 
	# ^^ read only header and first entry | returns pandas Index
	# *NO INDEX* --> this will include the 'RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL' in the .columns 
bool_gs_not_found = ~cols2read_reconstituted_genesets_matrix_with_rowname_symbol.isin(tmp_header_cols) # inverted boolean
index_tmp_not_found = cols2read_reconstituted_genesets_matrix_with_rowname_symbol[bool_gs_not_found]
if not index_tmp_not_found.empty:
	# Remember: cols2read_reconstituted_genesets_matrix_with_rowname_symbol is a pd.Series() object
	print "ERROR: Could not find all gene sets from the DEPICT gene set enrichment file in the DEPICT reconstituted gene set matrix."
	print "Number of gene sets from enrichment file that could not be found in the DEPICT reconstituted gene set matrix: {}".format(sum(bool_gs_not_found))
	print "List of gene sets not found: ", cols2read_reconstituted_genesets_matrix_with_rowname_symbol[bool_gs_not_found].values
	print "Please check that you are using the correct DEPICT reconstituted gene set matrix or that the DEPICT gene set enrichment file is not corrupted."
	print "The DEPICT version that generated the enrichment file should be the same version as the DEPICT reconstituted gene set matrix."
	print "Will exit the program..."
	sys.exit(0)
else:
	print "file check of file_reconstituted_genesets_matrix passed"


### Reading data
time_start = time.time()
print "Started reading file_reconstituted_genesets_matrix: {}. This may take a few minutes...".format(file_reconstituted_genesets_matrix)
df_reconstituted_genesets = pd.read_csv(file_reconstituted_genesets_matrix, sep="\t", usecols=cols2read_reconstituted_genesets_matrix_with_rowname_symbol, compression='gzip')
	# OBS 1: if an element (e.g. GeneSet_XYZ) in the argument of "usecols" is not in header Pandas will throw an error --> "ValueError: 'GeneSet_XYZ' is not in list"
		# ^^ We can check which columns could not be found
	# Hint: instead of index_col=0, you can use "df_reconstituted_genesets.set_index('-', inplace=True)"
	# Pandas version >=0.16 infer compression automatically: compression : {'gzip', 'bz2', 'infer', None}, default 'infer'
		# For on-the-fly decompression of on-disk data. If 'infer', then use gzip or bz2 if filepath_or_buffer is a string ending in '.gz' or '.bz2', respectively, and no decompression otherwise. Set to None for no decompression.
df_reconstituted_genesets.set_index(RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL, inplace=True) # setting the gene names as index col
df_reconstituted_genesets.index.name = "ENSG" # just to make it nice
time_elapsed = time.time() - time_start
print "Elapsed time for reading DEPICT reconstituted genesets matrix: {:.2f} sec".format(time_elapsed)
print "Dimension of df_reconstituted_genesets: ", df_reconstituted_genesets.shape
# OBSERVATION: KEGG_DRUG_METABOLISM_CYTOCHROME_P450 does not exists in the "file_reconstituted_genesets_matrix" file BUT ONLY in the "file_genesetenrichment"

# ### CHECK of missing genesets 
# # *OBS* THIS IS ACTUALLY NOT NECESSARY because Pandas will throw a "ValueError: XXX is not in list" if a gene set is in col2read, but *NOT* in the geneset matrix
# #assert(df_reconstituted_genesets.shape[1]==len(cols2read_reconstituted_genesets_matrix)) # ## REQUIRE that all columns of interest where loaded correctly --> maybe this is too much.

##################################################################################################
######################## Safety check for genesetID_network (geneset inset) ######################
##################################################################################################

if genesetID_network != "": # only check if argument is supplied
	if not df_reconstituted_genesets.columns.isin([genesetID_network]).any(): # we need to make sure the genesetID is in the enrichment file
		raise Exception("Value of genesetID_network argument '{}' is either not contained in 1) the DEPICT reconstituted geneset matrix or; 2) the FDR significant gene sets in the enrichment file. Please ensure you specified a valid identifier".format(genesetID_network))

##################################################################################################
###################################### Affinity Propagation ######################################
##################################################################################################

################## Calculate similarity matrix ##################
# using with correlation as similarity measurement

### Numpy approach --> significantly faster than Pandas
#numpy.cov(m) --> not "normalized". Remember: there is a difference between *covariance* and *correlation*
#numpy.corrcoef(X) --> correlation coefficients. X: A 1-D or 2-D array containing multiple variables and observations.  Each row of x represents a variable, and each column a single observation of all those variables. "
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.corrcoef.html
#matrix_corr = np.corrcoef(df_reconstituted_genesets, rowvar=0)
	# ^^ USE THIS INSTEAD OF TRANSPOSING: if *rowvar* is ZERO, then each column represents a variable, while the rows contain observations.
	
### Pandas approach
# http://pandas.pydata.org/pandas-docs/stable/computation.html#correlation
# df.corr() # Compute pairwise correlation of columns, excluding NA/null values
#matrix_corr = df_reconstituted_genesets.corr(method='pearson', min_periods=df_reconstituted_genesets.shape[0]) # --> Note that "min_periods=df_reconstituted_genesets_t.shape[0]" should not be needed, but it ensures that no "NaN" values are present to give weird results

### METHOD IN USE
print "Calculating correlation matrix between gene sets... This will be used as similarity matrix for Affinity Propagation"
#matrix_corr = np.corrcoef(df_reconstituted_genesets, rowvar=0) # NUMPY approach | slighty faster than Pandas
matrix_corr = df_reconstituted_genesets.corr(method='pearson', min_periods=df_reconstituted_genesets.shape[0]) # PANDAS approach | Note that "min_periods=df_reconstituted_genesets_t.shape[0]" should not be needed, but it ensures that no "NaN" values are present to give weird results
print "Dimension of correlation matrix: [{} x {}]".format(matrix_corr.shape[0], matrix_corr.shape[1])

################## Running AP ##################

#sklearn.cluster.AffinityPropagation(damping=0.5, max_iter=200, convergence_iter=15, copy=True, preference=None, affinity='euclidean', verbose=False)
af_obj = AffinityPropagation(affinity = 'precomputed', max_iter=10000, convergence_iter=1000) # using almost only default parameters
print "Affinity Propagation parameters:"
for param, val in af_obj.get_params().items():
	print "\t{}: {}".format(param, val)
print "Perfoming Affinity Propagation.."
af = af_obj.fit(matrix_corr)
n_iter = af.n_iter_
print "Affinity Propagation done"
print "Number of iterations used: {}".format(n_iter)

### Saving labels and centers
cluster_centers_indices = af.cluster_centers_indices_  # array, shape (n_clusters, n_features) | cluster center (boolean)s ("exemplars")
													   # cluster_centers_indices take on values in the range {0...n_samples-1}
labels = af.labels_ # array, shape (n_samples,) | Get the "labels"/assignments of each data point to a cluster index
					# labels take on values in the range: {0...n_clusters-1}



### Display some stats
n_clusters = len(cluster_centers_indices)
print('Estimated number of clusters: %d' % n_clusters)

######################################################################################################
###################################### Creating main data frame ######################################
######################################################################################################

df = add_cluster_results_to_data_frame(df_genesetenrichment)


##########################################################################################
###################################### Write output ######################################
##########################################################################################

#####################################################
################## Cluster results ##################
#####################################################

write_genesetenrichment_cluster_result_file(df)


#############################################################
####################### Network table #######################
#############################################################

if flag_network_plot:
	if genesetID_network:
		selected_clusterID = df.ix[genesetID_network,'Cluster ID'] # returns integer/scalar
		df_network_table = df[df['Cluster ID']==selected_clusterID] # extract observations with selected Cluster ID
		write_network_table_file(df_network_table, file_out=file_network_table)
	else:
		if node_selection == "cluster_center":
			df_network_table = df[df["Cluster center (boolean)"]==True]
		elif node_selection == "cluster_min_pval":
			df_network_table = df[df["Cluster minimum P value (boolean)"]==True]
		elif node_selection == "all":
			df_network_table = df
		else:
			raise Exception("Got unsupported node_selection argument")
		write_network_table_file(df_network_table, file_out=file_network_table)


#####################################################
################## Attribute file ###################
#####################################################

if flag_network_plot: # full mode (*NOT* "cluster only" mode)
	### Generate df_node_attributes (+more?)
	df_node_attributes = generate_node_attributes(df)
	
	### Write attribute file
	write_node_attribute_file(df_node_attributes)

	### Write and run cytoscape script
	write_cytoscape_script()
	run_cytoscape_script()



###################################### FINISH ######################################
time_script_elapsed = time.time() - time_script_start
print "RUNTIME: {:.1f} sec ({:.1f} min)".format(time_script_elapsed,time_script_elapsed/60)
print "====== network_plot.py is finished ======"

