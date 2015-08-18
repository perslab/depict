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

######################################  ######################################

#/usr/libexec/java_home -V

###################################### USAGE ######################################

### Full reconstituted
#python network_plot.py --file_genesetenrichment /Users/pascaltimshel/Dropbox/0_Work/DEPICT/DEPICT_scripts_PT/data_network/EA2_EduYears_pooled_Nweighted_single_gc_genesetenrichment.txt --file_reconstituted_genesets_matrix /Users/pascaltimshel/Dropbox/0_Work/DEPICT/DEPICT_scripts_PT/data_reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz --node_selection cluster_center
### Test set
#python network_plot.py --file_genesetenrichment /Users/pascaltimshel/Dropbox/0_Work/DEPICT/DEPICT_scripts_PT/data_network/EA2_EduYears_pooled_Nweighted_single_gc_genesetenrichment.txt --file_reconstituted_genesets_matrix /Users/pascaltimshel/Dropbox/0_Work/DEPICT/DEPICT_scripts_PT/data_reconstituted_genesets/DEPICT_matrix_reconstituted_genesets_EA2_EduYears_297.txt --node_selection cluster_center --genesetID_network GO:0048813

### Test cytoscape script (standalone)
#/Applications/Cytoscape_v3.2.1/cytoscape.sh -S /Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src/network_plot_cytoscape_script.txt




###################################### CONSTANTS ######################################
COLS2READ_GENESETENRICHMENT = ["Original gene set ID", "Original gene set description", "Nominal P value", "False discovery rate"]

RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL = "-" # update this symbol if the header symbol of the gene names changes


NETWORK_CORRELATION_CUTOFF = 0.3 # GREATHER OR EQUAL TO | this parameter is used to determine the cutoff of when to draw edges between nodes.


CYTOSCAPE_EXECUTABLE = "/Applications/Cytoscape_v3.2.1/cytoscape.sh" # path to cytoscape shell launcher (shell script)
CYTOSCAPE_STYLE = "/Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src/network_plot_CytoscapeStyle_v1.xml" # path to XML file with cytoscape style
	# OBS: note that the cytoscape style name must match the one in this script. (E.g. DEPICT-style-v1)
	# The style name can be found in the top part of the XML file: <visualStyle name="DEPICT-style-v1">

#######################################################################################
###################################### FUNCTIONS ######################################
#######################################################################################

def ParseArguments():
	""" Function to parse commandline arguments """
	arg_parser = argparse.ArgumentParser(description="DEPICT network plot")
	arg_parser.add_argument("--file_genesetenrichment", help="DEPICT geneset enrichment file [tab seperated]", required=True)
	arg_parser.add_argument("--file_reconstituted_genesets_matrix", help="DEPICT reconstituted geneset matrix file [tab seperated]. The file can be compressed (.gz) or uncompressed.", required=True)
	arg_parser.add_argument("--out", help="""Output path (absolute or relative) INCLUDING a file-prefix. 
											 E.g. /Users/john/DEPICT_results/my_phenotype_fileprefix 
											 *Default output path* is the current working directory (executing terminal path) with the file prefix 'network_plot'  """)
	

	### Consider using subparses
	# subparsers = arg_parser.add_subparsers(dest='subcommand', #---> 'dest' is where the varible is save: args.subcommand/arg_parser.subcommand
	# 								   title='subcommands in this script',
	# 								   description='valid subcommands. set subcommand after main program required arguments',
	# 								   help='You can get additional help by writing <program-name> <subcommand> --help')
	## Subparsers
	#arg_subparser_network_all_enriched_genesets = subparsers.add_parser('network_all_enriched_genesets')
	#arg_subparser_network_selected_enriched_genesets = subparsers.add_parser('network_selection_enriched_genesets')


	arg_parser.add_argument("--node_selection", help="""Choose which nodes to visualize in cytoscape network. 
													The option determines how the network_table file is written. 
													cluster_center: plot cluster representatives genesets (N_nodes=N_clusters); 
													cluster_min_pval: plot within-cluster mimimum p-value genesets (N_nodes=N_clusters); 
													all: plot all significant gene set (N_nodes=N_significant_genesets)""", 
													choices=['cluster_center', 'cluster_min_pval', 'all'], default='cluster_center')
	
	#TODO: allow a parameter for the correlation cutoff
	
	arg_parser.add_argument("--genesetID_network", help="""Argument value must be a valid genesetID ('Original gene set ID'). 
															If this argument is supplied, then an additional geneset network will be created for the specified genesetID """)
	#TODO: check the validity of the genesetID
	#TODO: implement a function to allow multiple genesetIDs as arguments


	#arg_parser.add_argument("--XXX", help="XXX; [default is false] ", action='store_true')
	#arg_parser.add_argument("--XXX", type=int, help="XXX", default=0)
	
	args = arg_parser.parse_args()
	return args


def add_cluster_results_to_data_frame(df):
	""" 
	Function will add columns with clustering results to data frame
	The data frame must contain the column 'Original gene set ID' (data frame for the DEPICT enrichments)
	"""

	df.set_index('Original gene set ID', drop=False, inplace=True) # do not delete columns after.
	df.index.name = 'idx' # just to be nice... not important

	### Initialyze columns with APPROPRIATE default values
	df['clusterID'] = np.nan # will become integer later
	df['cluster_center'] = False
	df['cluster_min_pval'] = False

	df['within_cluster_min_pval_gsID'] = np.nan
	df['within_cluster_min_pval'] = np.nan


	################## Assigning cluster labels to data frame ##################
	### using Pandas match
	#http://pandas.pydata.org/pandas-docs/stable/comparison_with_r.html#match
	#http://stackoverflow.com/a/15866201
	#pd.Series(pd.match(s,[2,4],np.nan)) # <-- *CONSIDER USING THIS*

	### Using index
	# *by indexing df according to "df_reconstituted_genesets.columns" we order to data frame correctly when assignning the labels*
	# IMPORTANT --> REQUIRES index to be 'Original gene set ID'
	#df.ix[['GO:0003712', 'GO:0003713'],:]
	df.ix[df_reconstituted_genesets.columns,'clusterID'] = labels

	################## Assigning cluster_center and within_cluster_min_pval ##################

	for k in range(n_clusters): # looping over clusters
		### REMEMBER ### 
		# We are mapping *FROM* column index/names in df_reconstituted_genesets, that was used for similarity calculations, *TO* df/df_genesetenrichment to get descriptions (and p-values) of Gene Sets.
		# We cannot just used the cluster_center_indices as index, because these are *NOT* in sync with the df/df_genesetenrichment
		
		### saving cluster centers
		cluster_center_ID = df_reconstituted_genesets.columns[cluster_centers_indices[k]] # cluster_center_ID --> a string
		df.ix[cluster_center_ID, 'cluster_center'] = True
		
		class_members_ID = df_reconstituted_genesets.columns[labels == k] # array of strings/genesetsIDs. selecting gene sets names that are members of the current cluster
		
		### saving class member with the lowest p-value (including the cluster center observation)
		idx_min_pval = df.ix[class_members_ID, "Nominal P value"].idxmin()
		df.ix[idx_min_pval, 'cluster_min_pval'] = True
		

		### saving more...
		df.ix[class_members_ID, 'within_cluster_min_pval'] = df.ix[idx_min_pval, "Nominal P value"] # scalar
		df.ix[class_members_ID, 'within_cluster_min_pval_gsID'] = df.ix[idx_min_pval, "Original gene set ID"] # scalar


	## safety check
	assert( sum(df["cluster_center"])==n_clusters )
	assert( sum(df["cluster_min_pval"])==n_clusters )

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
			c) convert string to centence case (all lowercase, except first letter)
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
		tmp = pattern_PPI.sub("protein complex", string) # strip pattern
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
			if corr >= NETWORK_CORRELATION_CUTOFF:
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
	### Columns
	# Original gene set ID
	# Original gene set description
	# Nominal P value
	# False discovery rate
	# clusterID
	# cluster_center
	# cluster_min_pval
	# within_cluster_min_pval_gsID
	# within_cluster_min_pval
	print "Writing geneset enrichment cluster result file..."
	df.sort(['clusterID','Nominal P value'], inplace=True)
	df.to_csv(file_genesetenrichment_cluster_result, sep="\t", index=False)
	print "Done"


def discretize_pval(pval, scale):
	# "scale" is a numpy array generated by linspace - increasing order.
	#*OBS*: we do not use the first [lowest] element of "scale"
	boundary_low = scale[1]
	boundary_medium = scale[2]
	boundary_high = scale[3]

	value_discretized = -99
	if 0 <= pval <= boundary_low:
		value_discretized = 1
	elif boundary_low < pval <= boundary_medium:
		value_discretized = 2
	elif boundary_medium < pval <= boundary_high:
		value_discretized = 3

	return value_discretized

def write_node_attribute_file(df):
	### Columns
	# Original gene set ID
	# Original gene set description
	# Nominal P value
	# False discovery rate
	# clusterID
	# cluster_center
	# cluster_min_pval
	# within_cluster_min_pval_gsID
	# within_cluster_min_pval
	# node_label_text
	# node_width
	# node_height
	# label_width
	# minuslogten_pval_discrete
	# within_cluster_min_pval_minuslogten
	# within_cluster_min_pval_discrete

	print "Writing node attribute file..."
	df.sort(['clusterID','Nominal P value'], inplace=True)
	df.to_csv(file_node_attribute, sep="\t", index=False)
	print "Done"
	
def write_cytoscape_script():
	"""
	Function to write out cytoscape script 

	*IMPORTANT*: Cytoscape does *NOT accept single quotes* - only double quotes.
	"""
	
	print "Writing cytoscape script..."

	f = open(file_cytoscape_script, 'w')

	### Import network
	f.write("""network import file file="{file}" firstRowAsColumnNames=true startLoadRow=2 indexColumnSourceInteraction=1 indexColumnTargetInteraction=2""".format(file=file_network_table)+"\n")

	### Import node attributes
	# REF: http://wiki.cytoscape.org/Cytoscape_3/UserManual/Attributes
	# *OBS*: tab seperated
	f.write(r"""table import file file="{file}" firstRowAsColumnNames=true startLoadRow=2 keyColumnIndex=1 delimiters="\t" """.format(file=file_node_attribute)+"\n")
		### ^ use raw string to inhibit python from setting a "	" character

	### Set layout
	f.write("""layout kamada-kawai EdgeAttribute=Pearson_correlation unweighted=false"""+"\n")
		# Kamada and Kawai (1988)
		# Network nodes are treated like physical objects that repel each other, such as electrons. 
		# The layout algorithm sets the positions of the nodes in a way that minimizes the sum of forces in the network
	# 1) [yWorks - not available through command line] yFiles Organic Layout: organic layout algorithm is a kind of *spring-embedded algorithm*
	# 2) layout force-directed EdgeAttribute=Pearson_correlation unweighted=false


	### Load and apply visual style
	# *OBS*: KEEP STYLE NAME UPDATED - must match name in .xml file
	f.write("""vizmap load file file="{file}" """.format(file=CYTOSCAPE_STYLE)+"\n")
	f.write("""vizmap apply styles="DEPICT-style-v1" """+"\n")
	
	### Set view to fit display
	f.write("""view fit content"""+"\n")

	### Export
	f.write("""view export OutputFile="{file_out}" options=PDF""".format(file_out=file_cytoscape_graphics)+"\n")
	f.write("""view export OutputFile="{file_out}" options=PNG""".format(file_out=file_cytoscape_graphics)+"\n")

	f.close()


def run_cytoscape_script():
	""" Function to run Cytoscape with the input script generated by this program """

	print "The program is now about to launch cytoscape ({}) and run the script".format(CYTOSCAPE_EXECUTABLE)

	list_existing_graphics = glob.glob(file_cytoscape_graphics+"*")
	if list_existing_graphics: # previous results files exists
		print "Any existing graphics will be overwritten/deleted..."
		for elem in list_existing_graphics:
			os.remove(elem)
			print "WARNING: found and deleted existing graphics file {}".format(elem)


	cmd = "{executable} -S {script}".format(executable=CYTOSCAPE_EXECUTABLE, script=file_cytoscape_script)
	print "Running command: {}".format(cmd)
	with open(os.devnull) as fnull:
		p = subprocess.Popen(cmd, stdout=fnull, stderr=None, shell=True) # stderr=None --> STDERR goes to the terminal output
	print "Waiting for process to finish. Please wait..."
	p.wait()
	print "Cytoscape done. Graphics are now created."

#######################################################################################
###################################### ARGUMENTS ######################################
#######################################################################################

args = ParseArguments()

file_genesetenrichment = args.file_genesetenrichment
file_reconstituted_genesets_matrix = args.file_reconstituted_genesets_matrix
node_selection = args.node_selection


out = args.out # complete pathname INCLUDING FILE PREFIX
if out is None:
	out = os.path.join(os.getcwd(), 'network_plot')
	# e.g. /Users/pascaltimshel/Dropbox/0_Projects/git/DEPICT/src/network_plot

### check that output dir is writable
if not os.access(os.path.dirname(out), os.W_OK):
	raise Exception("Output path: {} is not writable. Please fix this...".format(os.path.dirname(out)))

print "Output file-prefix: {}".format(out)

genesetID_network = args.genesetID_network # if not argument is not specified in commandline, the value will be None
if genesetID_network:
	if not isinstance(genesetID_network, str): # just to be sure... MAYBE OVERKILL.
		raise Exception("Value of genesetID_network argument is not of type string.")
	print "genesetID_network argument given: {}".format(genesetID_network)


################## Setting output filenames ##################
file_network_table = out + ".network_table.txt"
file_network_table_genesetID = "{out}_genesetID-{gsID_clean}.network_table.{ext}".format(out=out, gsID_clean=re.sub(r'[^\w]', '', genesetID_network), ext="txt") # regex: stripping all symbols from string. keeping only alpha-numeric characters.

file_node_attribute = out + "_nodeattributes.attrs.txt" # *OBS* the file CANNOT be loaded by Cytoscape if the extension is ".attrs"
file_genesetenrichment_cluster_result = out + "_cluster_result.txt"

file_summary = out + "_summary.txt" # file contains the discretizing range mapping and MORE

file_cytoscape_script = out + "_cytoscape_script.txt"
file_cytoscape_graphics = out + "network_graphic" # OBS: no extension

#######################################################################################
###################################### READ DATA ######################################
#######################################################################################


################## Geneset Enrichment file ##################
df_genesetenrichment = pd.read_csv(file_genesetenrichment, sep="\t", usecols=COLS2READ_GENESETENRICHMENT) # usecols: either column names or position numbers
print "Read gene enrichment file"
#print df_genesetenrichment.head()

### Subset data based on FDR
df_genesetenrichment = df_genesetenrichment[(df_genesetenrichment['False discovery rate']=="<0.01") | (df_genesetenrichment['False discovery rate']=="<0.05")]
print "Found {} number of gene sets with 'FDR<0.01' or 'FDR<0.5'. These gene sets will be used for further analysis".format(df_genesetenrichment.shape[0])

################## Reconstituted geneset matrix ##################
## Obs: compressed via .gz
## Dimensions: [genes x genesets] --> [19987 x 14463]
## Mem required for full matrix --> 2.3 GB
## *ASSUMPTIONS*: the first column with genes names should have the "cols2read_reconstituted_genesets_matrix_with_rowname_symbol" in the header!

### Columns to read
cols2read_reconstituted_genesets_matrix = df_genesetenrichment['Original gene set ID'] # cols to read gs from enrichment file. 
cols2read_reconstituted_genesets_matrix_with_rowname_symbol = pd.Series([RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL]).append(cols2read_reconstituted_genesets_matrix) # pushing "rownames" symbol header to the cols2read - OTHERWISE IT WILL NOT BE READ.

### Reading data
time_start = time.time()
print "Started reading file_reconstituted_genesets_matrix. This may take a few minutes..."
df_reconstituted_genesets = pd.read_csv(file_reconstituted_genesets_matrix, sep="\t", usecols=cols2read_reconstituted_genesets_matrix_with_rowname_symbol)
	# OBS 1: if an element (e.g. GeneSet_XYZ) in the argument of "usecols" is not in header Pandas will throw an error --> "ValueError: 'GeneSet_XYZ' is not in list"
		# ^^ We can check which columns could not be found
	# Hint: instead of index_col=0, you can use "df_reconstituted_genesets.set_index('-', inplace=True)"
	# Pandas infer compression automatically
df_reconstituted_genesets.set_index(RECONSTITUTED_GENESETS_MATRIX_ROWNAME_SYMBOL, inplace=True) # setting the gene names as index col
df_reconstituted_genesets.index.name = "ENSG" # just to make it nice
time_elapsed = time.time() - time_start
print "Elapsed time for reading DEPICT reconstituted genesets matrix: {:.2f} sec".format(time_elapsed)
print "Dimension of df_reconstituted_genesets: ", df_reconstituted_genesets.shape
# OBSERVATION: KEGG_DRUG_METABOLISM_CYTOCHROME_P450 does not exists in the "file_reconstituted_genesets_matrix" file BUT ONLY in the "file_genesetenrichment"

### CHECK of missing genesets
#assert(df_reconstituted_genesets.shape[1]==len(cols2read_reconstituted_genesets_matrix)) # ## REQUIRE that all columns of interest where loaded correctly --> maybe this is too much.
bool_gs_not_found = ~cols2read_reconstituted_genesets_matrix.isin(df_reconstituted_genesets) # inverted boolean
print "Number of gene sets from enrichment file that could not be found in DEPICT matrix: {}".format(sum(bool_gs_not_found))
df_tmp_not_found = cols2read_reconstituted_genesets_matrix[bool_gs_not_found]
if not df_tmp_not_found.empty:
	print "List: ", cols2read_reconstituted_genesets_matrix[bool_gs_not_found]

##################################################################################################
######################## Safety check for genesetID_network (geneset inset) ######################
##################################################################################################

if genesetID_network is not None: # only check if argument is supplied
	if not df_reconstituted_genesets.columns.isin([genesetID_network]).any(): # we need to make sure the genesetID is in the enrichment file
		raise Exception("Value of genesetID_network argument '{}' is either not contained in DEPICT reconstituted geneset matrix or the enrichment file. Please ensure you specified a valid identifier".format(genesetID_network))

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
cluster_centers_indices = af.cluster_centers_indices_  # array, shape (n_clusters, n_features) | cluster centers ("exemplars")
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

#############################################################
################## Network table - primary ##################
#############################################################

if node_selection == "cluster_center":
	df_network_table = df[df["cluster_center"]==True]
elif node_selection == "cluster_min_pval":
	df_network_table = df[df["cluster_min_pval"]==True]
elif node_selection == "all":
	df_network_table = df
else:
	raise Exception("Got unsupported node_selection argument")
write_network_table_file(df_network_table, file_out=file_network_table)


#######################################################################
################## Network table - genesetID_network ##################
#######################################################################

if genesetID_network:
	selected_clusterID = df.ix[genesetID_network,'clusterID'] # returns integer/scalar
	df_network_table = df[df['clusterID']==selected_clusterID] # extract observations with selected clusterID
	write_network_table_file(df_network_table, file_out=file_network_table_genesetID)


#####################################################
################## Cluster results ##################
#####################################################

write_genesetenrichment_cluster_result_file(df)

#####################################################
################## Attribute file ###################
#####################################################

df_node_attributes = df # copy
### Set Cytoscape label
df_node_attributes['node_label_text'] = df_node_attributes['Original gene set description'].map(set_cytoscape_node_label_text)
### set node and label dimension. TODO: consider using CYTOSCAPE for positioning labels: "Label Force-Directed Layout: The Label Force-Directed Layout is an automatic label repositioning" 
#df_node_attributes['node_width'], df_node_attributes['node_height'], df_node_attributes['label_width'] = zip(*df_node_attributes['Original gene set description'].map(set_cytoscape_node_height_and_width))
df_node_attributes['node_width'], df_node_attributes['node_height'], df_node_attributes['label_width'] = zip(*df_node_attributes['node_label_text'].map(set_cytoscape_node_height_and_width))
df_node_attributes['minuslogten_pval_discrete'] = -np.log10(df_node_attributes['Nominal P value']).round(0) # returns integer Series

### discretizing the "within_cluster_min_pval"
df_node_attributes['within_cluster_min_pval_minuslogten'] = -np.log10(df_node_attributes['within_cluster_min_pval']).round(0)
tmp_scale_min = df_node_attributes['within_cluster_min_pval_minuslogten'].min()
tmp_scale_max = df_node_attributes['within_cluster_min_pval_minuslogten'].max()
tmp_scale = np.linspace(tmp_scale_min, tmp_scale_max, num=4, endpoint=True).round(0)
	#Returns num evenly spaced samples, calculated over the interval [start, stop ]. The endpoint of the interval is by default INCLUDED




df_node_attributes['within_cluster_min_pval_discrete'] = df_node_attributes['within_cluster_min_pval_minuslogten'].apply(discretize_pval, scale=tmp_scale)
	#SEE http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.apply.html#pandas.Series.apply


### Write summary file
print "Writing summary file containing information about the 'Cytoscape Node Fill Color Mapping.'"
print "Please check out the summary mapping file, to be able to construct legends for the network. Path to file: {}".format(file_summary)
with open(file_summary, 'w') as f:
	f.write("### Cytoscape Node Fill Color Mapping ###\n")
	f.write("attribute mapping variable: within_cluster_min_pval_discrete. This is the within cluster minium pval, transformed into a discrete variable in the interval [1-3].")
	f.write("min(within_cluster_min_pval_discrete): {}\n".format(tmp_scale_min))
	f.write("max(within_cluster_min_pval_discrete): {}\n".format(tmp_scale_max))
	f.write("ranges 1=[{}-{}]; 2=]{}-{}]; 3=]{}-{}]\n".format(0,tmp_scale[1], tmp_scale[1],tmp_scale[2], tmp_scale[2],tmp_scale[3]))
	f.write("color 1: P-val <= {} \n".format(tmp_scale[1]))
	f.write("color 2: P-val <= {} \n".format(tmp_scale[2]))
	f.write("color 3: P-val <= {} \n".format(tmp_scale[3]))


### Write attribute file
write_node_attribute_file(df_node_attributes)



##################################################################################################
############################### Write and run cytoscape script ###################################
##################################################################################################

write_cytoscape_script()


run_cytoscape_script()

