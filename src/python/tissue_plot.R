############### SYNOPSIS ###################
# CLASS: script
# PURPOSE: generate tissue enrichment plot for DEPICT results

### === DESCRIPTION === ###
# This script reads DEPICT tissue enrichment results and generates plots.
# The script automatically detects whether the format is in "Gene Network" or "GTEx" format.

### === USAGE === ###
### Recommended invocation of script (prints to STDOUT)
# Rscript tissue_plot.R <file_tissue_enrichment> [output_prefix]
# Example call #1: Rscript tissue_plot.R ../example/ldl_teslovich_nature2010_tissueenrichment.txt ldl_teslovich
# Example call #2: Rscript tissue_plot.R ../example/ldl_teslovich_nature2010_tissueenrichment_gtex.txt ldl_teslovich

### Alternative invocation of script (writes STDOUT output to file). (Remember to quote the arguments):
# R CMD BATCH --no-save --no-restore '--args <file_tissue_enrichment>' tissue_plot.R <R_OUTPUT_FILENAME>.out

### === ARGUMENT (input) === ### 
# file_tissue_enrichment (required): 
# - an absolute or relative filepath to the DEPCIT tissue enrichment file.
# - the file can be either in "GTEx" or "Gene Network" format.
# output_prefix (optional):
# - a prefix to the output filenames. 
# - E.g. if output_prefix is 'ldl_teslovich', the output files would be tissue_plot_ldl_teslovich_PLOTSPECIFICSUFFIX.pdf

### === OUTPUT === ### 
# The script writes .pdf file plots to the DIRECTORY FROM WHERE THE SCRIPT IS CALLED (working directory of terminal session).
#TODO: write more information about the .pdf files


### === STEPS performed in script - behind the scenes === ###
# 1) Read DEPICT tissue enrichment results.
# 2) Auto-detect format type
# 3) Process files and build plot variables
# 4) Construct plots via ggplot2
# 5) Export plots to pdf files

################## LOAD PACKAGES and ARGUMENTS #####################
#suppressPackageStartupMessages()
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid)) # needed for "unit" function

rm(list=ls())
#options(warn=-1) # turn off warnings globally

#options(echo=TRUE) # Remove this comment if for a "behind-the-scenes" output level
print(sprintf("WORKING DIRECTORY (output for plots): %s", getwd()))

######################################################
################# GLOBAL variables ###################

flag.genenetwork.plot_extended_mesh_categories <- TRUE # boolean value | determines if mesh categories "Cells", "Tissues" and "dummy_mixed_categories" will be plotted
n.blank_spacers <- 4 # integer | number of empty bars between "groups"

### =========== INTERNAL PARAMETERS =========== ###
colnames.genenetwork <-  make.names(c("MeSH term","Name","MeSH first level term","MeSH second level term","Nominal P value","False discovery rate"))
colnames.gtex <- make.names(c("Name","Nominal P value","False discovery rate"))


# OBS: these variable is used for GROUPING and SORTING the data frame
group_variable.genenetwork.general <- make.names("MeSH first level term") # col_no = 3
group_variable.genenetwork.cells_and_tissues <- make.names("MeSH second level term") # col_no = 4 | *THIS IS USED TO SET CORRECT GROUP IN function.genenetwork.set_variables_and_split()*
group_variable.gtex <- make.names("Name")


######################################################
##################### FUNCTIONS ######################

function.input_procces <- function(df, cols2keep) {
  ### USE: UNIVERSIAL
  ### Subsetting data ###
  df.res <- subset(df, select=cols2keep)
  ### SORTING
  #df.res <- df.res[order(df.res[,sort.colname]),] # IMPORTANT: sorting data frame
  # ^ <--- no longer needed to sort here. The functions function.*.set_variables* does this!
  return(df.res)
}

function.genenetwork.set_variables_and_split <- function(df) {
  ### USE: Gene Network
  ### STEPS: 
  # 1) split df into "* System", "Cells", "Tissues" and "dummy_mixed_categories"
  # 2) clean "group" names
  # 3) adding new variables
  # 4) return *list* of data frames
  
  df.res <- df
  ### Set "significance" variable
  df.res$FDR.significant <- with(df.res, ifelse( (False.discovery.rate=="<0.01"|False.discovery.rate=="<0.05"), TRUE, FALSE))
  
  ### Set "group" variable [DO THIS BEFORE PARTITIONING]
  df.res$group <- df.res[,group_variable.genenetwork.general]
  
  ### Matching categories [grepl returns a logical vector]
  bool.system <- grepl("System",df.res[,group_variable.genenetwork.general],ignore.case=T)
  bool.cells <- grepl("Cells",df.res[,group_variable.genenetwork.general],ignore.case=T)
  bool.tissues <- grepl("Tissues",df.res[,group_variable.genenetwork.general],ignore.case=T)
  bool.mixed <- !(bool.system | bool.cells | bool.tissues) # "uncategorized"
  ##TODO: check for mutual exclusivity
  
  ### Partitioning into data frames
  df.system <- subset(df.res, bool.system)
  df.cells <- subset(df.res, bool.cells)
  df.tissues <- subset(df.res, bool.tissues)
  df.mixed <- subset(df.res, bool.mixed)
  
  ### *Cleaning names in "group" variable*  [AFTER PARTITIONING]
  df.system[,"group"] <- with(df.system, sub("(\\s*Systems\\s*)|(\\s*System\\s*)", "", group, ignore.case=TRUE, perl=TRUE)) #|(Systems\\s+)
  # e.g. Nervous System --> Nervous
  # note that there exists a mesh category called "Hemic and Immune Systems". Therefore we need to match "Systems". Remember to place that capture group first in the pattern
  # sub("(\\s*Systems\\s*)|(\\s*System\\s*)", "X", df.system$MeSH.first.level.term, ignore.case=TRUE, perl=TRUE) #|(Systems\\s+)
  
  #***OBS*** --> the group_variable.genenetwork.cells/group_variable.genenetwork.tissues is used INSTEAD OF GROUP
  # ^ this is done because the PARTITIONING is done on the "MeSH first level term" (group_variable.genenetwork)
  df.cells[,"group"] <- sub("\\s*Cells\\s*", "", df.cells[,group_variable.genenetwork.cells_and_tissues], ignore.case=TRUE, perl=TRUE)
  #"Cells": Stem Cells --> Stem; Antigen-Presenting Cells --> Antigen-Presenting
  df.tissues[,"group"] <- sub("\\s*Tissue\\s*", "", df.tissues[,group_variable.genenetwork.cells_and_tissues], ignore.case=TRUE, perl=TRUE)
  #"Tissue": Lymphoid Tissue --> Lymphoid; Epithelium -(*IRREGULAR*)-> Epithelium
  
  
  ### ***SORTING***  
  # It is very important to sort the data frame(s). The barplot breaks/labels rely on the df's being sorted.
  df.system <- df.system[order(df.system[,"group"]),]
  df.cells <- df.cells[order(df.cells[,"group"]),]
  df.tissues <- df.tissues[order(df.tissues[,"group"]),]
  df.mixed <- df.mixed[order(df.mixed[,"group"]),]
  
  list.dfs <- list(df.system,df.cells,df.tissues,df.mixed)
  return(list.dfs)
}

function.gtex.set_variables <- function(df) {
  ### USE: GTEx
  ### STEPS: adding new variables
  df.res <- df
  ### Set "significance" variable
  df.res$FDR.significant <- with(df.res, ifelse( (False.discovery.rate=="<0.01"|False.discovery.rate=="<0.05"), TRUE, FALSE))
  ### Set "group"/tissue variable
  tmp.strsplit <- strsplit(as.character(df.res[,group_variable.gtex]), " - ") # ALTERNATIVELY: #tmp.strsplit <- with(df.res, strsplit(as.character(Name), " - "))
  df.res$group <- sapply(tmp.strsplit, "[[", 1) # or use unlist() | or "lapply(strsplit(XX," - "), function(x) x[1])"
  
  ### ADDITIONAL VARIABLES | not needed
  #df.res <- df.res %>% group_by(group) %>% mutate(group.width=n()) # calculating the number of observations in each group ("group width")
  
  ### ***SORTING*** 
  df.res <- df.res[order(df.res[,group_variable.gtex]),] # IMPORTANT: sorting data frame
  
  return(df.res)
}

function.add.spacers <- function(df) {
  ### USE: UNIVERSIAL
  ### Function to INSERT spacers (dummy rows) in data frame between "groups"
  ### *IMPORTANT*:we are relying the data frame being (correctly) SORTED!!
  df.tissue_enrichment.spaced <- data.frame()
  group.previously <- df$group[1]
  for (i in 1:nrow(df)) {
    ### *OBS*: df MUST be sorted at this point!
    group.now <- df$group[i]
    if ( group.now != group.previously ) { # add "spacer" (blank) rows between "groups"
      df.dummy <- data.frame(Name="dummy",Nominal.P.value=1,False.discovery.rate=">=0.20")
      for (j in 1:n.blank_spacers) {
        # use dplyr::bind_rows() to fill unmatched columns with NA [rbind() will complain]
        df.tissue_enrichment.spaced <- suppressWarnings(dplyr::bind_rows(df.tissue_enrichment.spaced, df.dummy))
      }
      group.previously <- df$group[i]
    }
    df.tissue_enrichment.spaced <- dplyr::bind_rows(df.tissue_enrichment.spaced, df[i,])
    # GTEx --> Name  Nominal.P.value  False.discovery.rate
  }
  return(df.tissue_enrichment.spaced)
}


function.plot.barplot <- function(df.tissue_enrichment, xlabel="MISSING X-LABEL") {
  ### USE: UNIVERSIAL
  ### Function creates a ggplot object [geom_bar()]
  ### INPUT: 
  # df.tissue_enrichment: a sorted "cleaned" data frame
  # xlabel: a character string used as xlabel
  ### OUTPUT: 
  # p: a ggplot object
  
  ### STEPs
  # 1) create breaks and labels for plot
  # 2) *call function* "function.add.spacers()"
  # 3) make plot and adjust theme
  
  
  ##########################################################
  ################## Preparing for PLOTTING ################
  
  ### Constructing labels and break positions for x-axis
  # *OBS* df.tissue_enrichment is used for input and NOT the "spaced" data frame
  df.breaks_and_labels <- df.tissue_enrichment %>% group_by(group) %>% summarize(group.width=n()) # *OBS*: this code relies on the "group" being in the correct order
  df.breaks_and_labels$order.group.numeric <- seq(0,nrow(df.breaks_and_labels)-1) # 0,1,...,3
  df.breaks_and_labels$group.width.cumsum <- with(df.breaks_and_labels, cumsum(group.width)) # cumsum()
  df.breaks_and_labels$group.width.cumsum.shift <- with(df.breaks_and_labels, c(0,group.width.cumsum[-length(group.width.cumsum)])) # *OBS*: "shifting" position by one (pop array). Inserting zero in first position | try diff(x)?
  df.breaks_and_labels$break_position <- with(df.breaks_and_labels, 0.5+(n.blank_spacers*order.group.numeric)+group.width.cumsum.shift+group.width/2) 
  
  ### Adding spacers to data frame | *OBS*: we are relying on CORRECT SORTING of the data frame
  df.tissue_enrichment.spaced <- function.add.spacers(df.tissue_enrichment)
  df.tissue_enrichment.spaced$order.numeric <- 1:nrow(df.tissue_enrichment.spaced) # Add numeric code | maps to x-axis
  #str(df.tissue_enrichment.spaced)
  
  ##########################################################
  ########################## PLOT ##########################
  
  # ===================== CORE PLOT ====================== # 
  p <- ggplot(df.tissue_enrichment.spaced)
  p <- p + geom_bar(aes(x=order.numeric, y=-log10(Nominal.P.value), fill=FDR.significant), stat="identity")
  p <- p + scale_fill_manual(name="", values=c("TRUE"="#F15A22","FALSE"="#606060",guide='legend')) # set colors for FDR significant bars
  p <- p + scale_x_continuous(breaks=df.breaks_and_labels$break_position, labels=df.breaks_and_labels$group) # add x-axis labels
  p <- p + labs(y=expression(-log[10](paste(italic(P)," value")))) # add axis titles
  p <- p + guides(fill=FALSE) # remove legend
  
  # =============== FINE TUNING - theme() =============== # 
  ### Do not out-comment individual lines in the below block of code. It should either be out-commented or un-commented all together
  p.theme <- theme_classic() # default is "theme_gray()" "theme_classic()" theme removes background and more [white background, eliminates background, gridlines, and chart border]. Note that it contains black colored axis lines. Therefore it is used as a template for further modifications
  # OR --> p <- p + theme_bw() + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank()) 
  # See also theme_classic(), theme_minimal() and theme_bw()
  # Try the "theme() %+replace% theme()" operator
  
  ### Setting axis line for y-axis (only) | The following did not work --> theme(axis.line.y=element_line(color="black"))
  p.theme$axis.line.x <- element_blank() # empty list of class "element_blank" and "element" | SAME as theme_bw()$axis.line
  #p.theme$axis.line.x <- p.theme$axis.line # copy list
  #p.theme$axis.line.x$colour <- 'black'
  p.theme$axis.line.y <- p.theme$axis.line # copy list
  p.theme$axis.line.y$colour <- 'black' # OBS: sensitive to spelling: "colour" and NOT "color"
  
  
  ### Setting axis tick marks for y-axis (only)
  p.theme$axis.ticks.x <- element_blank()
  p.theme$axis.ticks.y <- p.theme$axis.ticks
  p.theme$axis.ticks.y$colour <- 'black' # redundant: theme_bw() has this by default
  # traditional approach --> theme(axis.ticks=element_blank())
  
  ### Adjusting length of tick marks and expansion
  p.theme$axis.ticks.length <- unit(0.15, "cm") # adjust distance from axis labels to axis | you may need to play with this a little
  # Note that setting different lengths for tickmarks for x and y does not work | that is, manually setting "axis.ticks.length.x"/"axis.ticks.length.x" does *NOT* work.
  # traditional approach --> theme(axis.ticks.length=unit(0, "cm"))
  p <- p + scale_y_continuous(expand = c(0, 0)) # removes expansion of y-axis. Now the y-axis starts at zero!
  
  ### Adjust x-axis labels (size and rotation)
  p.theme <- p.theme + theme(axis.text.x=element_text(angle=35, hjust=1, size=rel(0.8))) # size=rel(1.15) consider using %+replace% operator
  
  # ===================== Combining theme + plot ====================== # 
  p <- p + p.theme # saving "p.theme" into plot
  
  # ===================== input_type specific ====================== # 
  p <- p + labs(x=xlabel)
  
  # ===================== return value ====================== # 
  return(p) # return ggplot object
}

function.plot.save <- function(p, filename.prefix="MISSING-FILENAME-PREFIX", filename.suffix="MISSING-FILENAME-SUFFIX") {
  ### USE: UNIVERSIAL
  ### Function exports a ggplot object to pdf file
  ### INPUT: 
  # p: a ggplot object
  ### OUTPUT
  # NONE
  # ======================== Save plot =============================== # 
  #p.filename <- sprintf("tissue_plot_%s.pdf", filename.suffix)
  p.filename <- sprintf("tissue_plot_%s_%s.pdf", filename.prefix, filename.suffix)
  suppressWarnings(ggsave(file=p.filename, width=8, height=4, units="in", dpi=300))
  print(sprintf("Saved plot %s", p.filename))
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # TAKEN FROM http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


######################################################
################## READ ARGUMENTS ####################
### 
### R Base cmd line arguments
# Accept command line argument | trailingOnly=TRUE means that only the users arguments are returned
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  filename.prefix = "no_prefix"
} else if (length(args)==2) {
  filename.prefix = args[2]
} else {
  cat("### Arguments received ###\n")
  for (i in seq_along(args)) {cat(sprintf("Argument #%s: %s\n", i, args[i]))}
  stop(sprintf("Received wrong number of arguments. Please only specify 1 or 2 arguments."))
}

### Save arguments ###
file.tissue_enrichment <- args[1]

### Validation of arguments ###
if( file.access(file.tissue_enrichment) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file.tissue_enrichment))
}


#########################################################################
################################## MAIN #################################
##################### READ DATA AND CALL FUNCTIONS ######################
#########################################################################

#setwd("/Users/pascaltimshel/Dropbox/0_Work/DEPICT/DEPICT_scripts_PT") # *TEMPORARY* - only for development
#file.tissue_enrichment <- "data/ea_tissueenrichment_gtex.txt" # *TEMPORARY* - only for development
#file.tissue_enrichment <- "data/ea_tissueenrichment_genenetwork.txt" # *TEMPORARY* - only for development

df.input.raw <- read.table(file.tissue_enrichment, header=T, sep="\t")

### Auto detection of tissue enrichment format: genenetwork vs gtex
if ( all(colnames(df.input.raw[,1:6]) == colnames.genenetwork) ) { # or use read.table(..., check.names=FALSE)
  ### Set flag
  print("Detected Gene Network file format")
  flag.input_type <- "genenetwork"
  
  df.tissue_enrichment <- function.input_procces(df=df.input.raw, cols2keep=colnames.genenetwork)
  ### Processing
  list.dfs <- function.genenetwork.set_variables_and_split(df.tissue_enrichment)
  ### Unpacking list
  df.system <- list.dfs[[1]]
  df.cells <- list.dfs[[2]]
  df.tissues <- list.dfs[[3]]
  df.mixed <- list.dfs[[4]]
  
  
  # Plot and export
  p.system <- function.plot.barplot(df.system, xlabel="Physiological system")
  function.plot.save(p.system, filename.prefix=filename.prefix, filename.suffix="genenetwork_system")
  
} else if ( all(colnames(df.input.raw[,1:3]) ==colnames.gtex) ) {
  ### Set flag
  print("Detected GTEx file format")
  flag.input_type <- "gtex"
  
  df.tissue_enrichment <- function.input_procces(df=df.input.raw, cols2keep=colnames.gtex)
  ### Updating df
  df.tissue_enrichment <- function.gtex.set_variables(df.tissue_enrichment)
  
  # Plot and export
  p.gtex <- function.plot.barplot(df.tissue_enrichment, xlabel="Tissues")
  function.plot.save(p.gtex, filename.prefix=filename.prefix, filename.suffix="gtex")
  
} else {
  stop("ERROR. Could not recognize tissue enrichment file format. Please check the file format.")
}


##########################################################################################################
##########################################################################################################

####################################################################
########################## Extended plots ##########################
if ( as.logical(flag.genenetwork.plot_extended_mesh_categories) & (flag.input_type == "genenetwork") ) {
  print("Plotting exteded files for Gene Networks")
  
  # CELLS | Plot and export
  p.cells <- function.plot.barplot(df.cells, xlabel="Cells")
  function.plot.save(p.cells, filename.prefix=filename.prefix, filename.suffix="genenetwork_cells")
  
  # TISSUES | Plot and export
  p.tissues <- function.plot.barplot(df.tissues, xlabel="Tissues")
  function.plot.save(p.tissues, filename.prefix=filename.prefix, filename.suffix="genenetwork_tissues")
  
  # MULTIPLOT | System, Cells, Tissues
  m.plot.layout <- matrix(c(1,1,2,3), nrow=2, byrow=TRUE)
  m.plot.filename <- sprintf("tissue_plot_%s_genenetwork_multiplot.pdf", filename.prefix)
  pdf(file=m.plot.filename, width=12, height=5)
  suppressWarnings(multiplot(p.system, p.cells, p.tissues, layout=m.plot.layout))
  print(sprintf("Saved multiplot %s", m.plot.filename))
  dummy <- dev.off() # graphics.off() 
  
}



