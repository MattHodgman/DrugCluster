# Command line arguments:
# 1. clean_data.fcs
# 2. number of meta clusters
# 3. flag to include method name as column
# 4. output directory
# 5. output file name for drug/cluster assignment


# paramenters we could edit
# number clusters
# iterations
# transformations?? this is some sort of normalization... which I dont think we want to do (already done?)
# compensation?


# install packages
if(!require('flowCore')) {install.packages('flowCore')}
if(!require('FlowSOM')) {install.packages('FlowSOM')}

# load libraries
library("flowCore")
library("FlowSOM")


# get arguments
args <- commandArgs(trailingOnly=TRUE)

# read in data
data <- read.FCS(args[1])

# get num cols
num_cols <- strtoi(keyword(data, '$PAR')[1])

# run FlowSOM, not sure what to put for maxMeta and not sure if it matters
# from the paper (FlowSOM: Using Self-Organizing Maps forVisualization and Interpretation ofCytometry Data): 
# "To choose the number of meta-clusters, onecan either use prior knowledge about the number of expectedcell types, or one can use the so called “elbow”-criterion."
fSOM <- FlowSOM(data, colsToUse=c(1:num_cols), nClus=as.integer(args[2]))

# get cluster assignments
Cluster <- GetClusters(fSOM)

# get raw input data and add clusters
data_raw <- cbind(Cluster, exprs(data))

# make clusters.csv
cells <- data_raw[,c('Drug','Cluster')]
if (as.logical(args[3])) { # inlcude method column
    Method <- rep(c('FastPG'),nrow(cells))
    cells <- cbind(cells, Method)
}
write.table(cells,file=paste(args[4], args[5], sep='/'), row.names=FALSE, quote=FALSE, sep=',') # write data to csv