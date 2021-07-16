# Command line arguments:
# 1. clean_data.fcs
# 2. number of meta clusters
# 3. flag to include method name as column
# 4. output directory
# 5. output file name for cell/cluster assignment
# 6. output file name for cluster mean feature values


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

# make cells.csv
cells <- data_raw[,c('CellID','Cluster')]
if (as.logical(args[3])) { # inlcude method column
    Method <- rep(c('FastPG'),nrow(cells))
    cells <- cbind(cells, Method)
}
write.table(cells,file=paste(args[4], args[5], sep='/'), row.names=FALSE, quote=FALSE, sep=',') # write data to csv

# make clusters.csv
clusterData <- aggregate(subset(data_raw, select=-c(CellID)), list(data_raw[,'Cluster']), mean) # group feature/expression data by cluster and find mean expression for each cluster, remove CellID column
clusterData <- subset(clusterData, select=-c(Group.1)) # remove group number column because is identical to community assignation number
if (as.logical(args[3])) { # inlcude method column
    Method <- rep(c('FastPG'),nrow(clusterData))
    clusterData <- cbind(clusterData, Method)
}
write.table(clusterData,file=paste(args[4], args[6], sep='/'),row.names=FALSE,quote=FALSE,sep=',') # write data to csv