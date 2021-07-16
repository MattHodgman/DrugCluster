# runFastPG.r runs the FastPG program (https://github.com/sararselitsky/FastPG), an R implementation of the Phenograph method, to cluster cells by the markers
# included in the input csv.
#
# Arguments:
#     1. the cleaned data input csv
#     2. the local neighborhood size (k)
#     3. the number of cpus to use in the k nearest neighbors part of clustering
#     4. output directory
#     5. output file name for cell/cluster assignment
#     6. output file name for cluster mean feature values
#     7. flag to include method name as a column
#
# Output: 
#     cells.csv - which contains the cell ID and cluster ID
#     clusters.csv - which contains the mean expression values for each marker, for each cluster


# get data and cluster it
args <- commandArgs(trailingOnly=TRUE) # required command line arguments order: {cleaned data csv} {k} {num_threads} {output dir}
data <- as.matrix(read.csv(file=args[1])) # load cleaned data into matrix
clusters <- FastPG::fastCluster(data=data, k=as.integer(args[2]), num_threads=as.integer(args[3])) # compute clusters
Cluster <- clusters$communities # get all cell community assignations (these are in the same order as cells in data)
data <- cbind(Cluster, data) # add community assignation to data

# make cells.csv
cells <- (data[,c('CellID','Cluster')]) # get just cell IDs and community assignations for export
if (as.logical(args[7])) { # inlcude method column
    Method <- rep(c('FastPG'),nrow(cells))
    cells <- cbind(cells, Method)
}
write.table(cells,file=paste(args[4], args[5], sep='/'),row.names=FALSE,quote=FALSE,sep=',') # write data to csv

# make clusters.csv
clusterData <- aggregate(subset(data, select=-c(CellID)), list(data[,'Cluster']), mean) # group feature/expression data by cluster and find mean expression for each cluster, remove CellID column
clusterData <- subset(clusterData, select=-c(Group.1)) # remove group number column because is identical to community assignation number
if (as.logical(args[7])) { # inlcude method column
    Method <- rep(c('FastPG'),nrow(clusterData))
    clusterData <- cbind(clusterData, Method)
}
write.table(clusterData,file=paste(args[4], args[6], sep='/'),row.names=FALSE,quote=FALSE,sep=',') # write data to csv

cat(clusters$modularity) # output modularity