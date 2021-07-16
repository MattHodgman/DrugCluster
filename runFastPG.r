# runFastPG.r runs the FastPG program (https://github.com/sararselitsky/FastPG), an R implementation of the Phenograph method, to cluster cells by the markers
# included in the input csv.
#
# Arguments:
#     1. the data input csv
#     2. the local neighborhood size (k)
#     3. the number of cpus to use in the k nearest neighbors part of clustering
#     4. output directory
#     5. output file name for drug/cluster assignment
#     6. flag to include method name as a column
#
# Output:
#     clusters.csv - which contains drug cluster assignments


# get data and cluster it
args <- commandArgs(trailingOnly=TRUE) # required command line arguments order: {data csv} {k} {num_threads} {output dir}
data <- read.csv(file=args[1])) # read data
Drug <- data$Drug
data <- subset(data, select = -c(Drug))
data <- as.matrix(data)
rownames(data) <- Drug
clusters <- FastPG::fastCluster(data=data, k=as.integer(args[2]), num_threads=as.integer(args[3])) # compute clusters
Cluster <- clusters$communities # get all drug community assignations (these are in the same order as drug in data)

# make -clusters.csv
drugs <- data.frame(Drug,Cluster)
if (as.logical(args[6])) { # inlcude method column
    Method <- rep(c('FastPG'),nrow(drugs))
    drugs <- cbind(drugs, Method)
}
write.table(drugs,file=paste(args[4], args[5], sep='/'),row.names=FALSE,quote=FALSE,sep=',') # write data to csv

cat(clusters$modularity) # output modularity