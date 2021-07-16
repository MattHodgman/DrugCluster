import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Cluster cell types using mcmicro marker expression data.')
    parser.add_argument('-i', '--input', help="Input CSV of mcmicro marker expression data for cells", type=str, required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved', type=str, required=False)
    parser.add_argument('-k', '--neighbors', help='the number of nearest neighbors to use when clustering. The default is 30.', default=30, type=int, required=False)
    parser.add_argument('-c', '--method', help='Include a column with the method name in the output files.', action="store_true", required=False)
    parser.add_argument('-a', '--algorithm', help='Which clustering algorithm to use. Options are FastPG, Scanpy, or FlowSOM.', type=str, required=True)
    args = parser.parse_args()
    return args


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1] # get filename from end of input path
    dataName = fileName[:fileName.rfind('.')] # get data name by removing extension from file name
    return dataName


'''
Write CELLS_FILE from leidenCluster() adata
'''
def writeCells(adata):
    cells = pd.DataFrame(adata.obs[CELL_ID].astype(int)) # extract cell IDs to dataframe
    cells[CLUSTER] = adata.obs[LEIDEN] # extract and add cluster assignments to cells dataframe

    # add in method column if requested
    if args.method:
        cells[METHOD] = SCANPY

    cells.to_csv(f'{output}/{cells_file}', index=False)


'''
Write CLUSTERS_FILE from leidenCluster() adata
'''
def writeClusters(adata):
    clusters = pd.DataFrame(columns=adata.var_names, index=adata.obs[LEIDEN].cat.categories)   
    clusters.index.name = CLUSTER # name indices as cluster column
    for cluster in adata.obs.leiden.cat.categories: # this assumes that LEIDEN = 'leiden' if the name is changed, replace it for 'leiden' in this line
        clusters.loc[cluster] = adata[adata.obs[LEIDEN].isin([cluster]),:].X.mean(0)
    
    # add in method column if requested
    if args.method:
        clusters[METHOD] = SCANPY

    clusters.to_csv(f'{output}/{clusters_file}')


'''
Cluster data using the Leiden algorithm via scanpy
'''
def leidenCluster():

    sc.settings.verbosity = 3 # print out information
    adata_init = sc.read(f'{output}/{data_file}', cache=True) # load in clean data

    # move CellID info into .obs
    # this assumes that 'CELL_ID' is the first column in the csv
    adata_init.obs[CELL_ID] = adata_init.X[:,0]
    adata = ad.AnnData(np.delete(adata_init.X, 0, 1), obs=adata_init.obs, var=adata_init.var.drop([CELL_ID]))

    # compute neighbors and cluster
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs=None) # compute neighbors, using the number of neighbors provided in the command line. Default is 30.
    sc.tl.leiden(adata, key_added = LEIDEN, resolution=1.0) # run leidan clustering. default resolution is 1.0

    # write cell/cluster information to 'CELLS_FILE'
    writeCells(adata)

    # write cluster mean feature expression to 'CLUSTERS_FILE'
    writeClusters(adata)


'''
Main.
'''
if __name__ == '__main__':
    args = parseArgs() # parse arguments

    # get input data file
    data_file = args.input

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    # constants
    CELL_ID = 'CellID' # column name holding cell IDs
    CLUSTER = 'Cluster' # column name holding cluster number
    LEIDEN = 'leiden' # obs name for cluster assignment
    METHOD = 'Method' # name of column containing the method for clustering
    SCANPY = 'Scanpy' # the name of this method
    
    # output file names
    data_prefix = getDataName(args.input) # get the name of the input data file to add as a prefix to the output file names
    clusters_file = f'{data_prefix}-clusters.csv' # name of output CSV file that contains the mean expression of each feaute, for each cluster
    cells_file = f'{data_prefix}-cells.csv' # name of output CSV file that contains each cell ID and it's cluster assignation
    
    # cluster using scanpy implementation of Leiden algorithm
    leidenCluster()