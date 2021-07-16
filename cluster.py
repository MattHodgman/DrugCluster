import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import subprocess
import os


'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Cluster drugs using the semantic distances between them.')
    parser.add_argument('-i', '--input', help="Input CSV of mcmicro marker expression data for cells", type=str, required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved', type=str, required=False)
    parser.add_argument('-k', '--neighbors', help='the number of nearest neighbors to use when clustering. The default is 30.', default=30, type=int, required=False)
    parser.add_argument('-c', '--method', help='Include a column with the method name in the output files.', action="store_true", required=False)
    parser.add_argument('-a', '--algorithm', help='Which clustering algorithm to use. Options are FastPG, Scanpy, or FlowSOM.', type=str, required=True)
    parser.add_argument('-n', '--num-metaclusters', help='number of clusters for meta-clustering.', type=int, required=False, default=10)
    args = parser.parse_args()
    return args


'''
Get the path to the directory where this script is located and return it.
'''
def get_path():
    full_path = os.path.realpath(__file__)
    path_list = full_path.split('/')[:-1]
    s = '/'
    return s.join(path_list)


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1] # get filename from end of input path
    dataName = fileName[:fileName.rfind('.')] # get data name by removing extension from file name
    return dataName


'''
Write clusters_file from leidenCluster() adata
'''
def writeClusters(adata):
    cells = pd.DataFrame(adata.obs[DRUG].astype(int)) # extract cell IDs to dataframe
    cells[CLUSTER] = adata.obs[LEIDEN] # extract and add cluster assignments to cells dataframe

    # add in method column if requested
    if args.method:
        cells[METHOD] = SCANPY

    cells.to_csv(f'{output}/{clusters_file}', index=False)


'''
Cluster data using the Leiden algorithm via scanpy
'''
def leidenCluster():

    sc.settings.verbosity = 3 # print out information
    adata_init = sc.read(f'{output}/{data_file}', cache=True) # load in clean data

    # move CellID info into .obs
    # this assumes that 'CELL_ID' is the first column in the csv
    adata_init.obs[DRUG] = adata_init.X[:,0]
    adata = ad.AnnData(np.delete(adata_init.X, 0, 1), obs=adata_init.obs, var=adata_init.var.drop([DRUG]))

    # compute neighbors and cluster
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs=10) # compute neighbors, using the number of neighbors provided in the command line. Default is 30.
    sc.tl.leiden(adata, key_added = LEIDEN, resolution=1.0) # run leidan clustering. default resolution is 1.0

    # write cluster assignments 'CLUSTERS_FILE'
    writeClusters(adata)


'''
Run an R script that converts the clean data csv into fcs format
'''
def convertToFSC():
    if args.verbose:
        print('Converting CSV to FCS...')

    path = get_path() # get the path where the r script is located

    # Build subprocess command
    r_script = ['Rscript', f'{path}/CSVtoFCS.r'] # use CSVtoFCS.r script

    # pass args
    r_args = [output, data_file]

    # Build subprocess command
    command = r_script + r_args

    # run it
    subprocess.run(command, universal_newlines=True)

    if args.verbose:
        print('Done.')


'''
Run an R script that runs flowSOM
'''
def runFlowSOM():
    if args.verbose:
        print('Running flowSOM...')

    path = get_path() # get the path where the r script is located

    r_script = ['Rscript', f'{path}/runFlowSOM.r'] # use FastPG.r script
    # pass input data file, k value, number of cpus to use for the k nearest neighbors part of clustering, output dir, cells file name, clusters file name
    r_args = [f'{output}/{data_fcs_file}', str(args.num_metaclusters), str(args.method), output, clusters_file]

    # Build subprocess command
    command = r_script + r_args

    # run it
    subprocess.run(command, universal_newlines=True)

    if args.verbose:
        print('Done.')


'''
Run an R script that runs FastPG. Scriptception.
'''
def runFastPG():
    if args.verbose:
        print('Running R script...')
    
    path = get_path() # get the path where the r script is located

    r_script = ['Rscript', f'{path}/runFastPG.r'] # use FastPG.r script
    # pass input data file, k value, number of cpus to use for the k nearest neighbors part of clustering, output dir, cells file name, clusters file name
    r_args = [f'{output}/{data_file}', str(args.neighbors), str(args.num_threads), output, clusters_file, str(args.method)]

    # Build subprocess command
    command = r_script + r_args

    # run R script and get modularity from stdout 
    modularity = subprocess.check_output(command, universal_newlines=True)

    if args.verbose:
        print(f'Modularity: {modularity}')
        print('Done.')


'''
Main.
'''
if __name__ == '__main__':
    args = parseArgs() # parse arguments
    
    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    # constants
    DRUG = 'Drug' # column name holding drug name
    CLUSTER = 'Cluster' # column name holding cluster number
    LEIDEN = 'leiden' # obs name for cluster assignment
    METHOD = 'Method' # name of column containing the method for clustering
    SCANPY = 'Scanpy' # the name of this method
    FASTPG = 'FastPG'
    FLOWSOM = 'FlowSOM'
    
    # file names
    data_prefix = getDataName(args.input) # get the name of the input data file to add as a prefix to the output file names
    data_file = f'{data_prefix}.csv'
    data_fcs_file = f'{data_prefix}.fcs' # name of output data CSV file
    clusters_file = f'{data_prefix}-clusters.csv' # name of output CSV file that contains the cluster assignment for each drug
    
    # Select clustering algorithm
    algorithm = args.algorithm
    if algorithm == FASTPG:
        runFastPG()
    elif algorithm == SCANPY:
        leidenCluster()
    elif algorithm == FLOWSOM:
        convertToFSC()
        runFlowSOM()