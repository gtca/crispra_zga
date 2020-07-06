#!/usr/local/bin/python

# 
# Perform PCA, tSNE, PAGA for the preprocess AnnData
# 

import os
import sys
import argparse
import yaml
import types

import numpy as np
import pandas as pd
import scanpy.api as sc

import matplotlib.pyplot as plt
import seaborn as sns

os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"


def write_log(statlogfile: str, stats):
	with open(statlogfile, 'a') as slf:
		yaml.dump(stats, slf, default_flow_style=False)


def analyse_with_scanpy(input_file: str, output_file: str,
						lovain_resolution: float,
						rank_genes_groups_method: str,
						log: types.LambdaType = lambda x: None 
						):

	# Load AnnData object
	adata = sc.read(input_file)

	sc.tl.pca(adata)

	sc.pp.neighbors(adata, n_neighbors=15)

	sc.tl.tsne(adata, random_state=2, n_pcs=40)
	sc.tl.umap(adata, random_state=1, spread=3.0)
	sc.tl.louvain(adata, resolution=1.0)
	sc.tl.rank_genes_groups(adata, 'louvain', method=rank_genes_groups_method)
	sc.tl.paga(adata)

	log({"total_clusters": len(adata.obs.louvain.unique())})

	adata.write(output_file)




if __name__ == "__main__":

	# Parse input arguments

	input_file       = snakemake.input[0]
	output_file      = snakemake.output[0]
	options, statlog = snakemake.config["analysis"]["transcriptome"]["scanpy"]["options"], snakemake.config["analysis"]["transcriptome"]["scanpy"]["statlog"]

	# Flush statlog file content if any
	open(statlog, 'w').close()
	
	log_func = lambda stats: write_log(statlog, stats) if statlog is not None else None

	analyse_with_scanpy(input_file, output_file, log=log_func, **options)

