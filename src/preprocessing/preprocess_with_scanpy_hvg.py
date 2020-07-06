#!/usr/local/bin/python

# 
# Preprocess 10X scRNA-seq data for all the runs
# 

import os
import sys
import argparse
import yaml
import types

import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns

os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"

from src.log import *

def filter_normalise(adata: AnnData,
					 min_counts_per_cell:  int,
					 max_counts_per_cell:  int,
				     min_genes_per_cell:   int,
				     max_genes_per_cell:   int,
				     min_cells_per_gene:   int,
				     max_mito_threshold:   float,
				     counts_per_cell_norm: float,
				     log: types.LambdaType = lambda x: None 
				     ):

	genes_per_cell = np.count_nonzero(sp.sparse.csr_matrix.todense(adata.X), axis=1)
	log({"genes_per_cell_raw": '[' + ','.join([str(e) for e in genes_per_cell.flatten().tolist()[0]]) + ']'})

	cells_per_gene = np.count_nonzero(sp.sparse.csr_matrix.todense(adata.X), axis=0)
	log({"cells_per_gene_raw": '[' + ','.join([str(e) for e in cells_per_gene.flatten().tolist()[0]]) + ']'})

	log({"cells_total_count_raw": adata.shape[0]})
	log({"genes_total_count_raw": adata.shape[1]})


	sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
	sc.pp.filter_cells(adata, max_counts=max_counts_per_cell)


	sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
	sc.pp.filter_cells(adata, max_genes=max_genes_per_cell)

	sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)

	mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
	adata.obs['n_counts'] = adata.X.sum(axis=1)

	# adata = adata[adata.obs['n_genes'] < 2000, :]
	adata = adata[adata.obs['percent_mito'] < max_mito_threshold, :]

	# DEPRECATED
	# NOTICE: saved before normalizing as in the tutorial
	# adata.raw = sc.pp.log1p(adata, copy=True)
	
	sc.pp.normalize_per_cell(adata, counts_per_cell_after=float(counts_per_cell_norm))

	# NOTICE: saved after normalizing
	adata.raw = sc.pp.log1p(adata, copy=True)

	log({"cells_total_count_filtered": adata.shape[0]})
	log({"genes_total_count_filtered": adata.shape[1]})

	return adata


def mnn_correct(adata,
	            log: types.LambdaType = lambda x: None 
			    ):
	# See the API here: https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.mnn_correct.html

	# Split adata object to multiple objects
	adatas = [adata[np.where([i == sample_id for i in adata.obs.sample_id.values])[0],:] for sample_id in adata.obs.sample_id.values.unique()]

	# Perform correction and concatenate objects
	print(f"Performing MNN correction for {len(adatas)} datasets...")
	mnn_result = sc.pp.mnn_correct(*adatas)
	return mnn_result[0]



def log_scale(adata,
	          regress_counts_and_mito: bool = False,
			  log: types.LambdaType = lambda x: None 
			  ):

	# Logarithmise the data
	adata = sc.pp.log1p(adata, copy=True)

	if regress_counts_and_mito:
		# Regress out effects of total counts
		sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

	# Scale to unit variance
	sc.pp.scale(adata, max_value=10)
	return adata



def hvg_log_scale(adata,
	              regress_counts_and_mito: bool = False,
				  log: types.LambdaType = lambda x: None 
				  ):

	# Select highly variable genes
	filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.01, max_mean=5, min_disp=0.5)

	# Logarithmise the data
	adata = sc.pp.log1p(adata, copy=True)

	# Subset highly variable genes
	adata = adata[:, filter_result.gene_subset]

	log({"cells_total_count_hvg": adata.shape[0]})
	log({"genes_total_count_hvg": int(adata.shape[1])})  # outputs object numpy scalar object without explicit conversion
	
	if regress_counts_and_mito:
		# Regress out effects of total counts
		sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

	# Scale to unit variance
	sc.pp.scale(adata, max_value=10)
	return adata


def preprocess_with_scanpy(input_file: str, 
						   output_file_hvg: str,
						   min_counts_per_cell:   int,
						   max_counts_per_cell:   int,
						   min_genes_per_cell:   int,
						   max_genes_per_cell:   int,
						   min_cells_per_gene:   int,
						   max_mito_threshold:   float,
						   counts_per_cell_norm: float,
						   log: types.LambdaType = lambda x: None 
						   ):
	# Load AnnData object
	adata = sc.read(input_file)

	# Filter cells and normalise the expression values
	adata = filter_normalise(adata, min_counts_per_cell, max_counts_per_cell, min_genes_per_cell, max_genes_per_cell, min_cells_per_gene, max_mito_threshold, counts_per_cell_norm, log=log)
	
	# Select only highly variable genes
	# Then log-transform the counts and scale variance
	adata = hvg_log_scale(adata, regress_counts_and_mito=True, log=log)

	print(f"Saving AnnData object to {output_file_hvg}...")
	adata.write(output_file_hvg)

	# DEPRECATED
	# Perform MNN correction
	# adata = mnn_correct(adata, log=log)
	# print(f"Saving AnnData object to {output_file_mnn}...")
	# adata.write(output_file_mnn)


if __name__ == "__main__":

	# Parse input arguments

	input_file          = snakemake.input[0]
	output_file_hvg     = snakemake.output["hvg"]
	options, statlog    = snakemake.config["preprocessing"]["transcriptome"]["scanpy"]["options"], snakemake.config["preprocessing"]["transcriptome"]["scanpy"]["statlog"]

	# Flush statlog file content if any
	open(statlog, 'w').close()

	log_func = lambda stats: write_log(statlog, stats) if statlog is not None else None

	preprocess_with_scanpy(input_file, 
						   output_file_hvg, 
						   log=log_func, 
						   **options)


