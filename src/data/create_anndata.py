# 
# Load 10X scRNA-seq data for all the runs
# 

import sys
import os
import glob
import types
from typing import List
import yaml

import numpy as np
import pandas as pd
import scanpy.api as sc


os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"


def write_log(statlogfile: str, stats):
	with open(statlogfile, 'a') as slf:
		yaml.dump(stats, slf, default_flow_style=False)


def create_anndata(input_dir: str, barcodes_file: str,
				   exclude_samples: List[str] = [],
				   log: types.LambdaType = lambda x: None 
				   ):

	# Load CellRanger output
	ad_list = []
	run_ids = []

	for run_id in os.listdir(input_dir):
	    if run_id[0] != 'S':  # all the run IDs start with S
	        continue
	    print(f"Loading {run_id}...     ", end='', flush=True)
	    adata = sc.read_10x_h5(input_dir + '/{}/outs/filtered_gene_bc_matrices_h5.h5'.format(run_id), genome='mm10')
	    adata.obs["sample_id"] = run_id
	    ad_list.append(adata)
	    run_ids.append(run_id)

	# Construct AnnData object
	adata = ad_list[0].concatenate(*ad_list[1:])
	adata.var_names_make_unique()

	del ad_list

	adata.var = adata.var[['gene_ids-0']]
	adata.var.columns = ['gene_ids']

	print("Constructed AnnData object with dimensions {}".format(str(adata.shape)))
	log({"samples_count_raw": len(run_ids)})
	log({"cells_total_count_raw": adata.shape[0]})
	log({"genes_total_count_raw": adata.shape[1]})


	# Read and add metadata
	print("Adding metadata from {}...".format(barcodes_file))

	meta = pd.read_table(barcodes_file, sep=",")

	meta.index = meta["barcode"]

	for new_obs in ["transfection_lane", "transfection", "pool", "index_seq", "contamination"]:
		adata.obs[new_obs] = list(meta.loc[adata.obs["sample_id"].values][new_obs])

	# Standardise cell names
	#   from: barcode-1-batch, where batch is 0, 1, etc. in the order of loading the sample files
	#     to: barcode-transfection-lane
	adata.obs.index = ['-'.join((a, b)) for a, b in zip([s.split('-')[0] for s in adata.obs.index], [s.replace('_', '-') for s in adata.obs.transfection_lane.values])]

	# Save AnnData file

	log({"samples_count": len(adata.obs.sample_id.unique())})
	log({"cells_total_count": adata.shape[0]})


	return adata

	


if __name__ == "__main__":

	# Parse input arguments

	input_dir        = snakemake.input[0]
	barcodes_file    = snakemake.input[1]
	output_file      = snakemake.output[0]
	options, statlog = snakemake.config["data"]["transcriptome"]["scanpy"]["options"], snakemake.config["data"]["transcriptome"]["scanpy"]["statlog"]

	# Flush statlog file content if any
	open(statlog, 'w').close()

	log_func = lambda stats: write_log(statlog, stats) if statlog is not None else None


	adata = create_anndata(input_dir, barcodes_file, log=log_func, **options)
	print("Saving expression to {}...".format(output_file))
	adata.write(output_file)

