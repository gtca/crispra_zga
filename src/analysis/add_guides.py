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
from anndata.base import AnnData
import scanpy as sc

os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"


def write_log(statlogfile: str, stats):
	with open(statlogfile, 'a') as slf:
		yaml.dump(stats, slf, default_flow_style=False)


def add_guides(adata: AnnData, guide_info: pd.DataFrame, barcodes_with_guides: pd.DataFrame,
			   min_sgrna_frac: float = 0.9, max_sgrna_frac_se: float = 0.1,
			   log: types.LambdaType = lambda x: None
			  ):
	
	log({"n_cells_total": adata.obs.shape[0]})
	print(f"Total cell count before assigning guides: {adata.obs.shape[0]}")

	# Add transfection_lane field, e.g. t1_l1
	barcodes_with_guides["transfection_lane"] = ['t' + '_l'.join(i.split('-')[1:]) for i in barcodes_with_guides.barcode.values]
	barcodes_with_guides.transfection_lane = pd.Categorical(barcodes_with_guides.transfection_lane)

	# Select barcodes assigned with certain confidence
	barcodes_with_guides = barcodes_with_guides[(barcodes_with_guides["sgrna_frac"] > min_sgrna_frac) & (barcodes_with_guides["sgrna_frac_se"] < max_sgrna_frac_se)]

	# Add guides to adata observations
	print("Assigning guides to cells...")
	adata.obs["sgrna"] = adata.obs.join(barcodes_with_guides.loc[:,["barcode", "sgrna"]].set_index("barcode"))["sgrna"].values

	# Filter cells
	print("Filtering cells...")
	adata = adata[adata.obs.sgrna.notna()]
	adata = adata[adata.obs.sgrna != "EMPTY"]

	log({"n_cells_total": adata.obs.shape[0]})
	print(f"Total cell count after assigning guides: {adata.obs.shape[0]}")

	# Join guide metadata
	print("Adding guide information...")
	joined_obs = adata.obs.join(guide_info.set_index("guide_sequence"), how="left", on="sgrna")
	for n in ["sgrna_id", "read_1", "read_2", "guide_length", "transcript_id", "gene_name"]:
		adata.obs[n] = joined_obs[n]
	del joined_obs

	return adata



if __name__ == "__main__":

	# Parse input arguments

	input_file        = snakemake.input[0]
	guides_info_file  = snakemake.input[1]
	guides_cells_file = snakemake.input[2]
	output_file       = snakemake.output[0]
	options, statlog  = snakemake.config["analysis"]["add_guides"]["scanpy"]["options"], snakemake.config["analysis"]["add_guides"]["scanpy"]["statlog"]

	# Flush statlog file content if any
	open(statlog, 'w').close()
	
	log_func = lambda stats: write_log(statlog, stats) if statlog is not None else None
	
	# Load files
	print("Reading files...")
	adata = sc.read(input_file)
	guide_info = pd.read_csv(guides_info_file)
	barcodes_with_guides = pd.read_csv(guides_cells_file)

	adata = add_guides(adata, guide_info, barcodes_with_guides, log=log_func, **options)

	print("Saving files...")
	adata.write(output_file)

	print(f"The AnnData object with guides assigned and their meta-information is saved to {output_file}")


