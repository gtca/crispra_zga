#!/usr/local/bin/python

# 
# Assign guides to cell barcodes from the amplicon readout
#
# NOTE: please note vector-specific sequences in parse_reads_with_guides()
# 

import os
import sys
import argparse
import yaml
import types
from typing import List

import numpy as np
import pandas as pd

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import Counter
import gzip
import Levenshtein
from itertools import islice

os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"

from src.log import *


# Calculate the standard error for the fraction estimate
def calculate_frac_se(s_plus, n_total):
	rate = (s_plus + 1.0) / (n_total + 2)
	se = np.sqrt(rate * (1 - rate) / n_total)
	return se


def get_guides_df(guides_file: str, log: types.LambdaType = lambda x: None) -> pd.DataFrame:
	guides_df = pd.read_csv(guides_file)
	guides = guides_df.guide_sequence.values
	guide_names = guides_df.sgrna_id.values
	n_guides = len(guides)
	
	log({"n_guides": n_guides})
	print(f"There are {n_guides} guides.")
	
	# Some guides have alternative length
	print(np.unique(guides_df.guide_length.values, return_counts = True))
	print(guides_df.loc[guides_df.guide_length == 23,:].guide_sequence.values)

	guides_df["seq"] = pd.Series([s[:20] for s in guides_df.guide_sequence])
	
	return guides_df


def parse_reads_with_guides(r3_file: str, guides_df: pd.DataFrame, correct_fragments: bool = True, log: types.LambdaType = lambda x: None):
	lane_name = os.path.basename(r3_file).replace("_R3.fastq.gz", "")
	guides = guides_df.seq.values
	n_guides = len(guides)

	titles = []
	fragms = []
	reads  = []

	print(f"Reading {r3_file}...")
	with gzip.open(r3_file, 'rt') as f:
		for title, seq, qual in FastqGeneralIterator(f):
			titles.append(title)
			fragms.append(seq[23:43])  # expected position of a guide RNA in a read
			reads.append(seq)
	log({"n_reads": len(fragms)})
	print(f"Read {len(fragms)} sequences for {lane_name}")

	# Count sgRNA fragments from reads
	c = Counter(fragms)
	c_common = c.most_common()
	
	# Calculate how many sgRNAs can be detected in the readout
	n_guides_detected = np.sum([j in [i[0] for i in c_common] for j in guides])
	log({"n_guides_detected": int(n_guides_detected)})
	print(f"Detected {n_guides_detected} guide RNAs (out of {n_guides} guides) in the amplicon library for {lane_name}")

	n_reads_not_assigned = np.sum([e[1] for e in c.most_common() if e[0] not in guides])
	n_reads_total = np.sum([e[1] for e in c_common])
	log({"n_reads_total": int(n_reads_total)})
	log({"n_reads_not_assigned": int(n_reads_not_assigned)})
	print(f"{n_reads_not_assigned} of {n_reads_total} reads ({n_reads_not_assigned * 100.0 / n_reads_total:.2f}%) have no guide assigned yet")


	# Correct sequencing mistaked in the guide RNAs fragments
	if correct_fragments:

		gs = guides_df.guide_sequence.values  # full-length guides (not clipped to 20nt)
		sgrnas = []

		# See vector sequences in the misc directory
		up  = 'CTTGTGGAAAGGACGAAACACCG'
		dn  = 'GTTTAAGAGCTAGGCCAACATGA'

		print("Correcting sequencing mistakes in guide RNA fragments for better assignment...")
		for i, e in enumerate(fragms):
			if i % 10000 == 0:
				print("Done {:.2f}%\r".format(i * 100.0 / len(fragms)), end='', flush=True)
			match = False
			if e in gs:
				sgrnas.append(e)
				match = True
				continue
			else:
				for r in gs:
					# Handle special cases (e.g. guides of length different from most other guides)
					if r == "AATGAACTACAATCCGGAATTGG":   # this guide is 3 nts longer
						if Levenshtein.distance(e + reads[i][43:46], r) <= 2 and Levenshtein.distance(reads[i][:23], up) < 5 and Levenshtein.distance(reads[i][46:66], dn) < 5:
							sgrnas.append(r)
							match = True
							break
					elif r == "GCTCAGCAGTGACCCTTATCTGG": # this guide is 3 nts longer
						if Levenshtein.distance(e + reads[i][43:46], r) <= 2 and Levenshtein.distance(reads[i][:23], up) < 5 and Levenshtein.distance(reads[i][46:66], dn) < 5:
							sgrnas.append(r)
							match = True
							break
					else:
						if Levenshtein.distance(e, r) <= 2 and Levenshtein.distance(reads[i][:23], up) < 5 and Levenshtein.distance(reads[i][43:66], dn) < 5:
							sgrnas.append(r)
							match = True
							break
			if not match:
				# Test if it's an empty vector
				if Levenshtein.distance(e, dn[:20]) < 5 and Levenshtein.distance(reads[i][:23], up) < 5:
					sgrnas.append('EMPTY')
				else:
					sgrnas.append('')
	else:
		print("WARNING: guide correction in the amplicon readout is turned off")
		sgrnas = fragms

	print("Assigned fragments for {} sequences".format(len(sgrnas)))
	print("Assigned guides for {} sequences".format(len([i for i in sgrnas if i != '' and i != 'EMPTY'])))


	return (sgrnas, fragms)
	

def parse_reads_with_barcodes(r1_file: str, guides_df: pd.DataFrame, barcodes_file: str, sgrnas: List[str], fragms: List[str], log: types.LambdaType = lambda x: None):
	lane_name = os.path.basename(r1_file).replace("_R1.fastq.gz", "")
	guides = guides_df.seq.values
	n_guides = len(guides)

	cb_titles = []
	cb_seqs = []

	print(f"Reading {r1_file}...")
	with gzip.open(r1_file, 'rt') as f:
		for title, seq, qual in FastqGeneralIterator(f):
			cb_titles.append(title)
			cb_seqs.append(seq[:16])
	log({"n_cellbarcodes": len(cb_seqs)})
	log({"n_cellbarcodes_unique": len(np.unique(cb_seqs))})
	print(f"Read {len(cb_seqs)} sequences with cell barcodes for {lane_name}")
	print(f"Found {len(np.unique(cb_seqs))} unique cell barcodes for {lane_name}")

	
	N = len(cb_seqs)

	assert N == len(fragms)
	assert N == len(sgrnas)

	# Read cell barcodes white list
	barcodes = pd.read_csv(barcodes_file, header=None).loc[:,0]
	barcodes = [e.split('-')[0] for e in barcodes]

	assert len(barcodes) == len(np.unique(barcodes))

	log({"n_barcodes_white_list": len(barcodes)})
	print(f"Loaded a white list with {len(barcodes)} cell barcodes")

	"""
	NOTE:
	- this step takes a lot of time
	- only work with barcodes that are present in the white list & detected, discard others
	"""

	cb_to_sgrna = {}
	cb_to_fragm = {}

	print("Counting guides per barcode...")
	for i, cb in enumerate(cb_seqs):
		if cb not in barcodes:
			continue
		if i % 100 == 0:
			print("Done {:.2f}%\r".format(i * 100.0 / N), end='', flush=True)
		cb_to_sgrna[cb] = cb_to_sgrna.get(cb, []) + [sgrnas[i]]
		cb_to_fragm[cb] = cb_to_fragm.get(cb, []) + [fragms[i]]

	print(f"Assigned R3 fragments to {len(cb_to_fragm.keys())} cell barcodes")
	log({"n_barcodes_with_r3_fragment_assigned": len(cb_to_fragm.keys())})


	tmp_mean = float(np.mean([len(v) for k, v in cb_to_fragm.items()]))
	print(f"Mean number of fragments assigned per barcode: {tmp_mean}")
	log({"mean_n_fragments_per_barcode": tmp_mean})

	cb_df = pd.DataFrame(columns=["barcode",
								  "sgrna", "sgrna_n_reads", "sgrna_frac", "sgrna_frac_se",
								  "fragment", "fragment_n_reads", "fragment_frac",
								  "n_reads", "n_sgrnas", "n_fragments",
								  "sgrna_snd", "sgrna_snd_n_reads", "sgrna_snd_frac"])

	for i, k in enumerate(cb_to_sgrna.keys()):
		if i % 100 == 0:
			print("Done {:.2f}%\r".format(i * 100.0 / len(cb_to_sgrna)), end='', flush=True)

		v_fr = cb_to_fragm[k]
		v    = cb_to_sgrna[k]

		sgrna_table = Counter(v).most_common()
		fragm_table = Counter(v_fr).most_common()

		sgrna_snd = None
		sgrna_snd_frac = None
		sgrna_snd_n_reads = None
		if len(sgrna_table) > 1:  # multiple guides are assigned to the barcode
			sgrna_snd         = sgrna_table[1][0]
			sgrna_snd_n_reads = sgrna_table[1][1]
			sgrna_snd_frac    = sgrna_table[1][1] * 1.0 / len(v)
		cb_df = cb_df.append({
			"barcode": k,

			"sgrna": sgrna_table[0][0],
			"sgrna_n_reads": sgrna_table[0][1],
			"sgrna_frac": sgrna_table[0][1] * 1.0 / len(v),
			"sgrna_frac_se": calculate_frac_se(sgrna_table[0][1], len(v)),

			"fragment": fragm_table[0][0],
			"fragment_n_reads": fragm_table[0][1],
			"fragment_frac": fragm_table[0][1] * 1.0 / len(v_fr),

			"n_reads": len(v),
			"n_sgrnas": len(sgrna_table),
			"n_fragments": len(fragm_table),

			"sgrna_snd": sgrna_snd,
			"sgrna_snd_n_reads": sgrna_snd_n_reads,
			"sgrna_snd_frac": sgrna_snd_frac
		}, ignore_index=True)



	cb_df = cb_df.sort_values("sgrna_n_reads", ascending=False)

	return cb_df



if __name__ == "__main__":

	# Parse input arguments

	list_of_guides_file = snakemake.input["list_of_guides"]
	r3_file = snakemake.input["r3_file"]
	r1_file = snakemake.input["r1_file"]
	barcodes_file = snakemake.input["barcodes_file"]
	output_file = snakemake.output[0]
	options, statlog = snakemake.config["preprocessing"]["amplicon"]["options"], snakemake.config["preprocessing"]["amplicon"]["statlog"]

	# Create a statlog file specifically for this lane
	lane_name = os.path.basename(r3_file).replace("_R3.fastq.gz", "")
	statlog = statlog.replace("LANE", lane_name)

	# Flush statlog file content if any
	open(statlog, 'w').close()

	log_func = lambda stats: write_log(statlog, stats) if statlog is not None else None

	# Read the list of guides
	guides_df = get_guides_df(list_of_guides_file, log=log_func)

	# Read the R3 part of reads (containing guide RNA sequences)
	sgrnas, fragms = parse_reads_with_guides(r3_file, guides_df, correct_fragments=True, log=log_func)

	# Read the R1 part of reads (containing cell barcodes)
	cb_with_guides = parse_reads_with_barcodes(r1_file, guides_df, barcodes_file, sgrnas, fragms, log=log_func)

	print(f"Saved output to {output_file}")
	cb_with_guides.to_csv(output_file, index=False)

