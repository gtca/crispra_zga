# CRISPRa single-cell screen

This is the code that has been used to analyze a CRISPR-activation screen with single-cell RNA-seq readout in order to find regulators of zygotic genome activation.

The main analysis involves assigning guide to cells as well as pre-processing scRNA-seq data.

## Data overview

Raw data is of two main types: global transcriptome read-out (scRNA-seq, 10X Genomics) and amplicon sequencing of the same libraries to get guide &rarr; cell assignment. They are to be located in the directory `data/raw/main/scrnaseq/`, with Cell Ranger output (a folder per sample) in the `transcriptome` subfolder and amplicon sequencing FASTQ files in the `amplicon` subfolder. This repository does not contain raw data due to its large size.
