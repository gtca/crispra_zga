# CRISPRa single-cell screen

This is the code that has been used to analyze a CRISPR-activation screen with single-cell RNA-seq readout in order to find regulators of zygotic genome activation.

The main analysis involves assigning sgRNAs to cells as well as pre-processing scRNA-seq data.

## Data overview

Raw data is of two main types: global transcriptome read-out (scRNA-seq, 10X Genomics) and amplicon sequencing of the same libraries to get guide &rarr; cell assignment. They are to be located in the directory `data/raw/main/scrnaseq/`, with Cell Ranger output (a folder per sample) in the `transcriptome` subfolder and amplicon sequencing FASTQ files in the `amplicon` subfolder. This repository does not contain raw data due to its large size.

Original tables such as a list of sgRNAs can be found in [`data/raw/main/tables`](data/raw/main/tables) folder.

## Processing

Main processing steps are described in the `workflow/Snakefile`, which is also symlinked in the root folder. Processing scripts are to be found in the `src` folder and its subfolders.

## Repeat elements quantification

FASTA files with repeat elements sequences per family are located in `data/external/repeats`. Due to their filesize, the files are not present in this repository and can be obtained [from figshare](https://figshare.com/s/24f9a37e19fed9258338).

FASTQ files with cell barcodes are created from the unmapped reads of the BAM files in the Cell Ranger output and are then mapped to those families. Only reads mapped to a respective family are kept for downstream quantification.
