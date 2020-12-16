# CRISPRa single-cell screen

This is the code that has been used to analyze a CRISPR-activation screen with single-cell RNA-seq readout in order to find regulators of zygotic genome activation. Peer-reviewed publication that describes this dataset [can be found here](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30201-5) (Open Access).

The main analysis involves assigning sgRNAs to cells as well as pre-processing scRNA-seq data.

## Data overview

Raw data is of two main types: global transcriptome read-out (scRNA-seq, 10X Genomics) and amplicon sequencing of the same libraries to get guide &rarr; cell assignment. They are to be located in the directory `data/raw/main/scrnaseq/`, with Cell Ranger output (a folder per sample) in the `transcriptome` subfolder and amplicon sequencing FASTQ files in the `amplicon` subfolder. This repository does not contain raw data due to its large size. Original data has been deposited on GEO under [GSE135621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135621). 

Filtered cell barcodes as defined by Cell Ranger are made available [here on figshare](https://figshare.com/s/50500c2216d7ba617876):

```
curl 'https://ndownloader.figshare.com/articles/13393214/versions/1' -o GSE135621_filtered_barcodes.zip
```

Original tables such as a list of sgRNAs can be found in [`data/raw/main/tables`](data/raw/main/tables) folder.

## Processing

Main processing steps are described in the `workflow/Snakefile`, which is also symlinked in the root folder. Processing scripts are to be found in the `src` folder and its subfolders.

## Repeat elements quantification

FASTA files with repeat elements sequences per family are located in `data/external/repeats`. Due to their filesize, the files are not present in this repository and can be obtained [from figshare](https://figshare.com/s/24f9a37e19fed9258338):

```
curl 'https://ndownloader.figshare.com/files/23166317' -o repeats_fa.tar.gz
```

FASTQ files with cell barcodes are created from the unmapped reads of the BAM files in the Cell Ranger output and are then mapped to those families. Only reads mapped to a respective family are kept for downstream quantification.

Scripts to map originally unmapped reads to repeat elements and to aggregate results are located in the [`src/data/repeats`](src/data/repeats) directory.
