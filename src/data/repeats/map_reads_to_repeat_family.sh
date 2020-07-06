#!/bin/bash

# 
# Launch in a loop for samples:
# for sample in $(ls -d data/raw/main/scrnaseq/transcriptome/*/outs); do sh src/data/repeats/map_reads_to_repeat_family.sh $sample; done
#

#
# Available families:
# LINE_L1 LINE_L2 LTR_ERV1 LTR_ERVK TRVLLTR_Emajor_satllite minor_satellite rRNA SINGE_Alu SINE_B2 SINE_B4 telomeric_repeats MERVL
#

families=$(ls data/external/repeats/*fa | xargs -n1 basename | sed 's/\..*$//')

sample=$1  # CellRanger SIGAXXX/outs folder path
unaligned_bam=$sample/unaligned.bam
unaligned_fq=$sample/unaligned.fq

# Extracted unaligned reads
echo "Selecting unmapped reads..."
samtools view -f4 -b $sample/possorted_genome_bam.bam > $unaligned_bam

# Convert to fastq
# NOTE: bam_to_fq_with_cb has to be compiled
echo "Getting fastq files..."
src/data/bam_to_fq_with_cb $unaligned_bam > $unaligned_fq

# Prepare the index for the genome
echo "Preparing the index for repeat family..."
for family in $(echo $families | xargs -n1); do
	echo Preparing the index for the repeat family "$family"
	bwa index -p data/external/repeats/"$family" -a is $(ls data/external/repeats/"$family"*fa)
	echo Processing data/external/repeats/"$family"...
done

# Align reads
echo "Aligning reads..."
for family in $(echo $families | xargs -n1); do
	echo Processing data/external/repeats/"$family"...
	bwa aln -t 4 data/external/repeats/"$family" $unaligned_fq > $sample/aligned_"$family".sai
	bwa samse data/external/repeats/"$family" $sample/aligned_"$family".sai $unaligned_fq > $sample/aligned_"$family".sam
	echo "Done. Check $sample/aligned_"$family".sam file."
done

# Take only mapped reads
echo "Filtering out mapped reads"
for family in $(echo $families | xargs -n1); do
	samtools view -F4 $sample/aligned_"$family".sam > $sample/aligned_"$family"_mapped.sam 
	echo "Done. Check $sample/aligned_"$family".sam file."
done

