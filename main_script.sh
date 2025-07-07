#!/bin/bash

# Script to execute the commands for the computational human genomics project.
# Ensure that the necessary tools and files are available in the specified paths before running this script.

##################
# 1. Initial Setup
##################

# Moving to the data directory
cd ./data

# Sorting the bam files
samtools sort control.bam > control.sorted.bam; samtools sort tumor.bam > tumor.sorted.bam

# Creating index files for the sorted bam files
samtools index control.sorted.bam; samtools index tumor.sorted.bam

##############################
# 2. Realignment Around Indels
##############################

# Control
# Creating a bed file with the regions of interest
java -jar ../tools/genome_analysis_TK.jar -T RealignerTargetCreator -R ../annotations/human_g1k_v37.fasta -I ./control/control.sorted.bam -o ./control/control.realigner.intervals -L ../annotations/captured_regions.bed
# Realigning the bam file around indels
java -jar ../tools/genome_analysis_TK.jar -T IndelRealigner -R ../annotations/human_g1k_v37.fasta -I control.sorted.bam -targetIntervals control.realigner.intervals -o control.sorted.realigned.bam -L ../annotations/captured_regions.bed

# Tumor
java -jar ../tools/genome_analysis_TK.jar -T RealignerTargetCreator -R ../annotations/human_g1k_v37.fasta -I tumor.sorted.bam -o tumor.realigner.intervals -L ../annotations/captured_regions.bed
java -jar ../tools/genome_analysis_TK.jar -T IndelRealigner -R ../annotations/human_g1k_v37.fasta -I tumor.sorted.bam -targetIntervals tumor.realigner.intervals -o tumor.sorted.realigned.bam -L ../annotations/captured_regions.bed

############################################
# 3. Base Quality Score Recalibration (BQSR)
############################################

# Control
# Building recalibration model
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./control.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -o control.recal.table -L ../annotations/captured_regions.bed
# Creating a new bam file using the input table generated from the previous command which has accurate base substitution, insertion and deletion quality scores
java -jar ../tools/genome_analysis_TK.jar -T PrintReads -R ../annotations/human_g1k_v37.fasta -I ./control.sorted.realigned.bam -BQSR control.recal.table -o control.sorted.realigned.recalibrated.bam -L ../annotations/captured_regions.bed --emit_original_quals
# Second pass to evaluate what the data looks like after recalibration
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./control.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -BQSR control.recal.table -o control.after_recal.table -L ../annotations/captured_regions.bed
# Generating csv file based on before and after recalibration tables
java -jar ../tools/genome_analysis_TK.jar -T AnalyzeCovariates -R ../annotations/human_g1k_v37.fasta -before control.recal.table -after control.after_recal.table -csv control.recal.csv

# Tumor
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./tumor.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -o tumor.recal.table -L ../annotations/captured_regions.bed
java -jar ../tools/genome_analysis_TK.jar -T PrintReads -R ../annotations/human_g1k_v37.fasta -I ./tumor.sorted.realigned.bam -BQSR tumor.recal.table -o tumor.sorted.realigned.recalibrated.bam -L ../annotations/captured_regions.bed --emit_original_quals
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./tumor.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -BQSR tumor.recal.table -o tumor.after_recal.table -L ../annotations/captured_regions.bed
java -jar ../tools/genome_analysis_TK.jar -T AnalyzeCovariates -R ../annotations/human_g1k_v37.fasta -before tumor.recal.table -after tumor.after_recal.table -csv tumor.recal.csv

# Executing R script to visualize the BQSR results
Rscript ../bqsr_plot_generator.R

######################################
# 4. Identifying & Removing Duplicates
######################################

# Control
# Removing duplicates
java -jar ../tools/picard.jar MarkDuplicates I=control.sorted.realigned.recalibrated.bam O=control.sorted.realigned.recalibrated.dedup.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=control.sorted.realigned.recalibrated.picard.log ASSUME_SORTED=true
# Reindexing the new deduplicated bam file
samtools index control.sorted.realigned.recalibrated.dedup.bam
# Generating flagstat report for the deduplicated bam file
samtools flagstat control.sorted.realigned.recalibrated.dedup.bam > control.sorted.realigned.recalibrated.dedup.flagstat.txt

# Tumor
java -jar ../tools/picard.jar MarkDuplicates I=tumor.sorted.realigned.recalibrated.bam O=tumor.sorted.realigned.recalibrated.dedup.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=tumor.sorted.realigned.recalibrated.picard.log ASSUME_SORTED=true
samtools index tumor.sorted.realigned.recalibrated.dedup.bam
samtools flagstat tumor.sorted.realigned.recalibrated.dedup.bam > tumor.sorted.realigned.recalibrated.dedup.flagstat.txt
