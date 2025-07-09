#!/usr/bin/env Rscript

#############################################################
# Performing Purity & Ploidy Estimation Using CLONETv2 & TPES
#############################################################

# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###################
# Loading Libraries
###################

library(data.table)
library(CLONETv2)
library(TPES)
library(ggplot2)

##############
# Loading Data
##############

control = fread("data/control.csv",data.table=F)
tumor = fread("data/tumor.csv",data.table=F)
seg.tb <- fread("data/SCNA.copynumber.called.seg",data.table=F)
snv.reads = fread("data/somatic.pm",data.table=F)

####################
# Data Preprocessing
####################

# Adding allele frequency (af) columns
control$af = control$altCount/control$totalCount
tumor$af = tumor$altCount/tumor$totalCount

# Preparing pileup data
pileup.control = control[,c(1,2,4,5,14,8)]
colnames(pileup.control) = c("chr","pos","ref","alt","af","cov")
pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")

# Preprocessing SNV reads
snv.reads = snv.reads[which(snv.reads$somatic_status=="Somatic"),]
snv.reads = snv.reads[,c("chrom","position","position","tumor_reads1","tumor_reads2")]
colnames(snv.reads) = c("chr","start","end","ref.count","alt.count")
snv.reads$sample = "Sample.1"

##################
# Computing Tables
##################

# Computing beta table
bt <- compute_beta_table(seg.tb, pileup.tumor, pileup.control)

# Computing ploidy, admixture and allele-specific CNA tables
pl.table <- compute_ploidy(bt)
adm.table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl.table)
allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
                                                                ploidy_table = pl.table, 
                                                                admixture_table = adm.table)

# Saving the beta, ploidy, admixture tables and allele-specific CNA table
write.table(bt, file="results/beta_table.txt", sep="\t", row.names=F, quote=F)
write.table(pl.table, file="results/ploidy_table.txt", sep="\t", row.names=F, quote=F)
write.table(adm.table, file="results/admixture_table.txt", sep="\t", row.names=F, quote=F)
write.table(allele_specific_cna_table, file="results/allele_specific_cna_table.txt", sep="\t", row.names=F, quote=F)

###################
# CLONETv2 Analysis
###################

# Generating the ploidy and admixture plot
check.plot <- check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl.table,
                                         admixture_table = adm.table)

# Plotting the results
print(check.plot)

# Saving the ploidy and admixture plot
ggsave("figures/ploidy_and_admixture_plot.png", plot = check.plot, width = 10, height = 6)

###############
# TPES Analysis
###############

# Performing purity estimation using TPES
purity_results <- TPES_purity(ID = "Sample.1", SEGfile = seg.tb,
                              SNVsReadCountsFile = snv.reads,
                              ploidy = pl.table,
                              RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

# Saving purity results
write.table(purity_results, file="results/TPES_purity_results.txt", sep="\t", row.names=F, quote=F)

# Generating TPES report and saving the purity plot
png("figures/TPES_purity_plot.png", width=2400, height=1800, res=300)
TPES_report(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

dev.off()
