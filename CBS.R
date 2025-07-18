#!/usr/bin/env Rscript

#####################################################################
# Circular Binary Segmentation (CBS) for Somatic Copy Number Analysis
#####################################################################

# Setting working directory
folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(folder)

# Loading libraries
library(DNAcopy)

# Loading SCNA copy number data
cn <- read.table(file.path(folder, "./data/SCNA.copynumber.called"), header = T)

# Creating a PDF to save plots
pdf(file.path(folder, "./figures/seg_plots.pdf"))

# Plotting raw log2 rations
plot(cn$raw_ratio, pch = ".", ylim = c(-2.5, 2.5))

# Plotting adjusted log2 ratios
plot(cn$adjusted_log_ratio, pch = ".", ylim = c(-2.5, 2.5))

# Constructing a CNA object using adjusted log2 ratios
CNA.object <- CNA(
    genomdat = cn$adjusted_log_ratio,
    chrom = cn$chrom,
    maploc = cn$chr_start, data.type = "logratio"
)

# Smoothing the CNA data to reduce noise
CNA.smoothed <- smooth.CNA(CNA.object)

# Performing segmentation using CBS algorithm
segs <- segment(
    CNA.smoothed,
    min.width = 2,
    undo.splits = "sdundo", # Undoes splits with SD difference < undo.SD
    undo.SD = 3,
    verbose = 1
)

# Plotting the segmentation results across all chromosomes
plot(segs, plot.type = "w")

# Closing the PDF device
dev.off()

# Extracting and saving the segmentation output
segs_output <- segs$output
write.table(segs_output, file = file.path(folder, "./data/SCNA.copynumber.called.seg"), row.names = F, col.names = T, quote = F, sep = "\t")
