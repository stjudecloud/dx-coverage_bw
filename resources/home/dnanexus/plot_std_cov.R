#!/usr/bin/env Rscript
#
# Plots standard coverage
# 
# arg1: File containing coverage values. This file should contain the following 
## columns: 
## col1: sample name
## col2: 10x coverage
## col3: 20x coverage
## col4: 30x coverage
# arg2: Sequence type (Ex: EXCAP, WGS, RNASEQ, etc.)
# arg3: Tumor cutoff (Ex: 30 for RNASEQ and 80 for others)
# arg4: Normal cutoff (Ex: 30 for RNASEQ and 70 for others)
# arg5: Output file

# Load required libraries
require(ggplot2)
require(reshape)

# Parse the arguments
args <- (commandArgs(TRUE))
data_file <- args[1]
type <- args[2]
tumor_cutoff <- as.numeric(args[3])
normal_cutoff <- as.numeric(args[4])
output_file <- args[5]

# Create data frame
cov <- read.table(data_file, sep="\t", header=T)
m.cov <- melt(cov)

# Plot the coverage and save it to a PNG file
if (type == "EXCAP" | type == "FREQEXCAP") {
  y_axis_label <- "% coding exon bases covered"
} else {
  y_axis_label <- "% exon bases covered"
}
tumor_cutoff_label <- paste("Tumor cutoff = ", tumor_cutoff, "%", sep="")
normal_cutoff_label <- paste("Normal cutoff = ", normal_cutoff, "%", sep="")
png(file=output_file, width=1000, height=500, units="px")
ggplot(m.cov, aes(Sample, value, fill=variable)) + 
  geom_bar(stat="identity", position="identity") +
  opts(axis.text.x=theme_text(angle=90)) +
  scale_x_discrete(expand=c(0.01,0.01)) +
  scale_y_continuous(y_axis_label, limits=c(0,100), expand=c(0,0)) +
  geom_hline(y=tumor_cutoff, color="firebrick2") +
  geom_hline(y=normal_cutoff, color="goldenrod") +
  annotate("text", x=1.3, y=tumor_cutoff+2.5, label=tumor_cutoff_label, color="gray40", fontface="italic", size=3) +
  annotate("text", x=1.3, y=normal_cutoff-2.5, label=normal_cutoff_label, color="gray40", fontface="italic", size=3) +
  scale_fill_discrete(name="Coverage depth", breaks=c("X10x", "X20x", "X30x"), labels=c(">=10x", ">=20x", ">=30x"))
dev.off()
