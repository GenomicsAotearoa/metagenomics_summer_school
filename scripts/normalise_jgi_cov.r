#!/usr/bin/env Rscript

# Description:
# This script normalises results from jgi_summarize_bam_contig_depths by read
# depth and scales values by average library size.
#
# Usage:
# normalise_jgi_cov.r <coverage_table> <library_size>
# Coverage table must be generated from jgi_summarize_bam_contig_depth
# Library size must have BAM file name as the first column and library size as
# second column.


# Read command-line arguments
args = commandArgs(trailingOnly=TRUE)

# Parse arguments
covFile <- args[1]
libFile <- args[2]

# Check arguments
if (length(args) != 2) {
  stop(
    "Requires coverage table and library size files.\nUsage: normalise_jgi_cov.r <coverage_table> <library_size>"
  )
}

# Print status
cat("Running normalise_jgi_cov.r with the following files:\n")
cat(paste("\tCoverage table:", args[1], "\n"))
cat(paste("\tLibrary sizes:", args[2], "\n"))

# Read data
cov_table <- read.delim(args[1])
libsize_table <- read.delim(args[2], header = F)

# Check library size file
libsize_class <- unlist(
  lapply(libsize_table, class)
)

if (!all(libsize_class == c("character", "integer"))) {
  stop(
    "Library size file in wrong format. First column should be bam filename, second column should be library size (integer)"
  )
}

# Check samples
samples_libsize <- sort(libsize_table[, 1])
samples_cov <- sort(grep(".*.bam$", names(cov_table), value = T))

diff_samples <- setdiff(samples_libsize, samples_cov)

if (length(diff_samples) > 0) {
  stop(
    paste("These sample names do not match, please input files:", diff_samples)
  )
}

# Wrangle data
cov_matrix <- as.matrix(cov_table[, samples_libsize])
rownames(cov_matrix) <- cov_table[, "contigName"]

libsize_vector <- libsize_table[order(libsize_table[, 1]), 2]
names(libsize_vector) <- sort(libsize_table[, 1])

# Normalize data
cat(paste("Average library size:", mean(libsize_vector), "\n"))

norm_cov <- (cov_matrix / libsize_vector) * mean(libsize_vector)

# Output data
out_table <- cbind(
  cov_table[, c("contigName", "contigLen")],
  round(norm_cov, 4)
)

outName <- paste0("normalised_", basename(args[1]))

cat(paste0("Writing output to ./", outName, "\n"))

write.table(out_table, file = outName, sep = "\t", quote = FALSE,
            row.names = FALSE)