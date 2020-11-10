#!/usr/bin/env python

# norm_jgi_cov_table.py
## Calculate normalised average fold coverage for contigs based on output from jgi_summarize_bam_contig_depths
## Per-sample read depths are retrieved from the standard error output of the prior bowtie2 read mapping

# Required parameters: 
## --cov_table (-c) cov_table.txt : coverage table output from jgi_summarize_bam_contig_depths
## --bowtie2_err (-e) mapping.err : output file of captured standard error output from bowtie2 read mapping

# Optional parameters:
## --output_path (-o) output_directory : Path for output directory. Default = './'


####################
# Importing libraries and setting input arguments
####################

# import packages
from argparse import ArgumentParser
import os
import pandas as pd
import numpy as np
import re

# Set arguments
parser = ArgumentParser()
parser.add_argument("-c", "--cov_table", dest="cov_table_in", 
                    help="Coverage table output from jgi_summarize_bam_contig_depths",
                    metavar='cov_table', required=True)
parser.add_argument("-e", "--bowtie2_err", dest="bowtie2_err_in", 
                    help="Output file of captured standard error output from bowtie2 read mapping",
                    metavar='bowtie2_err', required=True)
parser.add_argument("-o", "--output_path", dest="out_path", 
                    help="Path/to/output/directory/. Default = current directory",
                    metavar='output_path', default='./')
args = parser.parse_args()


####################
# Main script
####################

print("\nRunning norm_jgi_cov_table.py\r\n")

## Calculate average reads/sample (library size)

# Generate empty list of read counts and add results for each sample via looping over the relevant lines in bowtie2_err
read_counts = []
with open(args.bowtie2_err_in, 'r') as readfile:
    for line in readfile:
        if "reads; of these:" in line:
            sample_read_count = int(re.sub(r'(\d+).*', r'\1', line))
            read_counts.append(sample_read_count)

# Calculate average library size
avg_read_counts = round(sum(read_counts) / len(read_counts))

## Calculate normalised coverage per sample and generate normalised count table

# covstats.txt files
# files = [entry.path for entry in filepath if '.covstats.txt' in entry.name]

# Read in coverage table
cov_table = pd.read_csv(args.cov_table_in, sep='\t')

# For each column that ends in .bam (the coverage values for each sample), normalise values by: coverage value / sample read count (from read_counts list) * avg_read_counts
i=0
for column in cov_table.columns[cov_table.columns.str.endswith('bam')]:
    cov_table.loc[:, column] = round((cov_table.loc[:, column] / read_counts[i]) * avg_read_counts, 4)
    i += 1

# Write output table to file
cov_table.to_csv(args.out_path+'normalised_'+args.cov_table_in, sep='\t', index=False)

print("Output:\r\n")
print(args.out_path+'normalised_'+args.cov_table_in+' : coverage table output from jgi_summarize_bam_contig_depths, normalised by average sample read depth.\r\n')
print("Completed norm_jgi_cov_table.py\r\n")


