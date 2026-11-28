#!/usr/bin/env python

# vcontact2_tax_predict.py
## Generate viral taxonomy predictions for each contig based on guilt-by-association-clustering results from vConTACT2.

# Optional parameters: 
## --vcontact2_results (-i) genome_by_genome_overview.csv : output file from vcontact2. Default = './genome_by_genome_overview.csv'
## --output_path (-o) output directory : Path for output directory. Default = './'


####################
# Importing libraries and setting input arguments
####################

#!/usr/bin/env python

# import packages
from argparse import ArgumentParser
import pandas as pd
import numpy as np
#import re

# Set arguments
parser = ArgumentParser()
parser.add_argument("-i", "--vcontact2_results", dest="infile", 
                    help="Output file from vcontact2: genome_by_genome_overview.csv. Default = ./genome_by_genome_overview.csv",
                    metavar='vcontact2_results', default='./genome_by_genome_overview.csv')
parser.add_argument("-o", "--output_path", dest="outpath", 
                    help="Path/to/output/directory/. Default = current directory",
                    metavar='output_path/', default='./')
args = parser.parse_args()


####################
# Main script
####################

print("\n####################\r\n")
print("Running vcontact2_tax_predict.py\r\n")

# Read in vcontact2 results
vcontact2_results = pd.read_csv(args.infile, index_col=0, dtype=object)

## Generate taxonomy predictions for each viral cluster
# For each taxonomic rank: Groupby by viral cluster (VC), generate list of unique taxonomies associated with that cluster (excluding 'Unassigned'), add <rank>_VC_predicted to each contig in that cluster.
vcontact2_results = vcontact2_results.join(vcontact2_results.groupby(by='VC')['Order'].apply(lambda s: list({x for x in s if x != "Unassigned"})), on='VC', rsuffix='_VC_predicted')
vcontact2_results = vcontact2_results.join(vcontact2_results.groupby(by='VC')['Family'].apply(lambda s: list({x for x in s if x != "Unassigned"})), on='VC', rsuffix='_VC_predicted')
vcontact2_results = vcontact2_results.join(vcontact2_results.groupby(by='VC')['Genus'].apply(lambda s: list({x for x in s if x != "Unassigned"})), on='VC', rsuffix='_VC_predicted')

# Replace instances of no predicted taxonomy for viral cluster (empty lists in predicted taxonomy column) with 'Unassigned'
vcontact2_results['Order_VC_predicted'] = vcontact2_results['Order_VC_predicted'].where(vcontact2_results['Order_VC_predicted'].str.len() > 0, 'Unassigned')
vcontact2_results['Family_VC_predicted'] = vcontact2_results['Family_VC_predicted'].where(vcontact2_results['Family_VC_predicted'].str.len() > 0, 'Unassigned')
vcontact2_results['Genus_VC_predicted'] = vcontact2_results['Genus_VC_predicted'].where(vcontact2_results['Genus_VC_predicted'].str.len() > 0, 'Unassigned')

# Subset only columns of interest.
vcontact2_results_sub  = vcontact2_results[['Genome', 'Order_VC_predicted', 'Family_VC_predicted', 'Genus_VC_predicted', 'VC', 'VC Subcluster', 'VC Status']]

# Write full table to file
vcontact2_results.to_csv(args.outpath+'genome_by_genome_overview_tax_predictions.txt', sep='\t', index=False)

# Write subset table of taxonomy predictions to file
vcontact2_results_sub.to_csv(args.outpath+'tax_predict_table.txt', sep='\t', index=False)


print("Output:\r\n")
print(args.outpath+"tax_predict_table.txt:\nTable of predicted taxonomy for all contigs\r\n")
print(args.outpath+"genome_by_genome_overview_tax_predictions.txt:\nReproduced genome_by_genome_overview table with the addition of predicted taxonomy for all contigs\r\n")
print("Completed vcontact2_tax_predict.py\r\n")
print("####################\r\n")
