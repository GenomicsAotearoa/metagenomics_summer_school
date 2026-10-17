#!/usr/bin/env python

'''tax_predict_vConTACT2.0.11.3.py

Generate viral taxonomy predictions for each contig based on guilt-by-association-clustering results from vConTACT2.

- Developed based on output from vConTACT2.0.11.3.


Optional parameters: 

--vcontact2_results (-i) genome_by_genome_overview.csv : output file from vcontact2. Default = './genome_by_genome_overview.csv'
--output_path (-o) output directory : Path for output directory. Default = './'
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np

parser = ArgumentParser()
parser.add_argument("-i", "--vcontact2_results", dest="infile", 
                    help="Output file from vcontact2: genome_by_genome_overview.csv",
                    metavar='vcontact2 genome_by_genome_overview.csv', required=True)
parser.add_argument("-o", "--output_path", dest="outpath", 
                    help="Path/to/output/directory. Default = current directory",
                    metavar='output_path', default='.')
args = parser.parse_args()

def main():
    print("\n--------------------\r\n")
    print("Running tax_predict.py\r\n")
    ## Read in vcontact2 results
    df = pd.read_csv(args.infile)
    df.columns = df.columns.str.replace(' ', '_')
    # Add subcluster column, and strip subcluster string ('_n') from VC column
    df['VC_Subcluster'] = df['VC']
    df['VC'] = df['VC'].str.replace(r'(VC_.*)_.*', r'\1')
    ## Generate taxonomy predictions for each viral cluster (one set for cluster, one set for subcluster)
    # For each taxonomic rank: Groupby viral cluster (VC), generate list of unique taxonomies associated with that cluster (excluding 'Unassigned'), add <rank>_VC_predicted to each contig in that cluster.
    # Replace instances of no predicted taxonomy for viral cluster (empty lists in predicted taxonomy column) with 'Unassigned'
    for taxa in ['Order', 'Family', 'Genus']:
        for VC_group in ['VC', 'VC_Subcluster']:
            df = df.join(df.groupby(by=VC_group)[taxa].apply(lambda s: list({x for x in s if x != "Unassigned"})), on=VC_group, rsuffix='_'+VC_group+'_predicted')
            df[taxa+'_'+VC_group+'_predicted'] = df[taxa+'_'+VC_group+'_predicted'].where(df[taxa+'_'+VC_group+'_predicted'].str.len() > 0, 'Unassigned')
    ## Add taking into account p-values or confidence scores?
    #
    ## Write out tables
    # genome_by_genome
    df.to_csv(args.outpath+'/genome_by_genome_overview_tax_predictions.tsv', sep='\t', index=False)
    # tax predictions
    subset_columns = ['Genome', 'Order_VC_predicted', 'Family_VC_predicted', 'Genus_VC_predicted', 'Order_VC_Subcluster_predicted', 'Family_VC_Subcluster_predicted', 'Genus_VC_Subcluster_predicted', 'VC', 'VC_Subcluster', 'VC_Status']
    df[subset_columns].to_csv(args.outpath+'/tax_predict_table.tsv', sep='\t', index=False)
    print("Output:\r\n")
    print(args.outpath+"/tax_predict_table.tsv:\nTable of predicted taxonomy for all contigs\r\n")
    print(args.outpath+"/genome_by_genome_overview_tax_predictions.tsv:\nReproduced genome_by_genome_overview table with the addition of predicted taxonomy for all contigs\r\n")
    print("Completed tax_predict.py\r\n")
    print("\n--------------------\r\n")

    
if __name__ == '__main__':
    main()

