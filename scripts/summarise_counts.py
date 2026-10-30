#!/usr/bin/env python

'''summarise_counts.py

Generate summary table of normalised coverage based on read mapping data, with additional functionality to run edgeR differential expression analyses (via summarise_counts.R).

Method:
- Extract input raw counts:
  - Extract from pileup.sh \*_fpkm.txt output: ContigID, length, read_count, fragments
  - Or extract read columns from feature counts output.
- Zero out read counts < set threshold (default = 10 reads)
- Re-calculate RPKM and FPKM based on the lib_norm option (normalised total library or mapped reads). 
  - Note: for consistency with pileup.sh output and comparability with RPKM values, FPKM is multiplied by 2.
  - Note: FPKM not calculated if --format 'featurecounts' (fragment counts not provided)
- TPM calculated
  - Take library size from \*_fpkm.txt header line: #READS (--format pileup) or from provided mapping file (--format featurecounts)
  - Note: Due to how TPM is calculated, this is normalised based on *mapped* reads, rather than incorporating total library size, regardless of the setting of --lib_norm
- TMM calculated (via R: edgeR) if sample mapping file provided
- edgeR pairwise differential expression analyses conducted if sample mapping file provided and --edgeR_out set
  - edgeR testing notes: 
    - differential expression tested via edgeR glmQLFit() and glmQLFTest()
    - design = model.matrix(~0+Groups)
    - group pairwise combinations generated via makeContrasts()
    - Each pairwise combination tested separately and appended to master summary table
    - FDR calculated independently for each tested pair
    - Outputs summary table of concatenated glmQLFTest results
- Genome-level count summaries calculated if --genome_mapping_file provided
- Output summary table of counts per sample for each metric

Required parameters:
--input (-i) : rpkm output file(s) from pileup (Takes wildcard to generate a list of input files for all samples, e.g. <-i '*rpkm.txt'>), or output file from featureCounts
--format (-f) : format of input counts file (options: 'pileup', 'featurecounts')

Optional parameters:
--lib_norm (-n) : Set whether to (re)calculate RPKM and FPKM based on total library size or mapped reads per sample. Options = ['total', 'mapped']. (Default = 'mapped')", default='mapped')
--count_threshold (-t) : Set threshold to zero out low count values < threshold (Default = 1)"
--read_counts (-r) : Output file name for table of read counts (per-sample total library, mapped read, and filtered mapped read (based on count_threshold) counts)
--sample_mapping_file (-s) : Sample mapping file (cols: sampleID = unique substrings; group = sample group). TMM calculated if sample mapping file provided. edgeR differential expression analysis conducted if sample mapping file provided and --edger_out is set.
--genome_mapping_file (-g) : genomeID (or magID or binID) to contigID mapping file (two columns: (genomeID or magID or binID); contigID). If provided, genome-level summaries are calculated and output as a separate file (<output_filename>_genomeSummary.tsv).
--genome_libSize_norm (-m) : Set whether to calculate genome-level summary read depth normalisation based on min or mean sample read depth. Options = ['min', 'mean']. (Default = 'min').
--edger_out (-e) : Output file name for edgeR pairwise differential expression analyses summary table. Requires both --mapping_file and --edger_out parameters to be provided. 
--output (-o) : Output file name for summary count table (default = 'normalised_summary_count_table.tsv')
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import re
from pathlib import Path
import subprocess
import os

parser = ArgumentParser()
parser.add_argument('-i', '--input', dest="input",
                    help="rpkm output file(s) from pileup (Takes wildcard to generate a list of input files for all samples, e.g. <-i *rpkm.txt>), or output file from featureCounts", required=True)
parser.add_argument('-f', '--format', dest="format",
                    help="Format of input counts file (options = ['pileup', 'featurecounts'])", required=True)
parser.add_argument('-n', '--lib_norm', dest="lib_norm",
                    help="Set whether to calculate normalisations based on total library size or mapped reads per sample. Options = ['total', 'mapped']. (Default = 'mapped')", default='mapped')
parser.add_argument('-t', '--count_threshold', dest="count_threshold",
                    help="Set threshold to zero out low count values < threshold (Default = 10)", default=1)
parser.add_argument('-r', '--read_counts', dest="read_counts",
                    help="Output file name for table of read counts (per-sample total library and mapped read counts)  (default = 'summary_read_counts.tsv')", default='summary_read_counts.tsv')
parser.add_argument('-s', '--sample_mapping_file', dest="sample_map",
                    help="Sample mapping file (cols: sampleID = unique substrings; group = sample group). TMM calculated if sample mapping file provided.", default=None)
parser.add_argument('-g', '--genome_mapping_file', dest="genome_map",
                    help="genomeID (or magID or binID) to contigID mapping file (two columns: (genomeID or magID or binID); contigID). If provided, genome-level summaries are calculated and output as a separate file (<output_filename>_genomeSummary.tsv).", default=None)
parser.add_argument('-m', '--genome_libSize_norm', dest="genome_libSize_norm",
                    help="Set whether to calculate genome-level summary read depth normalisation based on min or mean sample read depth. Options = ['min', 'mean']. (Default = 'min')", default='min')
parser.add_argument('-e', '--edger_out', dest="edger_out",
                    help="Output file name for edgeR pairwise differential expression analyses summary table. Requires both --sample_mapping_file and --edger_out parameters to be provided.", default=None)
parser.add_argument('-o', '--output', dest="output",
                    help="Output file name for summary count table (default = 'normalised_summary_count_table.tsv')", default='normalised_summary_count_table.tsv')
args = parser.parse_args()


def norm_from_featurecounts():
    # Read in input file
    input_df = pd.read_csv(args.input, sep='\t', skiprows=1, low_memory=False)
    # Extract counts columns names
    counts_cols = [names for names in input_df.columns if names not in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']]
    # Establish empty summary_count_table and read_counts DataFrames
    summary_count_table = pd.DataFrame()
    read_counts = pd.DataFrame()
    # loop through each counts column in the featurecounts input file
    # For each iteration:
    ## zero out counts below the count_treshold
    ## extract sample full library size from mapping file and calculate mapped reads counts (select for lib_norm based on args.lib_norm)
    ## calculate RPKM, TPM
    ## merge with summary_count_table DataFrame
    ## concatenate read counts onto read_counts
    print('Calculating normalisation: RPKM')
    print('Calculating normalisation: TPM')
    for counts_col in counts_cols:
        # Extract sampleID from counts_col
        sampleID = re.sub(r".*/(.*).bam", r'\1', counts_col)
        # Subset featurecounts input file for sample
        sample_counts_df = input_df[['Geneid','Chr', 'Length', counts_col]].copy()
        # update col names in sample_counts_df
        sample_counts_df.rename(columns = {'Geneid': 'geneID', 'Chr': 'contigID', counts_col: 'counts'}, inplace = True)
        # Save library size (extracted from sample mapping file, if provided)
        map_df = pd.read_csv(args.sample_map, sep='\t')
        if 'lib.size' not in map_df.columns:
            if args.lib_norm == 'total':
                print('Error: lib_norm = "total", but "lib.size" not included in sample mapping file. Add lib.size column to sample mapping file and re-run, or select lib_norm = "mapped". Exiting.')
                exit()
            else:
                libSize = 'Not provided'
        else:
            libSize = map_df[map_df['sampleID'].str.contains(sampleID, na=False)]['lib.size'].values[0]
        # Calculate mapped read count
        mapped_libSize = sample_counts_df['counts'].sum()
        # Zero out counts below the count_treshold
        sample_counts_df['counts'] = sample_counts_df['counts'].where(sample_counts_df['counts']>=int(args.count_threshold), other=0)
        # Calculate filtered mapped read count
        mapped_libSize_filt = sample_counts_df['counts'].sum()
        # Set lib_norm value based on args.lib_norm setting
        if args.lib_norm == 'mapped':
            lib_norm = mapped_libSize_filt
        elif args.lib_norm == 'total':
            lib_norm = libSize
        else:
            print('Error: --lib_norm setting invalid (must be one of: ["mapped", "total"]). Exiting.')
            exit()
        # Re-calculate RPKM (based on args.lib_norm setting and mapped read threshold filtereing)
        sample_counts_df['RPKM'] = round((sample_counts_df['counts']/(int(lib_norm)/1000000))/(sample_counts_df['Length']/1000), 4)
        # n.b. FPKM calculation omitted from featurecounts data (fragments counts not provided)
        # Calculate TPM (Note: Due to how TPM is calculated, this is normalised based on *mapped* reads, rather than incorporating total library size, regardless of the setting of args.lib_norm)
        rpk_tmp = sample_counts_df['counts']/(sample_counts_df['Length']/1000)
        sample_counts_df['TPM'] = rpk_tmp/(rpk_tmp.sum()/1000000)
        # Rename 'Reads' column as 'featurecounts_counts'
        sample_counts_df.rename(columns = {'counts': 'featurecounts_counts'}, inplace = True)
        # Add sampleID prefix to column headers (except 'contigID' and 'Length')
        sample_counts_df.rename(columns={col: str(sampleID)+'_'+col for col in sample_counts_df.columns if col not in ['geneID', 'contigID', 'Length']}, inplace = True)
        # Merge into summary_count_table
        if summary_count_table.empty:
            summary_count_table = sample_counts_df.copy()
        else:
            sample_counts_df.drop(columns='Length', inplace = True)
            summary_count_table = summary_count_table.merge(sample_counts_df, left_on=['geneID', 'contigID'], right_on=['geneID', 'contigID'], how = 'outer')
        # Add read counts to read_counts DataFrame
        if read_counts.empty:
            read_counts = pd.DataFrame([[sampleID, libSize, mapped_libSize, mapped_libSize_filt]],
                                       columns=['SampleID', 'Total_library_size', 'Mapped_reads', 'Filtered_mapped_reads'])
        else:
            read_counts = read_counts.append(pd.DataFrame([[sampleID, libSize, mapped_libSize, mapped_libSize_filt]], columns=['SampleID', 'Total_library_size', 'Mapped_reads', 'Filtered_mapped_reads']), ignore_index=True, sort=False)
    # Rearrange summary_count_table column order
    summary_count_table = summary_count_table[['geneID', 'contigID', 'Length'] + \
                                              [col for col in summary_count_table.columns if 'counts' in col] + \
                                              [col for col in summary_count_table.columns if 'RPKM' in col] + \
                                              [col for col in summary_count_table.columns if 'TPM' in col]
                                             ]
    # Write out summary_count_table to csv
    summary_count_table.to_csv(args.output, sep='\t', index=False)
    # If args.read_counts is set, write out read_counts to csv
    read_counts.to_csv(args.read_counts, sep='\t', index=False)


def norm_from_pileup():
    # Generate list of input files and save path to input files (if input files located in a subdirectory)
    in_files = []
    if '/' in args.input:
        in_path = args.input.rsplit('/', 1)[0] + '/'
    else:
        in_path = ''
    for path in Path(in_path).rglob(args.input.rsplit('/', 1)[-1]):
        in_files.append(path.name)
    in_files.sort()
    # Establish empty summary_count_table and read_counts DataFrames
    summary_count_table = pd.DataFrame()
    read_counts = pd.DataFrame()
    # loop through each input file
    # For each iteration:
    ## zero out counts below the count_treshold
    ## extract sample full library size and calculate mapped reads counts (select for lib_norm based on args.lib_norm)
    ## (re)calculate RPKM, FPKM, TPM
    ## rename headers based on input file name (proxy for sample name)
    ## merge with summary_count_table DataFrame
    ## concatenate read counts onto read_counts
    print('Calculating normalisation: RPKM')
    print('Calculating normalisation: FPKM')
    print('Calculating normalisation: TPM')
    for in_file in in_files:
        # Save sampleID from filename 
        # if wildcard supplied, strip off the wildcard-matched suffix segment)
        if '*' in args.input:
            sampleID = re.sub(rf"{args.input.rsplit('/', 1)[-1].strip('*')}", '', in_file)
        else:
            sampleID = in_file
        # Input sample file (skip 4 x extra header rows)
        input_df = pd.read_csv(in_path+in_file, sep='\t', skiprows=4, low_memory=False)
        # Rename '#Name' to 'contigID', and 'Coverage' to 'Average_fold_coverage'
        input_df.rename(columns = {'#Name': 'contigID', 'Coverage': 'Average_fold_coverage'}, inplace = True)
        # Drop columns
        input_df.drop(columns=['Bases', 'RPKM', 'FPKM'], inplace = True)
        # Save library size (extracted from #Reads header of input file)
        with open(in_path+in_file, 'r') as read_file:
            for i, line in enumerate(read_file):
                if '#Reads' in line:
                    libSize = re.sub('[^0-9]','', line)
                    break
        # Calculate mapped read count
        mapped_libSize = input_df['Reads'].sum()
        # Zero out counts below the count_treshold
        input_df['Reads'] = input_df['Reads'].where(input_df['Reads']>=int(args.count_threshold), other=0)
        # Calculate filtered mapped read count
        mapped_libSize_filt = input_df['Reads'].sum()
        # Set lib_norm value based on args.lib_norm setting
        if args.lib_norm == 'mapped':
            lib_norm = mapped_libSize_filt
        elif args.lib_norm == 'total':
            lib_norm = libSize
        else:
            print('Error: --lib_norm setting invalid (must be one of: ["mapped", "total"]). Exiting.')
            exit()
        # Re-calculate RPKM (based on args.lib_norm setting and mapped read threshold filtereing)
        input_df['RPKM'] = round((input_df['Reads']/(int(lib_norm)/1000000))/(input_df['Length']/1000), 4)
        # Re-calculate FPKM (based on args.lib_norm setting and mapped read threshold filtereing)
        input_df['FPKM'] = round(2*((input_df['Frags']/(int(lib_norm)/1000000))/(input_df['Length']/1000)), 4)  
        # Calculate TPM (Note: Due to how TPM is calculated, this is normalised based on *mapped* reads, rather than incorporating total library size, regardless of the setting of args.lib_norm)
        rpk_tmp = input_df['Reads']/(input_df['Length']/1000)
        input_df['TPM'] = rpk_tmp/(rpk_tmp.sum()/1000000)
        # Rename 'Reads' column as 'pileup_counts'
        input_df.rename(columns = {'Reads': 'pileup_counts'}, inplace = True)
        # Add sampleID prefix to column headers (except 'contigID' and 'Length')
        input_df.rename(columns={col: str(sampleID)+'_'+col for col in input_df.columns if col not in ['contigID', 'Length']}, inplace = True)
        # Merge into summary_count_table
        if summary_count_table.empty:
            summary_count_table = input_df.copy()
        else:
            input_df.drop(columns='Length', inplace = True)
            summary_count_table = summary_count_table.merge(input_df, left_on='contigID', right_on='contigID', how = 'outer')
        # Add read counts to read_counts DataFrame
        if read_counts.empty:
            read_counts = pd.DataFrame([[sampleID, libSize, mapped_libSize, mapped_libSize_filt]],
                                       columns=['SampleID', 'Total_library_size', 'Mapped_reads', 'Filtered_mapped_reads'])
        else:
            read_counts = read_counts.append(pd.DataFrame([[sampleID, libSize, mapped_libSize, mapped_libSize_filt]], columns=['SampleID', 'Total_library_size', 'Mapped_reads', 'Filtered_mapped_reads']), ignore_index=True, sort=False)
    # Rearrange summary_count_table column order
    summary_count_table = summary_count_table[['contigID', 'Length'] + \
                                              [col for col in summary_count_table.columns if 'counts' in col] + \
                                              [col for col in summary_count_table.columns if 'coverage' in col] + \
                                              [col for col in summary_count_table.columns if 'RPKM' in col] + \
                                              [col for col in summary_count_table.columns if 'FPKM' in col] + \
                                              [col for col in summary_count_table.columns if 'TPM' in col]
                                             ]
    # Write out summary_count_table to csv
    summary_count_table.to_csv(args.output, sep='\t', index=False)
    # If args.read_counts is set, write out read_counts to csv
    read_counts.to_csv(args.read_counts, sep='\t', index=False)


def norm_comments():
    if args.edger_out is not None and args.sample_map is not None:
        print('NOTE: Read count normalisation (RPKM, FPKM, and TMM) and edgeR differential expression analyses calculated based on ' + args.lib_norm + ' read counts per sample (--libnorm ' + args.lib_norm + ')')
        print('NOTE: Due to how TPM is calculated, this is normalised based on mapped read counts, regardless of the setting of --lib_norm)')
    elif args.sample_map is not None:
        print('NOTE: Read count normalisation (RPKM, FPKM, and TMM) calculated based on based on ' + args.lib_norm + ' read counts per sample (--libnorm ' + args.lib_norm + ')')
        print('NOTE: Due to how TPM is calculated, this is normalised based on mapped read counts, regardless of the setting of --lib_norm)')
    else: 
        print('NOTE: Read count normalisation (RPKM, and FPKM) calculated based onbased on ' + args.lib_norm + ' read counts per sample (--libnorm ' + args.lib_norm + ')')
        print('NOTE: Due to how TPM is calculated, this is normalised based on mapped read counts, regardless of the setting of --lib_norm)')
    print("\n")


def R_edgeR_and_genomeSummary():
    ## NOTE: if runtime becomes an issue, consider updating the to_csv (python) and read_tsv (R) to use feather to convert between the two
    # If --sample_mapping_file and --edger_out provided, calculate TMM (via edgeR in R) and run pairwise differential expression (via edgeR in R)
    # If --genome_mapping_file provided, calculate count summaries to genome level and output genome summary table
    if args.sample_map is not None and args.edger_out is not None:
            print('Running pairwise differential expression analyses (edgeR)')
            print('NOTE: Pairwise differential expression tested via edgeR glmQLFit() and glmQLFTest() (design = model.matrix(~0+Groups))')
            print('NOTE: Pairwise differential expression calculated based on ' + args.lib_norm + ' read counts per sample (--libnorm ' + args.lib_norm + ')')
            print('NOTE: Each pairwise combination tested separately and appended to master summary table (FDR calculated independently for each tested pair)\n')
            if args.genome_map is not None:
                print('Summarising counts to genome level')
                print('NOTE: read count normalisation calculated based on ' + args.lib_norm + ' reads per sample (--libnorm ' + args.lib_norm + ')')
                print('NOTE: read count normalisation calculated based on ' + args.genome_libSize_norm + ' sample read depth (--genome_libSize_norm ' + args.genome_libSize_norm + ')\n')
                subprocess.call(['summarise_counts.R', args.output, args.format, args.sample_map, args.lib_norm, args.read_counts, args.edger_out, args.genome_map, args.genome_libSize_norm], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
            else:
                subprocess.call(['summarise_counts.R', args.output, args.format, args.sample_map, args.lib_norm, args.read_counts, args.edger_out, 'none', args.genome_libSize_norm], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    # If --sample_mapping_file provided, calculate TMM (via edgeR in R)
    # If --genome_mapping_file provided, calculate count summaries to genome level and output genome summary table
    elif args.sample_map is not None:
        if args.genome_map is not None:
            print('Summarising counts to genome level')
            print('NOTE: read count normalisation calculated based on ' + args.lib_norm + ' reads per sample (--libnorm ' + args.lib_norm + ')')
            print('NOTE: read count normalisation calculated based on ' + args.genome_libSize_norm + ' sample read depth (--genome_libSize_norm ' + args.genome_libSize_norm + ')\n')
            subprocess.call(['summarise_counts.R', args.output, args.format, args.sample_map, args.lib_norm, args.read_counts, 'none', args.genome_map, args.genome_libSize_norm], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        else: 
            subprocess.call(['summarise_counts.R', args.output, args.format, args.sample_map, args.lib_norm, args.read_counts, 'none', 'none', args.genome_libSize_norm], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    # If --genome_mapping_file provided, calculate count summaries to genome level and output genome summary table
    elif args.genome_map is not None:
        print('Summarising counts to genome level')
        print('NOTE: read count normalisation calculated based on ' + args.lib_norm + ' reads per sample (--libnorm ' + args.lib_norm + ')')
        print('NOTE: read count normalisation calculated based on ' + args.genome_libSize_norm + ' sample read depth (--genome_libSize_norm ' + args.genome_libSize_norm + ')\n')
        subprocess.call(['summarise_counts.R', args.output, args.format, 'none', args.lib_norm, args.read_counts, 'none', args.genome_map, args.genome_libSize_norm], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


def output_comments():
    print("\nOutput:\n")
    print(args.read_counts + " : Summary table of per-sample read counts (per-sample total library, mapped read, and filtered mapped read counts)\n")
    print(args.output + " : Summary table of per-sample count data for each contig (incl. length, raw read counts (pileup counts or featurecounts), RPKM, FPKM (if --format = pileup), TPM, TMM (if mapping file provided))\n")        
    if args.edger_out is not None and args.sample_map is not None:
        print(args.edger_out + ' : Summary table of edgeR pairwise differential expression analyses\n')   
    if args.genome_map is not None:
        print(os.path.splitext(args.output)[0] + '_genomeSummary.tsv : Summary table of per-sample counts summarised to genome- (or MAG- or bin-) level. Summaries include: raw counts, counts normalised to "genome length" (sum of all contigs in genome); counts normalised to "genome length" and sample read count (mapped reads, or total library size, based on the setting of --lib_norm "mapped")\n')


def main(): 
    print("\n--------------------\n")
    print("Running summarise_counts.py\n\n")
    # Calculate normalisations
    if args.format.lower() == 'featurecounts'.lower():
        norm_from_featurecounts()
    elif args.format.lower() == 'pileup'.lower():
        norm_from_pileup()
    else:
        print('Error: Invalid input file format (-f) specificed. Must be one of: ["pileup", "featurecounts"]. Exiting.')
        exit()
    # summarise_counts.R: If options provided, calculate TMM (edgeR) and/or pairwise comparisons (edgeR) and/or genome-level count summaries
    if args.sample_map is not None:
        print('Calculating normalisation: TMM (edgeR)')
    norm_comments()
    R_edgeR_and_genomeSummary()
    # END
    output_comments()
    print("\nCompleted summarise_counts.py\n")
    print("--------------------\n")


if __name__ == '__main__':
    main()

