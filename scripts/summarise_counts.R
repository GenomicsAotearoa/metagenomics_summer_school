#!/usr/bin/env Rscript

# NOTE: this Rscript takes input variables passed to it via the python script summarise_counts.py
# Input variable (args) order: [args.output, args.format, args.sample_map, args.lib_norm, args.read_counts, args.edger_out, args.genome_map, args.genome_libSize_norm]

## Load arguments
args = commandArgs(trailingOnly=TRUE)

## Load packages (install and load if not already installed)
# Main packages
dynamic_require <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE , repos = "http://cran.us.r-project.org", quiet = TRUE)
      require( i , character.only = TRUE )
    }
  }
}
dynamic_require( c("dplyr" , "tibble" , "readr", "tidyr", "fuzzyjoin", "stringr", "matrixStats") )
# edgeR
if(!require(edgeR)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("edgeR")
    library(edgeR)
}
# EDAseq
if(!require(EDASeq)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("EDASeq")
    library(EDASeq)
}

## import read counts summary table
summary_count_table_R <- read_tsv(args[1])

## If --sample_mapping_file provided, calculate TMM normalisation
## IF --edgeR_out also provided, also run edgeR pairwise differential expression analyses 
if (args[3] != 'none') {
    ### Prepare data for edgeR
    # Generate empty list for edgeR calculations
    edgeR_data = list()
    # Read in raw counts (counts) from summary count table
    if (args[2] == 'pileup') {
        edgeR_data$Counts = summary_count_table_R %>%
            select(contigID, contains("counts")) %>%
            column_to_rownames(var="contigID")
    } else if (args[2] == 'featurecounts') {
        edgeR_data$Counts = summary_count_table_R %>%
            select(geneID, contains("counts")) %>%
            column_to_rownames(var="geneID")
    }
    # Extract groups from sample mapping file: Read in sample mapping file, match sample mapping file sampleIDs with 'Reads' sampleIDs in counts table above, extract groups in same order as 'Reads' sampleIDs
    mapping_file <- read_tsv(args[3])
    Groups <- data.frame(countsID = names(edgeR_data$Counts)) %>% 
        fuzzy_left_join(mapping_file, by = c("countsID" = names(mapping_file[1])), match_fun = str_detect) %>%
        pull(group) %>%
        as.factor()
    # Build DGE object
    edgeR_data$DGE = DGEList(counts=edgeR_data$Counts, group=Groups)
    # Update library size with total library size if set in python call (--lib_norm total) (otherwise leave as defulat (mapped reads))
    if (args[4] == 'total') {
        edgeR_data$DGE$samples$lib.size <- read_tsv(args[5]) %>%
            pull(Total_library_size)
    }
    # Calculate normalisation factors
    edgeR_data$DGE <- calcNormFactors(edgeR_data$DGE, method="TMM")
    ### Calculate TMM normalised counts and write out updated summary count table
    # Calculate normalised counts
    counts_df_TMM <- cpm(edgeR_data$DGE)
    if (args[2] == 'pileup') {
        counts_df_TMM <- as.data.frame(counts_df_TMM) %>%
            setNames(gsub("_pileup_counts","_TMM",names(.))) %>%
            rownames_to_column(var = "contigID") 
    } else if (args[2] == 'featurecounts') {
        counts_df_TMM <- as.data.frame(counts_df_TMM) %>%
            setNames(gsub("_featurecounts_counts","_TMM",names(.))) %>%
            rownames_to_column(var = "geneID")     
    }
    # Join with summary_count_table_R 
    if (args[2] == 'pileup') {
        summary_count_table_R <- summary_count_table_R %>%
            left_join(counts_df_TMM, by = "contigID")
    } else if (args[2] == 'featurecounts') {
        summary_count_table_R <- summary_count_table_R %>%
            left_join(counts_df_TMM, by = "geneID")
    }
    # Write out
    write_tsv(summary_count_table_R, args[1])
    ### If --edger_out provided to summarise_counts.py, run edgeR pairwise differential expression analyses and write out edgeR analyses summary table
    if (args[6] != 'none') {
        # Generate design
        # n.b. The design without intercept (~0+groups) does not use the first sample as the reference, therefore you can make all possible combinations via the contrast function
        edgeR_data$design = model.matrix(~0+Groups)
        # Estimate dispersion
        edgeR_data$DGE = estimateDisp(edgeR_data$DGE, edgeR_data$design)
        # GLM
        edgeR_data$DGE_fit = glmQLFit(edgeR_data$DGE, edgeR_data$design)
        # Generate all pairwise comparison options from sample mapping file groups
        pairwise_groups_list = combn(unique(Groups), 2, simplify = TRUE)
        # Create list of contrasts to pass as parameters to makeContrasts via do.call
        pairwise_contrasts <- list()
        for (i in 1:ncol(pairwise_groups_list)) {
            pairwise_contrasts[paste0(pairwise_groups_list[,i][1], '_vs_', pairwise_groups_list[,i][2])] <- paste0('Groups', pairwise_groups_list[,i][1], '-Groups', pairwise_groups_list[,i][2])
        }
        # Add levels=edgeR_data$design to contrasts list
        pairwise_contrasts <- append(pairwise_contrasts, list(levels = edgeR_data$design))
        # Create contrast object of all pairwise comparisons options
        edgeR_data$group_contrasts <- do.call(makeContrasts, pairwise_contrasts)
        # Create output for glmQLFTest summary
        edgeR_data$glmQLFTest_summary <- NULL
        # Test each contrast
        for (i in 1:ncol(edgeR_data$group_contrasts)){
            # run test
            tmp.glmQLFTest = glmQLFTest(edgeR_data$DGE_fit, contrast = edgeR_data$group_contrasts[,i])
            # Modify stats output table
            tmp.glmQLFTest$table = tmp.glmQLFTest$table %>% 
                rownames_to_column('Gene') %>%
                mutate(Contig = gsub('_\\d+$', '', Gene) ) %>%
                mutate(FDR = p.adjust(PValue, method='BH')) %>%
                mutate(Comparison = dimnames(edgeR_data$group_contrasts)$Contrasts[i]) %>%
                arrange(PValue)
            # Bind to summary table
            edgeR_data$glmQLFTest_summary <- bind_rows(edgeR_data$glmQLFTest_summary, tmp.glmQLFTest$table) 
        }
    # Write out edgeR glmQLFTest summary table
    write_tsv(edgeR_data$glmQLFTest_summary, args[6])
    }
}

## If --genome_contig_mapping_file provided to summarise_counts.py, also summarise at the genome level, incl: summed lengths and read counts; summarised counts normalised by genome length; summarised counts normalised by genome length and library size (mapped reads or total library size, depending on --lib_norm)
if (args[7] != 'none') {
    # Read in genome-to-contig mapping file/lookup table
    genome2contig_df <- read_tsv(args[7])
    # Check binID, magID, or genomeID column present (update to genomeID header for downstream use)
    if ("binID" %in% colnames(genome2contig_df)) {
        genome2contig_df <- genome2contig_df %>%
            dplyr::rename(genomeID = binID)
    }
    if ("magID" %in% colnames(genome2contig_df)) {
        genome2contig_df <- genome2contig_df %>%
            dplyr::rename(genomeID = magID)
    } 
    # Check that genomeID and contigID are both labelled correctly
    if((!"genomeID" %in% colnames(genome2contig_df)) || (!"contigID" %in% colnames(genome2contig_df))) {
        print('Error: Invalid column names provided to --genome_contig_map. Columns must be: (genomeID or magID or binID) and contigID')
    }
    # Join counts with binID lookup table
    count_table_wBins <- summary_count_table_R %>%
        select(matches("geneID"), contigID, Length, contains("counts")) %>%
        left_join(genome2contig_df)
    # group_by genomeID and summarise length and per-sample counts
    count_table_wBins_summarised <- count_table_wBins %>%
        group_by(genomeID) %>%
        summarise(genome_length = sum(Length), across(contains("counts"), ~ sum(.x)))
    # normalise to 'genome length' (per kilobase, or per base?)
    count_table_wBins_summarised <- count_table_wBins_summarised %>%
        mutate(across(contains("counts"),  ~ round(.x / (genome_length/1000), 4), .names = "{col}_lengthNorm"))
    # normalise to mean library size
    read_counts_df <- read_tsv(args[5])
    if (args[4] == 'total') {
        if (args[8] == 'min') {
            libSize_multiplier <- min(read_counts_df[['Total_library_size']])
        } else if (args[8] == 'mean') {
            libSize_multiplier <- mean(read_counts_df[['Total_library_size']])
        } else {
            print('Error: Invalid option passed to --genome_libSize_norm. Valid options = ["min", "mean"]. Using the default setting (min)')
            libSize_multiplier <- min(read_counts_df[['Total_library_size']])
        } 
    } else if (args[4] == 'mapped') {
        if (args[8] == 'min') {
            libSize_multiplier <- min(read_counts_df[['Filtered_mapped_reads']])
        } else if (args[8] == 'mean') {
            libSize_multiplier <- mean(read_counts_df[['Filtered_mapped_reads']])
        } else {
            print('Error: Invalid option passed to --genome_libSize_norm. Valid options = ["min", "mean"]. Using the default setting (min)')
            libSize_multiplier <- min(read_counts_df[['Filtered_mapped_reads']])
        } 
    }
    for (sampleID in read_counts_df[['SampleID']]) {
        sample_libSize <- read_counts_df %>% 
            filter(SampleID == sampleID) %>%
            pull(Total_library_size)
        count_table_wBins_summarised <- count_table_wBins_summarised %>%
            mutate(across((matches(sampleID) & matches('_lengthNorm')), ~ round(((.x / sample_libSize) * libSize_multiplier), 2), .names = "{col}_libNorm"))
    }
    # Write out count table
    write_tsv(count_table_wBins_summarised, paste0(tools::file_path_sans_ext(args[1]), '_genomeSummary.tsv'))
}

