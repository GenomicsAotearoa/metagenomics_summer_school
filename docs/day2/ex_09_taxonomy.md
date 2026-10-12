# Assigning taxonomy to refined prokaryotic bins

!!! info "Objectives"

    * [Assigning taxonomy to the refined bins](#assigning-taxonomy-to-the-dereplicated-bins)

---

## Assigning taxonomy to the dereplicated bins

It is always valuable to know the taxonomy of our binned MAGs, so that we can link them to the wider scientific literature. In order to do this, there are a few different options available to us:

1. Extract 16S rRNA gene sequences from the MAGs and classify them
1. Annotate each gene in the MAG and take the consensus taxonomy
1. Use a profiling tool like `Kraken`, which matches pieces of DNA to a reference database using *k*-mer searches
1. Identify a core set of genes in the MAG, and use these to compute a species phylogeny

For this exercise, we will use the last option in the list, making use of the `GTDB-TK` software (available on [github](https://github.com/Ecogenomics/GTDBTk)) to automatically identify a set of highly conserved, single copy marker genes which are diagnostic of the bacterial (120 markers) and archaeal (122 markers) lineages. Briefly, `GTDB-TK` will perform the following steps on a set of bins.

!!! quote ""
    1. Attempt to identify a set of 120 bacterial marker genes, and 122 archaeal marker genes in each MAG.
    1. Based on the recovered numbers, identify which domain is a more likely assignment for each MAG
    1. Create a concatenated alignment of the domain-specific marker genes, spanning approximately 41,000 amino acid positions
    1. Filter the alignment down to approximately 5,000 informative sites
    1. Insert each MAG into a reference tree create from type material and published MAGs
    1. Scale the branch lengths of the resulting tree, as described in [Parks et al.](https://www.ncbi.nlm.nih.gov/pubmed/30148503), to identify an appropriate rank to each branch event in the tree
    1. Calculate ANI and AAI statistics between each MAG and its nearest neighbours in the tree
    1. Report the resulting taxonomic assignment, and gene alignment

This can all be achieved in a single command, although it must be performed through a slurm script due to the high memory requirements of the process.

For the following exercises, we will be working in `8.prokaryotic_taxonomy/`.

Create a new script

!!! terminal-2 "Create script named `gtdbtk.sl`"

    ```bash
    nano gtdbtk.sl
    ```

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      gtdbtk
    #SBATCH --partition     milan
    #SBATCH --time          01:00:00
    #SBATCH --mem           140GB
    #SBATCH --cpus-per-task 24
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Load modules
    module purge
    module load GTDB-Tk/2.4.0-foss-2023a-Python-3.11.6
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.prokaryotic_taxonomy
    
    # Run GTDB-Tk
    gtdbtk classify_wf -x fna --cpus $SLURM_CPUS_PER_TASK \
                       --keep_intermediates \
                       --skip_ani_screen \
                       --genome_dir dastool_bins/ \
                       --out_dir gtdbtk_out/
    ```

!!! terminal-2 "Submit the script"

    ```bash
    sbatch gtdbtk.sl
    ```

As usual, lets look at the parameters here

|<div style="width:200px">Parameter</div>|Description|
|:---|:---|
|`classify_wf`|Specifies the sub-workflow from `GTDB-TK` that we wish to use|
|`-x`|Specify the file extension for MAGs within our input directory.<br>Default is *.fna*, but it's always good practice to specify it anyway|
|`--cpus`|Number of threads/CPUs to use when finding marker genes, and performing tree insertion operations|
|`--keep_intermediates`|Keep intermediate outputs|
|`--genome_dir`|Input directory containing MAGs as individual FASTA files|
|`--out_dir`|Output directory to write the final set of files|

Before submitting your job, think carefully about which set of MAGs you want to classify. You could either use the raw `DAS_Tool` outputs in the `../6.bin_refinement/dastool_out/_DASTool_bins/` folder, the renamed set of bins in the `../6.bin_refinement/example_data_unchopped/` folder (as we have here), the set of curated bins in the `filtered_bins/` folder, or your own set of refined bins. Whichever set you choose, make sure you select the correct input folder and extension setting as it may differ from the example here.

When the task completes, you will have a number of output files provided. The main ones to look for are `gtdbtk.bac120.summary.tsv` and `gtdbtk.arch122.summary.tsv` which report the taxonomies for your MAGs, split at the domain level. These file are only written if MAGs that fall into the domain were found in your data set, so for this exercise we do not expect to see the `gtdbtk.arch122.summary.tsv` file.

If you are interested in performing more detailed phylogenetic analysis of the data, the filtered multiple sequence alignment (MSA) for the data are provided in the `gtdbtk.bac120.msa.fasta` and `gtdbtk.arch122.msa.fasta` files.

Have a look at your resulting taxonomy. The classification of your MAGs will be informative when addressing your research goal for this workshop.

---
