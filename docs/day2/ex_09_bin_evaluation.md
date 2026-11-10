# Evaluating prokaryotic bins using CheckM

!!! info "Objectives"

    - [Evaluating bins using `CheckM`](#evaluating-bins-using-checkm)

---

## Evaluating bins using `CheckM`

Now that we have our dereplicated set of bins, it is a good idea to determine estimates of their completeness (how much of the genome was recovered) and contamination (how many contigs we believe have been incorrectly assigned to the bin). For organisms that lack a reference genome there is no **_definitive_** way to do this, but the tool `CheckM` provides a robust estimate for these statistics by searching each of your bins for a number of highly conserved, single copy genes. The number of markers depends on whether or not you are working with bacterial (120 markers) or archaeal (122 markers) genomes, but `CheckM` is able to determine which set is more appropriate for each of your bins as it runs.

There are several characteristics of the `CheckM` marker set worth noting:

**Highly conserved, single copy markers**

The marker sets used in `CheckM` were chosen because they are present in at least 95% of bacterial/archaeal genomes, and are single copy in ≥97% genomes tested. This means that if a gene is missing from a genome, it is likely due to incompleteness in either the original assembly or the binning approach. Similarly, if a marker is observed more than once in a bin it is likely the result of over-clustering of the data.

**Genes are considered as co-located clusters**

Rather than test the raw presence/absence of genes in the marker sets, the genes are organised into operon-like structures where genes known to be co-located are placed together. This is advantageous for two reasons

1. These co-located groups are distributed around the prokaryotic genome, so estimates are not biased by lucky/unlucky recovery of a gene hotspot
1. `CheckM` can account for how complete each individual gene cluster is, rather than just whether or not genes are present

**Lineage-specific duplications and losses can be identified**

As part of determining the correct marker set to use for each bin (bacterial or archaeal), `CheckM` uses a set of 43 conserved prokaryotic markers to insert each bin into a guide tree to estimate the phylogeny of the bin. There are several lineages which are known to have lost particular markers, or to have acquired a additional copies, and if `CheckM` places a bin into one of these lineages it can adjust its completeness/contamination estimates accordingly.

For the following exercises, we will be working in `5.bin_quality_and_taxonomy/`. Our results from the previous section (`../4.binning/dastool_out/_DASTool_bins/`) have been copied to `5.bin_quality_and_taxonomy/dastool_bins/`.

We will need to run `CheckM` under a slurm script. This is because the tree placement process requires a large amount of memory to perform, independently of the size of your data set. A basic script for submitting a `CheckM` job would be as follows:

Create a new script

!!! terminal-2 "Create script named `checkm.sl`"
    
    ```bash
    nano checkm.sl
    ```

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder."

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      CheckM
    #SBATCH --partition     milan
    #SBATCH --time          00:20:00
    #SBATCH --mem           40GB
    #SBATCH --cpus-per-task 8
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Load modules
    module purge
    module load CheckM/1.2.5-foss-2026-Python-3.14.4
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.bin_quality_and_taxonomy/
    
    # Run CheckM
    checkm lineage_wf --threads $SLURM_CPUS_PER_TASK \
                      --pplacer_threads $SLURM_CPUS_PER_TASK \
                      --extension fa  --tab_table \
                      --file checkm.tsv \
                      dastool_bins/ checkm_out/
    ```

!!! terminal-2 "Submit the script as a slurm job"

    ```bash
    sbatch checkm.sl
    ```

The breakdown of parameters is as follows

|<div style="width:200px">Parameter</div>|Function|
|:---|:---|
|`lineage_wf`|Specifies with mode of `CheckM` to run. This is the most common to use, but several others exist|
|`--threads`|Number of threads to use for the initial marker gene detection and clustering|
|`--pplacer_threads`|Number of threads to use when inserting bins into the guide tree|
|`--extension`|The *FASTA* file extension to look for in the input folder. Default is *.fna*|
|`--tab_table`|If this parameter is present, a summary table will be written for the `CheckM` run|
|`--file`|The name of the output file for the summary|
|`dastool/_DASTool_bins/`|The location of the bins to test|
|`checkm_out/`|The location to write intermediate and output files|

!!! tip "Increasing `--pplacer_threads`"

    Increasing this parameter results in a *linear* increase in memory requirement - seting it to 10 means that `CheckM` will need about 10 times more memory than with a single thread

When your job completes, we will download the summary file and examine it.

!!! note "CheckM2"

    An independent update to CheckM was made in late 2022 in the form of a new software [CheckM2](https://github.com/chklovski/CheckM2) (see [publication](https://doi.org/10.1038/s41592-023-01940-w)). It uses machine learning to estimate completeness and contamination. Computationally, it uses less resources than CheckM and can be used on putative genomes with reduced genomes or unusual biology. Keep in mind that because both software use different methods for assessing genome completeness, the output statistics will be different (see [discussion here for a case example](https://github.com/chklovski/CheckM2/issues/32)) and it does not estimate strain heterogeneity.
    
    If you would like to try it out, here's a script analogous to that we ran using CheckM:

    !!! terminal "code"

        ```bash linenums="1"
        #!/bin/bash -e
        #SBATCH --account       nesi02659
        #SBATCH --job-name      CheckM2
        #SBATCH --partition     milan
        #SBATCH --time          00:10:00
        #SBATCH --mem           10GB
        #SBATCH --cpus-per-task 8
        #SBATCH --error         %x_%j.err
        #SBATCH --output        %x_%j.out

        # Load modules
        module purge
        module load CheckM2/1.0.1-Miniconda3

        # Working directory
        cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.bin_quality_and_taxonomy/

        checkm2 predict --threads $SLURM_CPUS_PER_TASK \
                        --extension .fa \
                        --input dastool_bins/ \
                        --output-directory checkm2_out/
        ```
