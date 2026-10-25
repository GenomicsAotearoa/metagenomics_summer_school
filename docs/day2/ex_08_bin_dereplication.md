# Binning II: Bin dereplication

!!! info "Objectives"

    - [Bin dereplication using `DAS_Tool` - Creating input tables](#bin-dereplication-using-das_tool---creating-input-tables)
    - [Bin dereplication using `DAS_Tool`](#bin-dereplication-using-das_tool)
    - [Evaluating bins using `CheckM`](#evaluating-bins-using-checkm)
    - [Dereplicating bins across multiple assemblies](#dereplicating-bins-across-multiple-assemblies)

---

## Bin dereplication using `DAS_Tool` - Creating input tables

As we discussed in the previous exercise, we have now generated two sets of bins from the same single assembly. With this mock data set we can see that `MetaBAT` recovered 12 bins, while `MaxBin` and `MetaDecoder` recovered 10. Note that we are aiming to recover prokaryote genomes using these binning tools, and 10 bacterial and archaeal genomes were used in the creation of this mock community. In our case here, these tools have likely recovered 10 bins each of the same genomes. The additional two bins identified by `MetaBAT` may be the result of noise introduced into the binning process by the viral contigs included in the data. Furthermore, it is not clear which tool has done a better job of recruiting contigs to each bin - we very rarely expect to see the complete genome recovered from these kinds of data, so while it is probably the case that an equivalent bin is present in the `MetaBAT`, `MaxBin` and `MetaDecoder` outputs, they will likely be of differing quality.

`DAS_Tool` is a program designed to analyse the bins in each of our binning sets and determine where these equivalent triplets exist and return the 'best' one. `DAS_Tool` does not use the actual bins, but a set of text files that link contigs to their corresponding bins in each of the bin sets. We can produce these files using `bash`.

For this exercise, we will continue working in the `5.binning/` directory.

### Create contig/bin tables - `MetaBAT`

For each of our binning tools, we need to extract the contigs assigned to each bin and create a single file that reports these as

!!! success ""

    ```
    Contig[tab]Bin
    ```

This can be done with a bit of `bash` scripting. There's quite a bit going on here, so we'll provide the full command, and then a step-by-step explanation of what's happening.

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal-2 "Navigate to working directory"

    ```bash
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/
    ```

!!! terminal-2 "Create `metabat_refs.tsv`"

    ```bash linenums='1'
    for bin_path in metabat/*.fa; do
        bin_name=$(basename ${bin_path} .fa)
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> metabat_refs.tsv
    done
    ```

You can check the contents of this file using `less` or `head`:

!!! terminal-2 "Inspect `metabat_association.txt`"

    ```bash
    head -n 5 metabat_refs.tsv
    ```

!!! circle-check "Terminal output"

    ```
    NODE_3_length_835083_cov_0.372452   metabat.10
    NODE_5_length_586102_cov_0.372632   metabat.10
    NODE_20_length_292112_cov_0.373169  metabat.10
    NODE_24_length_265402_cov_0.371105  metabat.10
    NODE_33_length_247431_cov_0.373159  metabat.10
    ```

We will now walk through the content of the "Create `metabat_refs.tsv`" code block above, breaking apart each individual step.


!!! terminal-2 "Line 1: Initiate a `for` loop"

    ```bash
    for bin_path in metabat/*.fa; do
    ```

This is the initial loop, returning each file in the `metabat/` folder that ends with the file extension *.fa*.

!!! terminal-2 "Line 2: Remove path prefix"

    ```bash
        bin_name=$(basename ${bin_path} .fa)
    ```

As we have previously seen, the `basename` command removes the path information from the input variable `bin_path` (i.e. - `metabat/metabat.1.fa` becomes `metabat.1.fa`) and assigns it to a new variable, `bin_name`. As an optional additional parameter, we can pass extra pieces of text to be removed from the variable, in this case the *.fa* extension.

!!! terminal-2 "Line 4: Subset information to create and update `metabat_refs.tsv`"

    ```bash
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> metabat_refs.tsv
    ```

We will break this line down step-by-step:

| <div style="width:230px">Command</div> | Description |
| :--- | :--- |
| `grep ">" ${bin_path}` | The `grep` command searches the input file `bin_path` for lines containing the `>` character, which is the *fastA* demarcation for a sequence name. |
| `sed 's/>//g'` | The first `sed` command replaces the `>` character with the empty character `''`, as we do not need this character in the final file. |
| `sed "s/$/\t${bin_name}/g"` | We now use `sed` again, this time to replace the `$` character. In many command line tools and software environments (including `R` and `python`) the characters `^` and `$` are used as shortcuts for the beginning and ending of a line, respectively. By using the character in this way, we are telling `sed` to replace the end-of-line with the text `\t${bin_name}`. `sed` will parse this text to mean the tab character followed by the content of the `bin_name` variable. The nature of `sed` is that it will automatically insert a new end-of-line. |
| `>> metabat_refs.tsv` | This is similar to the stdout redirection we have previously used, but the double use of the `>` character means that we are appending our text to the end of the file `metabat_refs.tsv`. Because we are looping through several files in this exercise, if we were to use the single `>` character then on each new *fastA* file read, the content of `metabat_refs.tsv` would be replaced.<br><br>It is important to note that because we are **_appending_** to the file, not **_replacing_** the contents, if you make a mistake in the command and need to re-run it, you will need to explicitly delete the `metabat_refs.tsv` file using `rm`, otherwise your new (correct) output will be pasted to the end of your old (incorrect) output.|

### Create contig/bin tables - `MaxBin`, `MetaDecoder`

The process for creating the `MaxBin` and `MetaDecoder` tables is basically the same: Because we have already normalised all bin filename extensions to *.fa*, we only have to change the names for the input directory and the output file.

!!! terminal-2 "Create `maxbin_refs.tsv`"

    ```bash linenums="1"
    for bin_path in maxbin/*.fa; do
        bin_name=$(basename ${bin_path} .fa)
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> maxbin_refs.tsv
    done
    ```

!!! terminal-2 "Create `mdecoder_refs.tsv`"

    ```bash linenums="1"
    for bin_path in metadecoder/*.fa; do
        bin_name=$(basename ${bin_path} .fa)
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> mdecoder_refs.tsv
    done
    ```


!!! warning "Warning for unbinned contigs"

    Both `MetaBAT` and `MaxBin` have the option to output unbinned contigs after binning completes. We have not used that parameter here, but if you do choose to enable it you will end up with another *FASTA* file in your output directory which you will need to avoid in the loops for creating `DAS_Tool` tables.

---

## Bin dereplication using `DAS_Tool`

We are now ready to run `DAS_Tool`. This can be done from the command line, as it does not take a particularly long time to run for this data set.

!!! terminal-2 "Run `DAS_Tool`" 

    ```bash linenums="1"
    # Remove modules to ensure a clean environment
    module purge

    # Load DAS Tool
    module load DAS_Tool/1.1.5-gimkl-2022a-R-4.2.1

    # Create DAS_Tool output directory
    mkdir -p dastool_out/

    # Run DAS_Tool
    DAS_Tool --threads 2  --write_bins \
             --search_engine diamond \
             --bins metabat_refs.tsv,maxbin_refs.tsv,mdecoder_refs.tsv \
             --labels MetaBAT,MaxBin,MetaDecoder \
             --contigs spades_assembly/spades_assembly.m1000.fna \
             --outputbasename dastool_out/
    ```

!!! circle-check "Terminal output"

    ```
    DAS Tool 1.1.5 
    Analyzing assembly 
    Predicting genes 
    Annotating single copy genes using diamond 
    Dereplicating, aggregating, and scoring bins 
    Writing bins
    ```

As usual, we will break down the parameters:

|<div style="width:200px">Parameter</div>|Description|
|:---|:---|
|`--threads`|Number of threads to use|
|`--write_bins`|An argument telling `DAS_Tool` whether or not to write out a new set of bins<br>This is recommended, because `DAS_Tool` can create slices of old bins based on marker composition (see [the paper](https://www.nature.com/articles/s41564-018-0171-1) for details)|
|`--search_engine diamond`|Specify whether to use `usearch`, `diamond`, or `BLAST` as the alignment tool for comparing gene sequences|
|`--bins`|A comma-separated list of the contig/bin files we wish to process|
|`--labels`|A comma-separated list of the binning tools used|
|`--contigs`|Path to the assembly used in binning|
|`--outputbasename`|Output directory for all files|

When `DAS_Tool` has completed, we will have a final set of bins located in the folder path `dastool_out/_DASTool_bins`. Have a look at the output and see which bins made it to the final selection. Did a single binning tool pick the best bins, or are the results a split between `MetaBAT`, `MaxBin` and `MetaDecoder`?

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
    #SBATCH --cpus-per-task 10
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Load modules
    module purge
    module load CheckM/1.2.5-foss-2026-Python-3.14.4
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/
    
    # Run CheckM
    checkm lineage_wf --threads $SLURM_CPUS_PER_TASK \
                      --pplacer_threads $SLURM_CPUS_PER_TASK \
                      --extension fa  --tab_table \
                      --file checkm.tsv \
                      dastool_out/_DASTool_bins/ checkm_out/
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
        #SBATCH --time          00:20:00
        #SBATCH --mem           50GB
        #SBATCH --cpus-per-task 10
        #SBATCH --error         %x_%j.err
        #SBATCH --output        %x_%j.out

        # Load modules
        module purge
        module load CheckM2/1.0.1-Miniconda3

        # Working directory
        cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/

        checkm2 predict --threads $SLURM_CPUS_PER_TASK \
                        --extension .fa \
                        --input dastool_out/_DASTool_bins/ \
                        --output-directory checkm2_out/
        ```

---

## Dereplicating bins across multiple assemblies

In this workshop, we have generated a set of putative MAGs by binning scaffolds taken from a *single co-assembly* of all four samples. Alternatively, we may have chosen to generate multiple assemblies (for example, mini-co-assemblies for each sample group, or individual assemblies for each sample). In this case, it would be necessary to work through the binning process for each assembly, and then conduct an additional dereplication step *across* the multiple assemblies to generate a single set of dereplicated bins for all assemblies. 

This is beyond the scope of this workshop (and unnecessary here, since we are working with a single co-assembly). For future reference for your own work, further information about how to dereplicate bins and viral contigs across multiple assemblies via `dRep` and `dedupe` has been provided as an appendix [here](../resources/1_APPENDIX_ex8_Dereplication.md).

---
