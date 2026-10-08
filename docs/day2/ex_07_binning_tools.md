# Binning I: Bin with multiple tools

!!! info "Objectives"

    * [Overview](#overview)
    * [Create initial bins using `MetaBAT`](#metabat)
    * [Create initial bins using `MaxBin`](#maxbin)
    * [Create initial bins using `MetaDecoder`](#metadecoder)

---

## Overview

With the mapping information computed in the last exercise, we can now perform binning. There are a multitude of good binning tools currently published, and each have their strengths and weaknesses. As there is no best tool for binning, the current strategy for binning is to use a number of different tools on your data, then use the tool `DAS_Tool` to evaluate all potential outcomes and define the best set of bins across all tools used.

In our own workflow, we often use the tools `MetaBAT` and `MaxBin` for binning, but there are many alternatives that are equally viable. In the interests of time, we are only going to demonstrate one alternative, `MetaDecoder`. We recommend that you experiment with several tools when conducting your own research:

1. [MetaBAT](https://bitbucket.org/berkeleylab/metabat)
1. [MaxBin](https://sourceforge.net/projects/maxbin)
1. [MetaDecoder](https://github.com/liu-congcong/MetaDecoder)
1. [COMEBin](https://github.com/ziyewang/COMEBin)
1. ...

---

## `MetaBAT`

`MetaBAT` binning occurs in two steps. First, the *bam* files from the last exercise are parsed into a tab-delimited table of the average coverage depth and variance per sample mapped. Binning is then performed using this table.

The *.bam* files can be passed in via either a user-defined order, or using wildcards.

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal-2 "Navigate to working directory"

    ```bash
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/
    ```

!!! terminal "code"

    ```bash
    module purge
    module load MetaBAT/2.17-GCC-12.3.0

    # Manual specification of files
    jgi_summarize_bam_contig_depths --outputDepth metabat.tsv sample1.bam sample2.bam sample3.bam sample4.bam

    # Wildcard
    jgi_summarize_bam_contig_depths --outputDepth metabat.tsv sample?.bam
    ```

??? tip "Bash wildcards"

    The question mark `?` is also a Bash wildcard. The difference between `?` and `*` is that:

    - `?` represents any character ***once***
    - `*` represents any character ***zero or more times***

    The reason we use a `?` here instead is because if you did the SPAdes contig mapping with the `for` loop AND with the array, you will capture each sample twice (i.e., `sample1.bam` and `sample1_a.bam`) with the asterisk `*`. For this exercise, we really want each sample file to be counted only once.

Both give the same result, although the sample order may vary.

??? tip ".tsv files"

    Tab separated value files with the extension `.tsv` will be displayed in tabular layout by Jupyter when opening the file from the Jupyter file browser.


We can then pass the table `metabat.tsv` into the `MetaBAT` binning tool.

Before we proceed, note that when you run `MetaBAT` on NeSI you will see the text `GIT-NOTFOUND` appear in your command line. This has no impact on the performance of the tool.

!!! terminal-2 "Run MetaBAT"

    ```bash
    metabat2 --numThreads 2 --seed 13 \
             --minContig 1500 \
             --inFile spades_assembly/spades_assembly.m1000.fna \
             --abdFile metabat.tsv \
             --outFile metabat/metabat
    ```

Note that many metagenomic software tools utilise stochastic processes and therefore make use of random numbers at some point. As a consequence, the results are usually not fully reproducible, i.e., different runs may result in slightly different bins. Often the software tools provide an option to supply a seed for the random number generator, which would be a machine generated random number otherwise. If you want to receive the same deterministic result with every run you will need to provide the same random number seed. For clustering with `MetaBAT2` this can be achieved with the `--seed` option.

Note also that we are specifying a minimum contig size of 1,500 bp, which is the lower limit allowed by the authors of `MetaBAT`. This is larger than the minimum threshold of 1,000 bp when we filtered the assembly, which means there are some assembled contigs which cannot be binned. Consider the choice of this parameter and your initial contig filtering carefully when binning your own data.

When specifying the output file, notice that we pass both a folder path (*metabat/*) and file name (*metabat*). The reason I do this is that `MetaBAT` writes its output files using the pattern

!!! success ""

    `[USER VALUE].[BIN NUMBER].fa`

If we only provided the path, without a file name prefix, `MetaBAT` would create output like the following:

!!! success "" 

    ```
    metabat/.1.fa
    metabat/.2.fa
    metabat/.3.fa
    ```

The problem with this is that on Linux systems, prefixing a file or folder name with a '.' character means the the file is hidden. This can lead to a lot of confusion when your binning job completes successfully but no files are visible!

---

## `MaxBin`

Like `MetaBAT`, `MaxBin` requires a text representation of the coverage information for binning. Luckily, we can be sneaky here and just reformat the `metabat.tsv` file into the format expected by `MaxBin`. We use `cut` to select only the columns of interest, which are the *contigName* and coverage columns, but not the *contigLen*, *totalAvgDepth*, or variance columns.

We can inspect the `metabat.tsv` file with `head` or `less` to identify the correct column indices for `cut`.

!!! terminal "code"

    ```bash
    less metabat.tsv
    ```

!!! circle-check "Terminal output: Snapshot of `metabat.tsv`"

    ```
    contigName      contigLen       totalAvgDepth   sample1.bam     sample1.bam-var sample2.bam     sample2.bam-var sample3.bam     sample3.bam-var sample4.bam     sample4.bam-var
    NODE_1_length_1221431_cov_0.752208      1.22143e+06     18.9178 3.55592 3.51379 10.865  11.3668 4.49688 4.51001 0       0
    NODE_2_length_871377_cov_1.172283       871377  29.4909 6.14223 6.1249  2.81105 2.77014 20.5376 20.6339 0       0
    NODE_3_length_835083_cov_0.372452       835083  9.2898  1.79405 1.74475 4.2584  4.19122 3.23735 3.23068 0       0
    NODE_4_length_803085_cov_0.518929       803085  13.054  2.63945 2.65518 2.53475 2.48676 7.87984 8.38717 0       0
    NODE_5_length_686960_cov_0.527949       686960  13.2164 2.64157 2.98377 2.59676 2.90338 7.97809 10.8339 0       0
    NODE_6_length_607162_cov_1.000759       607162  25.0853 16.7701 17.576  2.759   2.74899 5.55618 5.69041 0       0
    ```

!!! terminal-2 "Create `maxbin.tsv` based on `metabat.tsv`"

    ```bash
    cut -f1,4,6,8,10 metabat.tsv > maxbin.tsv
    ```

Generally speaking, this pattern of first, fourth, then [n + 2]<sup>th</sup> should work for any number of mapping files, although we always recommend that you check and confirm before you continue.

This table is then passed to `MaxBin`. Unlike the case with `MetaBAT`, if we want to direct the output files into a folder, we must create that folder in advance.

!!! terminal-2 "Create a new script to submit as a slurm job"

    ```bash
    nano maxbin_clustering.sl
    ```

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      maxbin_clustering
    #SBATCH --time          00:05:00
    #SBATCH --partition     milan
    #SBATCH --mem           10GB
    #SBATCH --cpus-per-task 10
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Load modules
    module purge
    module load MaxBin/2.2.7-GCC-11.3.0-Perl-5.34.1
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/
    
    # Output directory
    mkdir -p maxbin/
    
    # Run MaxBin
    run_MaxBin.pl -thread $SLURM_CPUS_PER_TASK -min_contig_length 1500 \
                  -contig spades_assembly/spades_assembly.m1000.fna \
                  -abund maxbin.tsv \
                  -out maxbin/maxbin
    
    # Normalise filenames
    rename ".fasta" ".fa" maxbin/*.fasta
    ```

!!! terminal-2 "Submit the script as a slurm job"

    ```bash
    sbatch maxbin_clustering.sl
    ```

!!! note "`MaxBin` runtime"

    This will take a bit longer to complete, as `MaxBin` uses gene prediction tools to identify the ideal contigs to use as the start of each bin.


---

## `MetaDecoder`

`MetaDecoder` binning occurs in three steps. The `MetaDecoder` software package contains all required tools. First, as with `MetaBAT`, the *bam* files are parsed into a tab-delimited table of average coverage depths. Then single-copy marker genes are mapped to the assembly to generate seed sequences for the clustering bins. Finally, clustering is performed using the coverage table, the *seed* file, and internally calculated compositions (k-mer frequencies).

Note that in order to get fully reproducible results with `MetaDecoder` the `--random_number` parameter should be set (cf. `--seed` for MetaBAT2).

Note also that we are specifying a minimum contig size of 2,000 bp, which is the lower limit allowed for `MetaDecoder`, the default value would have been 2,500 bp.

!!! terminal-2 "Create a new script to submit as a slurm job"

    ```bash
    nano metadecoder_clustering.sl
    ```


!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      metadecoder_clustering
    #SBATCH --time          00:05:00
    #SBATCH --partition     milan
    #SBATCH --mem           2GB
    #SBATCH --cpus-per-task 4
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Load modules and SW environment
    module purge
    module load Python/3.11.6-foss-2023a 
    source /nesi/project/nesi02659/MGSS_2026/resources/tools/MetaDecoder_v1.2.2_python3.11-venv/bin/activate
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/
    
    # I/O
    src="spades_assembly/spades_assembly.m1000.fna"
    dstdir="metadecoder"

    # Create destination dir
    mkdir -p "$dstdir"

    # Run MetaDecoder pipeline
    metadecoder coverage --threads $SLURM_CPUS_PER_TASK \
            --bam *.bam \
            --output "mdecoder.tsv"

    metadecoder seed --threads $SLURM_CPUS_PER_TASK \
            --fasta "$src" \
            --output "${dstdir}/mdecoder.seed"

    metadecoder cluster --threads $SLURM_CPUS_PER_TASK \
            --min_sequence_length 2000 --random_number 13 \
            --fasta "$src" \
            --coverage "mdecoder.tsv" \
            --seed "${dstdir}/mdecoder.seed" \
            --output "${dstdir}/metadecoder"

    # Normalise filenames and finish
    rename ".fasta" ".fa" "$dstdir"/*.fasta
    deactivate 
    ```

!!! terminal-2 "Submit the script as a slurm job"

    ```bash
    sbatch metadecoder_clustering.sl
    ```



---
