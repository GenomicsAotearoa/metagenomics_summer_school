# Binning with multiple tools

!!! info "Objectives"

    * [Overview](#overview)
    * [Create initial bins using `MetaBAT`](#metabat)
    * [Create initial bins using `MaxBin`](#maxbin)

---

## Overview

With the mapping information computed in the last exercise, we can now perform binning. There are a multitude of good binning tools currently published, and each have their strengths and weaknesses. As there is no best tool for binning, the current strategy for binning is to use a number of different tools on your data, then use the tool `DAS_Tool` to evaluate all potential outcomes and define the best set of bins across all tools used.

In our own workflow, we use the tools `MetaBAT`, `MaxBin`, and `CONCOCT` for binning, but there are many alternatives that are equally viable. In the interests of time, we are only going to demonstrate the first two tools. However, we recommend that you experiment with some of the following tools when conducting your own research.

1. [Tetra-ESOM](https://github.com/tetramerFreqs/Binning)
1. [VAMB](https://github.com/RasmussenLab/vamb)

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
    module load MetaBAT/2.15-GCC-11.3.0

    # Manual specification of files
    jgi_summarize_bam_contig_depths --outputDepth metabat.txt sample1.bam sample2.bam sample3.bam sample4.bam

    # Wildcard
    jgi_summarize_bam_contig_depths --outputDepth metabat.txt sample?.bam
    ```

??? tip "Bash wildcards"

    The question mark `?` is also a Bash wildcard. The difference between `?` and `*` is that:

    - `?` represents any character ***once***
    - `*` represents any character ***zero or more times***

    The reason we use a `?` here instead is because if you did the SPAdes contig mapping with the `for` loop AND with the array, you will capture each sample twice (i.e., `sample1.bam` and `sample1_a.bam`) with the asterisk `*`. For this exercise, we really want each sample file to be counted only once.

Both give the same result, although the sample order may vary.

We can then pass the table `metabat.txt` into the `MetaBAT` binning tool.

Before we proceed, note that when you run `MetaBAT` on NeSI you will see the text `vGIT-NOTFOUND` appear in your command line. This has no impact on the performance of the tool.

!!! terminal-2 "Run MetaBAT"

    ```bash
    metabat2 -t 2 -m 1500 \
             -i spades_assembly/spades_assembly.m1000.fna \
             -a metabat.txt \
             -o metabat/metabat
    ```

Note here that we are specifying a minimum contig size of 1,500 bp, which is the lower limit allowed by the authors of `MetaBAT`. This is larger than the minimum threshold of 1,000 bp when we filtered the assembly, which means there are some assembled contigs which cannot be binned. Consider the choice of this parameter and your initial contig filtering carefully when binning your own data.

When specifying the output file, notice that we pass both a folder path (*metabat/*) and file name (*metabat*). The reason I do this is that `MetaBAT` writes its output files using the pattern

!!! success ""

    `[USER VALUE].[BIN NUMBER].fa`

If we only provided the path, without a file name prefix, `MetaBAT` would create output like the following:

!!! success "" 

    ```bash
    metabat/.1.fa
    metabat/.2.fa
    metabat/.3.fa
    ```

The problem with this is that on Linux systems, prefixing a file or folder name with a '.' character means the the file is hidden. This can lead to a lot of confusion when your binning job completes successfully but no files are visible!

---

## `MaxBin`

Like `MetaBAT`, `MaxBin` requires a text representation of the coverage information for binning. Luckily, we can be sneaky here and just reformat the `metabat.txt` file into the format expected by `MaxBin`. We use `cut` to select only the columns of interest, which are the *contigName* and coverage columns, but not the *contigLen*, *totalAvgDepth*, or variance columns.

We can inspect the `metabat.txt` file with `head` or `less` to identify the correct column indices for `cut`.

!!! terminal "code"

    ```bash
    less metabat.txt
    ```

!!! circle-check "Terminal output: Snapshot of `metabat.txt`"

    ```
    contigName      contigLen       totalAvgDepth   sample1.bam     sample1.bam-var sample2.bam     sample2.bam-var sample3.bam     sample3.bam-var sample4.bam     sample4.bam-var
    NODE_1_length_1221431_cov_0.752208      1.22143e+06     18.9178 3.55592 3.51379 10.865  11.3668 4.49688 4.51001 0       0
    NODE_2_length_871377_cov_1.172283       871377  29.4909 6.14223 6.1249  2.81105 2.77014 20.5376 20.6339 0       0
    NODE_3_length_835083_cov_0.372452       835083  9.2898  1.79405 1.74475 4.2584  4.19122 3.23735 3.23068 0       0
    NODE_4_length_803085_cov_0.518929       803085  13.054  2.63945 2.65518 2.53475 2.48676 7.87984 8.38717 0       0
    NODE_5_length_686960_cov_0.527949       686960  13.2164 2.64157 2.98377 2.59676 2.90338 7.97809 10.8339 0       0
    NODE_6_length_607162_cov_1.000759       607162  25.0853 16.7701 17.576  2.759   2.74899 5.55618 5.69041 0       0
    ```

!!! terminal-2 "Create `maxbin.txt` based on `metabat.txt`"

    ```bash
    cut -f1,4,6,8,10 metabat.txt > maxbin.txt
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
                  -abund maxbin.txt \
                  -out maxbin/maxbin
    ```

!!! terminal-2 "Submit the script as a slurm job"

    ```bash
    sbatch maxbin_clustering.sl
    ```

!!! note "`MaxBin` runtime"

    This will take a bit longer to complete, as `MaxBin` uses gene prediction tools to identify the ideal contigs to use as the start of each bin.

---
