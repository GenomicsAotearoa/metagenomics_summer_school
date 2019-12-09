# Introduction to binning

### Objectives

* Overview
* Create initial bins using `MetaBAT`
* Create initial bins using `MaxBin`

---

### Overview

With the mapping information computed in the last exercise, we can now perform binning. There are a multitude of good binning tools currently published, and each have their strengths and weaknesses. As there is no best tool for binning, the current strategy for binning is to use a number of different tools on your data, then use the tool `DAS_Tool` to evaluate all potential outcomes and define the best set of bins across all tools used.

In our own workflow, we use the tools `MetaBAT`, `MaxBin`, and `CONCOCT` for binning, but there are many alternatives that are equally viable. In the interests of time, we are only going to demonstrate the first two tools. However, we recommend that you experiement with some of the following tools when conducting your own research.

1. [GroopM](http://ecogenomics.github.io/GroopM/)
1. [Tetra-ESOM](https://github.com/tetramerFreqs/Binning)
1. [VAMB](https://github.com/RasmussenLab/vamb)

---

### *MetaBAT*

`MetaBAT` binning occurs in two steps. First, the *bam* files from the last exercise are parsed into a tab-delimited table of the average coverage depth and variance per sample mapped. Binning is then performed using this table.

The *.bam* files can be passed in in either a user-defined order, or using wildcards.
```bash
module load MetaBAT/2.13-GCC-7.4.0

# Manual specification of files
jgi_summarize_bam_contig_depths --outputDepth metabat.txt sample1.bam sample2.bam sample3.bam sample4.bam

# Wildcard
jgi_summarize_bam_contig_depths --outputDepth metabat.txt sample*.bam
```

Both give the same result, although the sample order may vary.

We can then pass the table `metabat.txt` into the `MetaBAT` binning tool.

Before we proceed, note that when you run `MetaBAT` on NeSI you will see the text `vGIT-NOTFOUND` appear in your command line. This has no impact on the performance of the tool.

```bash
metabat2 -t 10 -m 1500 \
         -i spades_assembly/spades_assembly.m1000.fna \
         -a metabat.txt \
         -o metabat/metabat
```

Note here that we are specifying a minimum contig size of 1,500 bp, which is the lower limit allowed by the authors of `MetaBAT`. This is larger than the minimum threshold of 1,000 bp when we filtered the assembly, which means there are some assembled contigs which cannot be binned. Consider the choice of this parameter and your initial contig filtering carefully when binning your own data.

When specifying the output file, notice that we pass both a folder path (*metabat/*) and file name (*metabat*). The reason I do this is that `MetaBAT` writes its output files using the pattern

`[USER VALUE].[BIN NUMBER].fa`

If we only provided the path, without a file name prefix, `MetaBAT` would create output like the following:

```bash
metabat/.1.fa
metabat/.2.fa
metabat/.3.fa
```

The problem with this is that on Linux systems, prefixing a file or folder name with a '.' character means the the file is hidden. This can lead to a lot of confusion when your binning job completes successfully but no files are visible!

---

### *MaxBin*

Like `MetaBAT`, `MaxBin` requires a text representation of the coverage information for binning. Luckily, we can be sneaky here and just reformat the `metabat.txt` file into the format expected by `MaxBin`. We use `cut` to select only the columns of interest, which are the *contigName* and coverage columns, but not the *contigLen*, *totalAvgDepth*, or variance columns.

We can inspect the `metabat.txt` file with `head` or `less` to identify the correct column indices for `cut`.

```bash
less metabat.txt

cut -f1,4,6,8,10 metabat.txt > maxbin.txt
```

Generally speaking, this pattern of first, fourth, then [n + 2]<sup>th</sup> should work for any number of mapping files, although we always reocmmend that you check and confirm before you continue.

This table is then passed to `MaxBin`. Unlike the case with `MetaBAT`, if we want to direct the output files into a folder, we must create that folder in advance.

```bash
module load MaxBin/2.2.6-gimkl-2018b-Perl-5.28.1

mkdir -p maxbin/
run_MaxBin.pl -thread 10 -min_contig_length 1500 \
              -contig spades_assembly/spades_assembly.m1000.fna \
              -abund maxbin.txt \
              -out maxbin/maxbin
```

This will take a bit longer to complete, as `MaxBin` uses gene prediction tools to identify the ideal contigs to use as the start of each bin.

---