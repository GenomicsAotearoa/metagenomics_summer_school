# Assembly

### Objectives

* Become familiar with the standard input files for **SPAdes** and **IDBA-UD**
* Understand the basic parameters that should be modified when using these assemblers
* Prepare an assembly job to run under slurm

All work for this exercise will occur in the *3.assembly/* directory.

---

### The standard input files for *SPAdes* and *IDBA-UD*

Although they both make use of the same types of data, both **SPAdes** and **IDBA-UD** have their own preferences for how sequence data is provided to them. To begin, we will look at the types of data accepted by **SPAdes**:

```bash
module load SPAdes/3.13.1-gimkl-2018b

spades.py -h
# ...
#Input data:
#-1      <filename>      file with forward paired-end reads
#-2      <filename>      file with reverse paired-end reads
#-s      <filename>      file with unpaired reads
#--mp<#>-1       <filename>      file with forward reads for mate-pair library number <#> (<#> = 1,2,..,9)
#--mp<#>-2       <filename>      file with reverse reads for mate-pair library number <#> (<#> = 1,2,..,9)
#--hqmp<#>-1     <filename>      file with forward reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)
#--hqmp<#>-2     <filename>      file with reverse reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)
#--nxmate<#>-1   <filename>      file with forward reads for Lucigen NxMate library number <#> (<#> = 1,2,..,9)
#--nxmate<#>-2   <filename>      file with reverse reads for Lucigen NxMate library number <#> (<#> = 1,2,..,9)
#--sanger        <filename>      file with Sanger reads
#--pacbio        <filename>      file with PacBio reads
#--nanopore      <filename>      file with Nanopore reads
#--tslr  <filename>      file with TSLR-contigs
#--trusted-contigs       <filename>      file with trusted contigs
#--untrusted-contigs     <filename>      file with untrusted contigs
```

At a glance, you could provide any of the following data types to **SPAdes** and have it perform an assembly:

1. Illumina paired-end sequencing data, either as standard library or Mate Pairs
1. Sanger sequences
1. PacBio reads
1. Oxford Nanopore reads
1. Pre-assembled scaffolds for guiding the assembly

Awkwardly, while **SPAdes** accepts multiple input libraries (i.e. samples) in a single assembly, this behaviour does not work with the **-meta** flag enabled. We therefore need to concatenate our four individual samples together ready for sequencing.

```bash
cat sample1_R1.fastq.gz sample2_R1.fastq.gz sample3_R1.fastq.gz sample4_R1.fastq.gz > for_spades_R1.fq.gz
cat sample1_R2.fastq.gz sample2_R2.fastq.gz sample3_R2.fastq.gz sample4_R2.fastq.gz > for_spades_R2.fq.gz
```

Note that these *fastQ* files are compressed, yet we can concatenate them together with the **cat** command regardless. This is a nice feature of *gz* files that is handy to remember.

By contrast, what does **IDBA-UD** accept?

```bash
module load IDBA/1.1.3-gimkl-2017a

idba_ud
# ...
# -o, --out arg (=out) output directory
# -r, --read arg       fasta read file (<=128)
# ...
# -l, --long_read arg  fasta long read file (>128)
```

'Short' or 'long' reads, and only a single file for each. This means that if we want to assembly our community data using **IDBA-UD** we will need to pool the paired-end data into a single, interleaved *fastA* file. Interleaved means that instead of having a pair of files that contain the separate forward and reverse sequences, the read pairs are in a single file in alternating order. For example

```bash
# Paired-end file, forward
>read1_1
...
>read2_1
...

# Paired-end file, forward
>read1_2
...
>read2_2
...

# Interleaved file
>read1_1
...
>read1_2
...
>read2_1
...
>read2_2
...
```

Fortunately, the **IDBA** set of tools comes with some helper scripts to achieve just this. Unfortunately we cannot apply this shuffling operation to compressed data, so we must decompress the data first.

```bash
module load pigz/2.4-GCCcore-7.4.0

for i in sample1 sample2 sample3 sample4;
do
  pigz --keep --decompress ${i}_R1.fastq.gz ${i}_R2.fastq.gz
  fq2fa --merge ${i}_R1.fastq ${i}_R2.fastq ${i}.fna
done

fq2fa --merge sample1_R1.fastq.gz sample1_R2.fastq.gz asample1.fa
cat sample1.fna sample2.fna sample3.fna sample4.fna > for_idba.fna
    
```

### Basic assembly parameters

For any assembler, there are a **_lot_** of parameters that can be fine-tuned depending on your data. As no two data sets are the same, it is almost impossiblet o predict which parameter combinations will yield the best outcome for your data set. That said, assembly can be quite a resource-intensive process and it is generally not practical to test every permutation of parameter value with your data. In genomics, the saying goes that the best assembly is the one that answers your question. As long as the data you are receiving is meaningful to the hypothesis you are seeking to address then your assembly is as good as it needs to be.

Generally speaking, assemblers are developed in a way where they run with default parameters that have been emperically showed to produce the best outcome **_on average_** across multiple data sets. For most purposes, there is not a lot of need to change these, but some parameters that we would always want to look at include:

1. *k*-mer sizes to be assembled over, and step size if using a range
1. Number of threads to use during assembly
1. Memory limit to prevent the assembler from using up all available RAM and forcing the computer to use its [swap space](https://web.mit.edu/rhel-doc/5/RHEL-5-manual/Deployment_Guide-en-US/ch-swapspace.html)

#### Setting the *k*-mer size

Depending on which assembler you are using the commands for chosing the *k*-mer sizes for assembly vary slightly, but they are recognisable between programs. In **SPAdes**, you can set the *k*-mer size using either

```bash
spades.py -k 21,33,55,77,99,121 ...

spades.py -k auto ...
```

Which amount to either specifying the *k*-mers ourselves, or letting **SPAdes** pick the sizes it thinks are best. For **IDBA-UD**, we select *k*-mer size using

```bash
idba_ud --mink 21 --maxk 121 --step 22
```

Unlike **SPAdes**, we do not have fine-scale control over the *k*-mer sizes used in assembly. We instead provide **IDBA-UD** with the first and last *k*-mer sizes to use, then specify the increment to use between these. In either case, it is important that we are always assembling using a *k*-mer of uneven (odd) length in order to avoid the creation of palindromic *k*-mers.

#### Specifying the number of threads

This is simple in either assembler:

```bash
spades.py -t 20 ...

idba_ud --num_threads 20 ...
```

The only thing to keep in mind is that these tools have different default behaviour. If no thread count is specified by the user, **SPAdes** will assemble with 16 threads. **IDBA-UD** will use all available threads which is can be problematic if you are using a shared compute environment that does not use a resource management system like slurm.

#### Setting a memory limit

By far, the worst feature of **SPAdes** is the high memory requirements for performing assembly. In the absence of monitoring, **SPAdes** will request more and more memory as it proceeds. If this requires more memory than is available on your computer your system will start to store memory to disk space. This is an extremely slow operation and can make your computer to effectively unusable. In managed environments such as NeSI a memory limit is imposed upon all running jobs, but if you are not using such a system you are advised to set a memory limit when executing **SPAdes**:

```bash
spades.py -m 400GB ...
```

No such parameter exists in **IDBA-UD**, although as this tool needs far less RAM than **SPAdes** you are less likely to need it.

## Preparing an assembly job for slurm

NeSI does not allow users to execute large jobs interactively on the terminal. Instead, the node that we have logged in to (*lander02*) has only a small fraction of the computing resources that NeSI houses. The *lander* node is used to write small command scripts, which are then deployed to the large compute nodes by a system called **slurm**. The ins and outs of working in **slurm** are well beyond the scope of this workshop (and may not be relevant if your institution uses a different resource allocation system). In this workshop, we will therefore only be showing you how to write minimal **slurm** scripts sufficient to achieve our goals. By the end of the workshop, you should have built up a small collection of **slurm** scripts for performing the cnessary stages of our workflow and with experience you will be able to modify these to suit your own needs.

To begin, we need to open a text file using the **nano** text editor. 

```bash
nano assembly.sl
```

