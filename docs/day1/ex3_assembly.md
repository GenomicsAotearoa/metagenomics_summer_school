# Assembly

!!! info "Objectives"

    * [Become familiar with the standard input files for `SPAdes` and `IDBA-UD`](#the-standard-input-files-for-spades-and-idba-ud)
    * [Understand the basic parameters that should be modified when using these assemblers](#basic-assembly-parameters)
    * [Prepare an assembly job to run under slurm](#preparing-an-assembly-job-for-slurm)

All work for this exercise will occur in the `3.assembly/` directory.

---

### The standard input files for *SPAdes* and *IDBA-UD*

Although they both make use of the same types of data, both `SPAdes` and `IDBA-UD` have their own preferences for how sequence data is provided to them. To begin, we will look at the types of data accepted by `SPAdes`:

```bash
module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5
```
```bash
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

At a glance, you could provide any of the following data types to `SPAdes` and have it perform an assembly:

1. Illumina paired-end sequencing data, either as standard library or Mate Pairs
1. Sanger sequences
1. PacBio reads
1. Oxford Nanopore reads
1. Pre-assembled scaffolds for guiding the assembly

Awkwardly, while `SPAdes` accepts multiple input libraries (i.e. samples) in a single assembly, this behaviour does not work with the `-meta` flag enabled, which is needed in our example to activate metagenome assembly mode. We therefore need to concatenate our four individual samples together ready for sequencing.

```bash
cat sample1_R1.fastq.gz sample2_R1.fastq.gz sample3_R1.fastq.gz sample4_R1.fastq.gz > for_spades_R1.fq.gz
cat sample1_R2.fastq.gz sample2_R2.fastq.gz sample3_R2.fastq.gz sample4_R2.fastq.gz > for_spades_R2.fq.gz
```

Note that these *fastQ* files are compressed, yet we can concatenate them together with the `cat` command regardless. This is a nice feature of *.gz* files that is handy to remember.

By contrast, what does `IDBA-UD` accept?

```bash
module purge
module load IDBA-UD/1.1.3-GCC-11.3.0
```
```
idba_ud
# ...
# -o, --out arg (=out) output directory
# -r, --read arg       fasta read file (<=128)
# ...
# -l, --long_read arg  fasta long read file (>128)
```

'Short' or 'long' reads, and only a single file for each. This means that if we want to assemble our community data using `IDBA-UD` we will need to pool the paired-end data into a single, interleaved *fastA* file. Interleaved means that instead of having a pair of files that contain the separate forward and reverse sequences, the read pairs are in a single file in alternating order. For example

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

Fortunately, the `IDBA` set of tools comes with some helper scripts to achieve just this. Unfortunately we cannot apply this shuffling operation to compressed data, so we must decompress the data first.

```bash
module load pigz/2.7

for i in sample1 sample2 sample3 sample4;
do
  pigz --keep --decompress ${i}_R1.fastq.gz ${i}_R2.fastq.gz
  fq2fa --merge ${i}_R1.fastq ${i}_R2.fastq ${i}.fna
done

cat sample1.fna sample2.fna sample3.fna sample4.fna > for_idba.fna
    
```

---

### Basic assembly parameters

For any assembler, there are a **_lot_** of parameters that can be fine-tuned depending on your data. As no two data sets are the same, it is almost impossible to predict which parameter combinations will yield the best outcome for your dataset. That said, an assembly can be quite a resource-intensive process and it is generally not practical to test every permutation of parameter values with your data. In genomics, the saying goes that the best assembly is the one that answers your question. As long as the data you are receiving is meaningful to the hypothesis you are seeking to address, then your assembly is as good as it needs to be.

Generally speaking, assemblers are developed in a way where they run with default parameters that have been empirically demonstrated to produce the best outcome **_on average_** across multiple data sets. For most purposes, there is not a lot of need to change these, but some parameters that we would always want to look at include:

1. *k*-mer sizes to be assembled over, and step size if using a range
1. Number of threads to use during assembly
1. Memory limit to prevent the assembler from using up all available RAM and forcing the computer to use its [swap space](https://web.mit.edu/rhel-doc/5/RHEL-5-manual/Deployment_Guide-en-US/ch-swapspace.html)

#### Setting the *k*-mer size

Depending on which assembler you are using, the commands for choosing the *k*-mer sizes for the assembly vary slightly, but they are recognisable between programs. In `SPAdes`, you can set the *k*-mer size using either

```bash
spades.py -k 21,33,55,77,99,121 ...

spades.py -k auto ...
```

The first command lets us specify the *k*-mers ourselves, or we are letting `SPAdes` automatically pick the most appropriate size. For `IDBA-UD`, we can select the *k*-mer size using

```bash
idba_ud --mink 21 --maxk 121 --step 22
```

Unlike `SPAdes`, we do not have fine-scale control over the *k*-mer sizes used in the assembly. We instead provide `IDBA-UD` with the first and last *k*-mer size to use, then specify the increment to use between these. In either case, it is important that we are always assembling using a *k*-mer of uneven (odd) length in order to avoid the creation of palindromic *k*-mers.

#### Specifying the number of threads

This is simple in either assembler:

```bash
spades.py -t 20 ...

idba_ud --num_threads 20 ...
```

The only thing to keep in mind is that these tools have different default behaviour. If no thread count is specified by the user, `SPAdes` will assemble with 16 threads. `IDBA-UD` will use all available threads, which can be problematic if you are using a shared compute environment that does not use a resource management system like slurm.

#### Setting a memory limit

By far, the worst feature of `SPAdes` is the high memory requirement for performing an assembly. In the absence of monitoring, `SPAdes` will request more and more memory as it proceeds. If this requires more memory than is available on your computer, your system will start to store memory to disk space. This is an extremely slow operation and can render your computer effectively unusable. In managed environments such as NeSI a memory limit is imposed upon all running jobs, but if you are not using such a system you are advised to set a memory limit when executing `SPAdes`:

```bash
spades.py -m 400GB ...
```

No such parameter exists in `IDBA-UD`, but it requires far less RAM than `SPAdes`, so you are less likely to need it.

---

### Preparing an assembly job for slurm

NeSI does not allow users to execute large jobs interactively on the terminal. Instead, the node that we have logged in to (*lander02*) has only a small fraction of the computing resources that NeSI houses. The *lander* node is used to write small command scripts, which are then deployed to the large compute nodes by a system called **slurm**. The ins and outs of working in slurm are well beyond the scope of this workshop (and may not be relevant if your institution uses a different resource allocation system). In this workshop, we will therefore only be showing you how to write minimal slurm scripts sufficient to achieve our goals. By the end of the workshop, you should have built up a small collection of slurm scripts for performing the necessary stages of our workflow and with experience you will be able to modify these to suit your own needs.

#### Submitting a *SPAdes* job to NeSI using slurm

To begin, we need to open a text file using the `nano` text editor. 

```bash
cd 3.assembly/

nano assembly_spades.sl
```

Into this file, either write or copy/paste the following commands:

```bash
#!/bin/bash -e
#SBATCH --account nesi02659
#SBATCH --job-name spades_assembly
#SBATCH --res SummerSchool
#SBATCH --time 00:50:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --error spades_assembly.err
#SBATCH --output spades_assembly.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load SPAdes/3.15.3-gimkl-2020a

spades.py --meta -k 33,55,77,99,121 -t 12 -1 for_spades_R1.fq.gz -2 for_spades_R2.fq.gz -o spades_assembly/
```

To save your file, use `Ctrl + O` to save the file, then `Ctrl + X` to exit `nano`. Going through those lines one by one;

|Slurm parameter|Function|
|:---|:---|
|**#!/bin/bash -e**|Header for the file, letting NeSI know how to interpret the following commands. The `-e` flag means that the slurm run will halt at the first failed command (rather than pushing through and trying to execute subsequent ones)|
|**--account ...**|The name of the project account to run the run under. You are provided with this when you create a project on NeSI|
|**--job-name ...**|The name of the job, to display when using the `squeue` command|
|**--res ...*|This is a parameter that we need due to a system reservation for this workshop. For your own work, ignore this flag|
|**--time ...**|Maximum run time for the job before it is killed|
|**--mem ...**|Amount of server memory to allocate to the job. If this is exceeded, the job will be terminated|
|**--ntasks ...**|Number of tasks to distribute the job across. For all tools we use today, this will remain as 1|
|**--cpus-per-task ...**|number of processing cores to assign to the job. This should match with the number used by your assembler|
|**--error ...**|File to log the [standard error](https://en.wikipedia.org/wiki/Standard_streams) stream of the program.<br>This is typically used to prove error reports, or to just inform the user of job progress|
|**--output ...**|File to log the [standard output](https://en.wikipedia.org/wiki/Standard_streams) stream of the program.<br>This is typically used to inform the user of job progress and supply messages|

The `module load` command needs to be invoked within your slurm script. It is also a good idea to explicitly set the path to your files within the job so that

1. There is no chance of having the job fail immediately because it cannot find the relevant files
1. When looking back through your slurm logs, you know where the data is meant to be

When executing the `SPAdes` command, there are a few parameters to note here:

|Parameter|Function|
|:---|:---|
|**--meta**|Activate metagenome assembly mode. Default is to assemble your metagenome using single genome assembly assumptions|
|**-k ...**|*k*-mer sizes for assembly. These choices will provide the output we will use in the **Binning** session, but feel free to experiment with these to see if you can improve the assembly|
|**-1 ..**|Forward reads, matched to their reverse partners|
|**-2 ...**|Reverse reads, matched to their forward partners|
|**-o ...**|Output directory for all files|

Note that we also prefix the command (`spades.py`) with the `srun` command. This is a command specific to slurm and allows NeSI to track the resource usage of the `SPAdes` job.

We don't explicitly set memory or thread counts for this job, simply for the sake of keeping the command uncluttered. The default memory limit of `SPAdes` (250 GB) is much higher than the 10 GB we have allowed our job here. If the memory cap was violated then both slurm and `SPAdes` will terminate the assembly. We have also left the number of threads at the default value of 16, which matches the number specified in the slurm header.

It is a good idea to match your number of threads request in the slurm script with what you intend to use with `SPAdes` because your project usage is calculated based off what you request in your slurm scripts rather than what you actually use. Requesting many unused threads simply drives your project down the priority queue. By contrast, requesting fewer threads than you attempt to use in the program (i.e. request 10 in slurm, set thread count to 30 in `SPAdes`) will result in reduced performance, as your `SPAdes` job will divide up jobs as though it has 30 threads, but only 10 will be provided. This is discussed [in this blog post](https://www.codeguru.com/cpp/sample_chapter/article.php/c13533/Why-Too-Many-Threads-Hurts-Performance-and-What-to-do-About-It.htm).

Once you are happy with your slurm script, execute the job by navigating to the location of your script and entering the command

```bash
sbatch assembly_spades.sl
```

You will receive a message telling you the job identifier for your assembly. Record this number, as we will use it in the next exercise.

#### Monitoring job progress

You can view the status of your current jobs using the command

```bash
squeue -u <login name>

  JOBID     USER ACCOUNT            NAME  ST REASON    START_TIME                TIME TIME_LEFT NODES CPUS
8744675  dwai012 ga02676 spades_assembly  PD Priority  2019-11-29T16:37:43       0:00 00:30:00      1   16
```

We can see here that the job has not yet begun, as NeSI is waiting for resources to come available. At the stage the `START_TIME` is an estimation of when the resources are expected to become available. When they do, the output will change to

```bash
  JOBID     USER ACCOUNT            NAME  ST REASON    START_TIME                TIME TIME_LEFT NODES CPUS
8744675  dwai012 ga02676 spades_assembly   R None      2019-11-29T16:40:00       1:20 00:28:40      1   16
```

Which allows us to track how far into our run we are, and see the remaining time for the job. The `START_TIME` column now reports the time the job actually began.

#### Submitting an *IDBA-UD* job to NeSI using slurm

To run an equivalent assembly with `IDBA-UD`, create a new slurm script as follows

```bash
nano idbaud_assembly.sl
```

Paste or type in the following:

```bash
#!/bin/bash -e
#SBATCH --account nesi02659
#SBATCH --job-name idbaud_assembly
#SBATCH --res SummerSchool
#SBATCH --time 00:35:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --error idbaud_assembly.err
#SBATCH --output idbaud_assembly.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load IDBA/1.1.3-gimkl-2017a

idba_ud --num_threads 8 --mink 33 --maxk 99 --step 22 -r for_idba.fna -o idbaud_assembly/
```

And submit the script as a slurm job:

```bash
sbatch idbaud_assembly.sl
```

Remember to record your job identification number.

---
