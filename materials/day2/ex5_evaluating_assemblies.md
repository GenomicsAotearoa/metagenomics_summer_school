# Evaluating the assemblies

### Objectives

* Evaluate the resource consumption of various assemblies
* Evaluate the assemblies
* Future considerations

---

### Evaluate the resource consumption of various assemblies

Check to see if your jobs from last night have completed. If you have multiple jobs running or queued, the easiest way to check this is to simply run the `squeue` command from yesterday.

```bash
squeue -u <user name>

#  JOBID     USER ACCOUNT            NAME  ST REASON    START_TIME                TIME TIME_LEFT NODES CPUS
```

Since there are no jobs listed, either everything running has completed or failed. To get a list of all jobs we have run in the last day, we can use the `sacct` command. By default this will report all jobs for the day but we can add a parameter to tell the command to report all jobs run since the date we are specifying.

```bash
sacct -S 2019-12-10

#         JobID         JobName     Elapsed     TotalCPU Alloc   MaxRSS      State
#-------------- --------------- ----------- ------------ ----- -------- ----------
#8744675        spades_assembly    00:15:16     01:41:43    10          COMPLETED
#8744675.batch  batch              00:15:16    00:00.547    10    3908K COMPLETED
#8744675.extern extern             00:15:16     00:00:00    10        0 COMPLETED
#8744675.0      spades.py          00:15:15     01:41:42    10 6260072K COMPLETED
#8744677        idbaud_assembly    00:11:06     01:27:42    10          COMPLETED
#8744677.batch  batch              00:11:06    00:00.477    10    3760K COMPLETED
#8744677.extern extern             00:11:06     00:00:00    10        0 COMPLETED
#8744677.0      idba_ud            00:11:05     01:27:42    10 2541868K COMPLETED
```

Each job has been broken up into several lines, but the main ones to keep an eye are the base JobID values, and the values suffixed with *.0*. The first of these references the complete job. The later (and any subsequent suffixes like *.1*, *.2*) are the individual steps in the script that were called with the `srun` command.

We can see here the time elapsed for each job, and the number of CPU hours used during the run. If we want a more detailed breakdown of the job we can use the `seff` command

```bash
seff 8744675

#Job ID: 8744675
#Cluster: mahuika
#User/Group: dwai012/dwai012
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 10
#CPU Utilized: 01:41:44
#CPU Efficiency: 66.64% of 02:32:40 core-walltime
#Job Wall-clock time: 00:15:16
#Memory Utilized: 5.97 GB
#Memory Efficiency: 29.85% of 20.00 GB
```

Here we see some of the same information, but we also get some information regarding how well our job used the resources we allocated to it. You can see here that my CPU and memory usage was not particularly efficient, in hindsight I could have request a lot less RAM and still had the job run to completion.

CPU efficiency is harder to interpret as it can be impacted by the behaviour of the program. For example, mapping tools like `bowtie` and `BBMap` can more or less use all of their threads, all of the time and achieve nearly 100% efficiency. More complicated processes, like those performed in `SPAdes` go through periods of multi-thread processing and periods of single-thread processing, drawing the average efficiency down.

---

### Evaluate the assemblies

Evaluating the quality of a raw metagenomic assembly is quite a tricky process. Since, by definition, our community is a mixture of different organisms, the genomes from some of these organisms assemble better than those of others. It is possible to have an assembly that looks 'bad' by traditional metrics that still yields high-quality genomes from individual species, and the converse is also true.

A few quick checks I recommend are to see how many contigs or scaffolds your data were assembled into, and then see how many contigs or scaffolds you have above a certain minimum length threshold. We will use `seqmagick` for performing the length filtering, and then just count sequence numbers using `grep`.

```bash
module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/4.evaluation/

seqmagick convert --min-length 1000 spades_assembly/spades_assembly.fna \
                                    spades_assembly/spades_assembly.m1000.fna
grep -c '>' spades_assembly/spades_assembly.fna spades_assembly/spades_assembly.m1000.fna
# spades_assembly/spades_assembly.fna:1913
# spades_assembly/spades_assembly.m1000.fna:1047

seqmagick convert --min-length 1000 idbaud_assembly/idbaud_assembly.fna \
                                    idbaud_assembly/idbaud_assembly.m1000.fna
grep -c '>' idbaud_assembly/idbaud_assembly.fna idbaud_assembly/idbaud_assembly.m1000.fna
# idbaud_assembly/idbaud_assembly.fna:4891
# idbaud_assembly/idbaud_assembly.m1000.fna:1901
```

If you have your own assemblies and you want to try inspect them in the same way, try that now. Note that the file names will be slightly different to the files provided above. If you followed the exact commands in the previous exercise, you can use the following commands.

```bash
seqmagick convert --min-length 1000 ../3.assembly/spades_assembly/scaffolds.fasta my_spades_assembly.m1000.fna

seqmagick convert --min-length 1000 ../3.assembly/idbaud_assembly/scaffold.fa my_idbaud_assembly.m1000.fna
```

*Note: The tool `seqtk` is also available on NeSI and performs many of the same functions as `seqmagick`. My choice of `seqmagick` is mostly cosmetic as the parameter names are more explicit so it's easier to understand what's happening in a command when I look back at my log files. Regardless of which tool you prefer, we strongly recommend getting familiar with either `seqtk` or `seqmagick` as both perform a lot of common *fastA* and *fastQ* file manipulations.*

As we can see here, the `SPAdes` assembly has completed with fewer contigs assembled than the `IDBA-UD`, both in terms of total contigs assembled and contigs above the 1,000 bp size. This doesn't tell us a lot though - has `SPAdes` managed to assemble fewer reads, or has it managed to assemble the sequences  into longer (and hence fewer) contigs? We can check this by looking at the N50/L50 of the assembly with `BBMap`.

```bash
module load BBMap/38.73-gimkl-2018b

stats.sh in=spades_assembly/spades_assembly.m1000.fna
```

This gives quite a verbose output:

```bash
A       C       G       T       N       IUPAC   Other   GC      GC_stdev
0.2525  0.2468  0.2461  0.2546  0.0019  0.0000  0.0000  0.4929  0.0973

Main genome scaffold total:             1047
Main genome contig total:               2797
Main genome scaffold sequence total:    33.635 MB
Main genome contig sequence total:      33.572 MB       0.189% gap
Main genome scaffold N/L50:             53/154.193 KB
Main genome contig N/L50:               104/75.816 KB
Main genome scaffold N/L90:             336/13.309 KB
Main genome contig N/L90:               839/4.343 KB
Max scaffold length:                    871.42 KB
Max contig length:                      779.486 KB
Number of scaffolds > 50 KB:            150
% main genome in scaffolds > 50 KB:     75.48%


Minimum         Number          Number          Total           Total           Scaffold
Scaffold        of              of              Scaffold        Contig          Contig
Length          Scaffolds       Contigs         Length          Length          Coverage
--------        --------------  --------------  --------------  --------------  --------
    All                  1,047           2,797      33,635,467      33,572,059    99.81%
   1 KB                  1,047           2,797      33,635,467      33,572,059    99.81%
 2.5 KB                    801           2,471      33,229,209      33,172,921    99.83%
   5 KB                    601           2,069      32,495,504      32,450,300    99.86%
  10 KB                    401           1,489      31,045,953      31,015,553    99.90%
  25 KB                    230             795      28,346,907      28,334,024    99.95%
  50 KB                    150             521      25,386,942      25,377,983    99.96%
 100 KB                     87             346      20,945,889      20,939,278    99.97%
 250 KB                     30             133      12,482,836      12,479,278    99.97%
 500 KB                      9              31       5,848,625       5,847,955    99.99%
```

But what we can highlight here is that the statistics for the `SPAdes` assembly, with short contigs removed, yielded an N50 of 104 kbp at the contig level. We will now compute those same statistics from the other assembly options

```bash
stats.sh in=spades_assembly/spades_assembly.fna

stats.sh in=idbaud_assembly/idbaud_assembly.m1000.fna
stats.sh in=idbaud_assembly/idbaud_assembly.fna
```

|Assembly|N50 (contig)|L50 (contig)|
|:---|:---:|:---:|
|**SPAdes** (filtered)|104 kbp|76 kbp|
|**SPAdes** (unfiltered)|106 kbp|76 kbp|
|**IDBA-UD** (filtered)|76 kbp|107 kbp|
|**IDBA-UD** (unfiltered)|83 kbp|101 kbp|

#### *Optional:* Evaluating assemblies using *MetaQUAST*

For more genome-informed evaluation of the assembly, we can use the `MetaQUAST` tool to view our assembled metagenome. This is something of an optional step because, like `QUAST`, `MetaQUAST` aligns your assembly against a set of reference genomes. Under normal circumstances we wouldn't know the composition of the metagenome that led to our assembly. In this instance determining the optimal reference genomes for a `MetaQUAST` evaluation is a bit of a problem. For your own work, the following tools could be used to generate taxonomic summaries of your metagenomes to inform your reference selection:

1. [Kraken2](https://ccb.jhu.edu/software/kraken2/) (DNA based, *k*-mer classification)
1. [CLARK](http://clark.cs.ucr.edu/) (DNA based. *k*-mer classification)
1. [Kaiju](http://kaiju.binf.ku.dk/) (Protein based, BLAST classification)
1. [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml) (DNA based, sequence alignment classification)
1. [MeTaxa2](https://microbiology.se/software/metaxa2/) or [SingleM](https://github.com/wwood/singlem) (DNA based, 16S rRNA recovery and classification)
1. [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) (DNA based, clade-specific marker gene classification)

A good summary and comparison of these tools (and more) was recently published by [Ye *et al.*](https://www.ncbi.nlm.nih.gov/pubmed/31398336).

However, since we **_do_** know the composition of the original communities used to build this mock metagenome, `MetaQUAST` will work very well for us today. In your `4.evaluation/` directory you will find a file called `ref_genomes.txt`. This file contains the names of the 9 genomes used to build these mock metagenomes. We will provide these as the reference input for `MetaQUAST`.

```bash
module load QUAST/5.0.2-gimkl-2018b

metaquast.py spades_assembly/spades_assembly.fna spades_assembly/spades_assembly.m1000.fna \
             idbaud_assembly/idbaud_assembly.fna idbaud_assembly/idbaud_assembly.m1000.fna \
             --references-list ref_genomes.txt --max-ref-number 10 -t 10
```

By now, you should be getting familiar enough with the console to understand what most of the parameters here refer to. The one parameter that needs explanation is the `--max-ref-number` flag, which we have set to 10. This caps the maximum number of reference genomes to be downloaded from NCBI which we do in the interest of speed. Since there are 10 species names in the file `ref_genomes.txt`, `MetaQUAST` will download one of each. If we increase the number we will start to get multiple references per name provided which is usually desirable.

We will now look at a few interesting assembly comparisons.

#### Brief summary of assemblies

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex5_fig1_shortsummary.PNG)

#### Comparison of NGA50 between assemblies

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex5_fig2_nga50.PNG)

#### Comparison of aligned contigs

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex5_fig3_contigsmatched.PNG)

#### Inspection of unaligned contigs

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex5_fig4_contigsunmatched.PNG)

---