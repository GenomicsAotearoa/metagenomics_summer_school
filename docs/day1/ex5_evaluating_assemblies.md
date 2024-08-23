# Assembly evaluation

!!! info "Objectives"

    * [Evaluating the resource consumption of various assemblies](#evaluating-the-resource-consumption-of-various-assemblies)
    * [Evaluating the assemblies using `BBMap`](#evaluating-the-assemblies-using-bbmap)
    * [*(Optional)* Evaluating assemblies using `MetaQUAST`](#optional-evaluating-assemblies-using-metaquast)

---

<center>
![image](../theme_images/eval_assembly.png){width="450"}
</center>

## Evaluating the resource consumption of various assemblies

Check to see if your assembly jobs have completed. If you have multiple jobs running or queued, the easiest way to check this is to simply run the `squeue` command.

!!! terminal-2 "Check job progress"

    ```bash
    squeue --me
    ```

!!! circle-check "Terminal output"

    ```
    JOBID         USER     ACCOUNT   NAME        CPUS MIN_MEM PARTITI START_TIME     TIME_LEFT STATE    NODELIST(REASON)    
    39035482      jboe440  nesi02659 spawner-jupy   2      4G interac 2023-08-31T1     7:47:42 RUNNING  wbn004      
    ```

If there are no jobs besides your Jupyter session listed, either everything running has completed or failed. To get a list of all jobs we have run in the last day, we can use the `sacct` command. By default this will report all jobs for the day but we can add a parameter to tell the command to report all jobs run since the date we are specifying.

!!! terminal-2 "Check progress for jobs started on specific date"

    ```bash
    sacct -S 2023-08-12
    ```

!!! circle-check "Terminal output"

    ```
    JobID           JobName          Alloc     Elapsed     TotalCPU  ReqMem   MaxRSS State      
    --------------- ---------------- ----- ----------- ------------ ------- -------- ---------- 
    38483216        spawner-jupyter+     2    07:45:01     00:00:00      4G          NODE_FAIL  
    38483216.batch  batch                2    07:45:01     00:00:00                  CANCELLED  
    38483216.extern extern               2    07:45:01     00:00:00                  CANCELLED  
    38485254        spades_assembly     12    00:14:38     01:56:40     10G          COMPLETED  
    38485254.batch  batch               12    00:14:38     01:56:40         7227872K COMPLETED  
    38485254.extern extern              12    00:14:38     00:00:00                0 COMPLETED
    ```

Each job has been broken up into several lines, but the main ones to keep an eye on are the base `JobID` values. 

??? circle-info "Using `srun`"

    If you use `srun`, the JobID will have values suffixed with *.0*. The first of these references the complete job. The later (and any subsequent suffixes like *.1*, *.2*) are the individual steps in the script that were called with the `srun` command.

We can see here the time elapsed for each job, and the number of CPU hours used during the run. If we want a more detailed breakdown of the job we can use the `nn_seff` command

!!! terminal-2 "Check job resource use"

    ```bash
    nn_seff 38485254
    ```

!!! circle-check "Terminal output"

    ```
    Cluster: mahuika
    Job ID: 38485254
    State: COMPLETED
    Cores: 6
    Tasks: 1
    Nodes: 1
    Job Wall-time:   48.8%  00:14:38 of 00:30:00 time limit
    CPU Efficiency: 132.9%  01:56:40 of 01:27:48 core-walltime
    Mem Efficiency:  68.9%  6.89 GB of 10.00 GB0
    ```

Here we see some of the same information, but we also get some information regarding how well our job used the resources we allocated to it. You can see here that my CPU and memory usage was somewhat efficient but had high memory efficiency. In the future, I can request less time and retain the same RAM and still had the job run to completion.

CPU efficiency is harder to interpret as it can be impacted by the behaviour of the program. For example, mapping tools like `bowtie` and `BBMap` can more or less use all of their threads, all of the time and achieve nearly 100% efficiency. More complicated processes, like those performed in `SPAdes` go through periods of multi-thread processing and periods of single-thread processing, drawing the average efficiency down.

---

## Evaluating the assemblies using `BBMap`

Evaluating the quality of a raw metagenomic assembly is quite a tricky process. Since, by definition, our community is a mixture of different organisms, the genomes from some of these organisms assemble better than those of others. It is possible to have an assembly that looks 'bad' by traditional metrics that still yields high-quality genomes from individual species, and the converse is also true.

A few quick checks I recommend are to see how many contigs or scaffolds your data were assembled into, and then see how many contigs or scaffolds you have above a certain minimum length threshold. We will use `seqmagick` for performing the length filtering, and then just count sequence numbers using `grep`.

These steps will take place in the `4.evaluation/` folder, which contains copies of our `SPAdes` and `IDBA-UD` assemblies.

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal "code"

    ```bash
    # Load seqmagick
    module purge
    module load seqmagick/0.8.4-gimkl-2020a-Python-3.8.2

    # Navigate to working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/4.evaluation/

    # Filter assemblies and check number of contigs
    seqmagick convert --min-length 1000 spades_assembly/spades_assembly.fna \
                                        spades_assembly/spades_assembly.m1000.fna
    grep -c '>' spades_assembly/spades_assembly.fna spades_assembly/spades_assembly.m1000.fna

    seqmagick convert --min-length 1000 idbaud_assembly/idbaud_assembly.fna \
                                        idbaud_assembly/idbaud_assembly.m1000.fna
    grep -c '>' idbaud_assembly/idbaud_assembly.fna idbaud_assembly/idbaud_assembly.m1000.fna
    ```

!!! circle-check "Terminal output"

    === "`SPAdes`"

        ```
        spades_assembly/spades_assembly.fna:1318
        spades_assembly/spades_assembly.m1000.fna:933
        ```

    === "`IDBA-UD`"

        ```
        idbaud_assembly/idbaud_assembly.fna:5057
        idbaud_assembly/idbaud_assembly.m1000.fna:1996
        ```

If you have your own assemblies and you want to try inspect them in the same way, try that now. Note that the file names will be slightly different to the files provided above. If you followed the exact commands in the previous exercise, you can use the following commands.

!!! terminal "code"

    ```bash
    seqmagick convert --min-length 1000 ../3.assembly/my_spades_assembly/scaffolds.fasta my_spades_assembly.m1000.fna

    seqmagick convert --min-length 1000 ../3.assembly/my_idbaud_assembly/scaffold.fa my_idbaud_assembly.m1000.fna
    ```

!!! note "Choice of software: sequence file manipulation"

    The tool `seqtk` is also available on NeSI and performs many of the same functions as `seqmagick`. My choice of `seqmagick` is mostly cosmetic as the parameter names are more explicit so it's easier to understand what's happening in a command when I look back at my log files. Regardless of which tool you prefer, we strongly recommend getting familiar with either `seqtk` or `seqmagick` as both perform a lot of common FASTA and FASTQ file manipulations.

As we can see here, the `SPAdes` assembly has completed with fewer contigs assembled than the `IDBA-UD`, both in terms of total contigs assembled and contigs above the 1,000 bp size. This doesn't tell us a lot though - has `SPAdes` managed to assemble fewer reads, or has it managed to assemble the sequences into longer (and hence fewer) contigs? We can check this by looking at the N50/L50 (see more information about this statistic [here](https://www.molecularecologist.com/2017/03/29/whats-n50/)) of the assembly with `BBMap`.

!!! terminal "code"

    ```bash
    # Load BBMap module
    module purge
    module load BBMap/39.01-GCC-11.3.0

    # Generate statistics for filtered SPAdes assembly
    stats.sh in=spades_assembly/spades_assembly.m1000.fna
    ```

This gives quite a verbose output:

!!! circle-check "Terminal output"

    ```bash
    A       C       G       T       N       IUPAC   Other   GC      GC_stdev
    0.2536  0.2466  0.2462  0.2536  0.0019  0.0000  0.0000  0.4928  0.0960

    Main genome scaffold total:             933
    Main genome contig total:               2710
    Main genome scaffold sequence total:    34.300 MB
    Main genome contig sequence total:      34.236 MB       0.186% gap
    Main genome scaffold N/L50:             52/158.668 KB
    Main genome contig N/L50:               107/72.463 KB
    Main genome scaffold N/L90:             302/15.818 KB
    Main genome contig N/L90:               816/4.654 KB
    Max scaffold length:                    1.221 MB
    Max contig length:                      1.045 MB
    Number of scaffolds > 50 KB:            151
    % main genome in scaffolds > 50 KB:     76.76%


    Minimum         Number          Number          Total           Total           Scaffold
    Scaffold        of              of              Scaffold        Contig          Contig  
    Length          Scaffolds       Contigs         Length          Length          Coverage
    --------        --------------  --------------  --------------  --------------  --------
        All                    933           2,710      34,299,647      34,235,702    99.81%
       1 KB                    933           2,710      34,299,647      34,235,702    99.81%
     2.5 KB                    745           2,458      33,980,511      33,921,524    99.83%
       5 KB                    579           2,142      33,383,109      33,329,777    99.84%
      10 KB                    396           1,605      32,059,731      32,022,009    99.88%
      25 KB                    237             936      29,559,828      29,540,698    99.94%
      50 KB                    151             593      26,330,017      26,317,988    99.95%
     100 KB                     91             411      22,108,846      22,100,263    99.96%
     250 KB                     29             141      12,338,782      12,335,701    99.98%
     500 KB                      7              38       5,611,200       5,610,890    99.99%
       1 MB                      1               2       1,221,431       1,221,421   100.00%
    ```

!!! danger "N50 and L50 in `BBMap`"

    Unfortunately, the N50 and L50 values generated by `stats.sh` are switched. N50 should be a length and L50 should be a count. The results table below shows the corrected values based on `stats.sh` outputs.

But what we can highlight here is that the statistics for the `SPAdes` assembly, with short contigs removed, yielded an N50 of 72.5 kbp at the contig level. We will now compute those same statistics from the other assembly options.

!!! terminal "code"

    ```bash
    stats.sh in=spades_assembly/spades_assembly.fna

    stats.sh in=idbaud_assembly/idbaud_assembly.m1000.fna
    stats.sh in=idbaud_assembly/idbaud_assembly.fna
    ```

|Assembly|N50 (contig)|L50 (contig)|
|:---|:---:|:---:|
|**SPAdes** (filtered)|72.9 kbp|107 |
|**SPAdes** (unfiltered)|72.3 kbp|108 |
|**IDBA-UD** (filtered)|103.9 kbp|82 |
|**IDBA-UD** (unfiltered)|96.6 kbp|88 |

## Sequence taxonomic classification using `Kraken2`

Most, if not all, of the time, we never know the taxonomic composition of our metagenomic assemblies *a priori*. For some environments, we can make good guesses (e.g., *Prochlorococcus* in marine samples, members of the Actinobacteriota in soil samples, various Bacteroidota in the gut microbiome, etc.). Here, we can use `Kraken2`, a k-mer based taxonomic classifier to help us interrogate the taxonomic composition of our samples. This is helpful if there are targets we might be looking for (e.g. working hypotheses, well characterised microbiome) or as a check on what we might be missing out on after binning (taught in [day 2](../day2/ex6_initial_binning.md) and [day 3](../day3/ex11_coverage_and_taxonomy.md)).

!!! terminal "code"

    ```bash linenums="1"
    nano kraken2.sl
    ```

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      kraken2
    #SBATCH --partition     milan
    #SBATCH --time          10:00
    #SBATCH --mem           80G
    #SBATCH --cpus-per-task 20
    #SBATCH --error         %x.%j.err
    #SBATCH --output        %x.%j.out

    # Load module
    module purge
    module load Kraken2/2.1.3-GCC-11.3.0

    # Point to database
    K2DB=/nesi/project/nesi02659/MGSS_2024/resources/databases/k2_standard_20240605

    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/4.evaluation

    # Run Kraken2
    kraken2 --threads $SLURM_CPUS_PER_TASK \
            --classified-out spades_assembly.m1000.k2_classified.fna \
            --unclassified-out spades_assembly.m1000.k2_unclassified.fna \
            --report spades_assembly.m1000.k2_report.txt \
            --output spades_assembly.m1000.k2_out \
            --db ${K2DB} \
            spades_assembly/spades_assembly.m1000.fna
    ```

The main output `spades_assembly.m1000.k2_out` is quite dense and verbose with columns indicated [here](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format). The report `spades_assembly.m1000.k2_report.txt` is more human-readable, with nicely spaced columns that indicate: 

1. Percentage reads mapped to taxon
2. Number of reads mapped to taxon
3. Number of reads directly assigned to taxon
4. Rank of taxon
5. NCBI taxonomy ID
6. Scientific name of taxon

We also get the sequences that were classified and unclassified that can be used for further analyses (e.g. coverage estimation) if required.

As `Kraken2` classifications are *k*-mer based, we can also classify reads. This can be helpful if trying to filter out reads that may belong to taxonomic classifications that you're not interested in. When using reads for classification, we can also estimate the abundance of reads that belong to those taxa using Bracken.

!!! terminal "code"

    ```bash
    nano kraken2_bracken.sl
    ```

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      kraken2_bracken
    #SBATCH --partition     milan
    #SBATCH --time          30:00
    #SBATCH --mem           80G
    #SBATCH --cpus-per-task 20
    #SBATCH --error         %x.%j.err
    #SBATCH --output        %x.%j.out

    # Load modules
    module purge
    module load \
        Kraken2/2.1.3-GCC-11.3.0 \
        Bracken/2.7-GCC-11.3.0

    # Point to database
    K2DB=/nesi/project/nesi02659/MGSS_2024/resources/databases/k2_standard_20240605

    # Create output directory
    mkdir -p read_classification/

    # Run Kraken2 and Bracken
    for r1 in ../3.assembly/sample?_R1.fastq.gz; do
      # Output basename
      outbase=read_classification/$(basename ${r1} _R1.fastq.gz)
      # Taxonomic classification of reads
      kraken2 --paired --threads $SLURM_CPUS_PER_TASK \
              --classified-out ${outbase}#.k2_classified.fq \
              --unclassified-out ${outbase}#.k2_unclassified.fq \
              --report ${outbase}.k2_report.txt \
              --output ${outbase}.k2_out \
              --db ${K2DB} \
              ${r1} ${r1/R1/R2}
      # Estimate taxa abundance
      bracken -d ${K2DB} -i ${outbase}.k2_report.txt -r 100 \
              -o ${outbase}.bracken
    done
    ```

!!! note "Remember to produce the `--report` as Bracken bases its estimation on the report"

The modified code tells `Kraken2` that the inputs `${r1}` and `${r1/R1/R2}` are paired-end reads. The `#` after `${base}` will be replaced with the read orientation (i.e., forward `_1` or reverse `_2`).

Outputs for this run are similar to that of the assembly, with the major difference being paired-end read (fastq) outputs for classified and unclassified sequences.

The Bracken outputs provides us with adjusted number (columns `added_reads` and `new_est_reads`) and fraction of reads that were assigned to each species identified by Kraken2 for downstream analyses.

??? warning "Database construction for `Kraken2` and `Bracken`"

    For this workshop (and most applications), the standard or other extensive [pre-built Kraken2 databases](https://benlangmead.github.io/aws-indexes/k2) are sufficient. However, if you require specific reference sequences to be present, the database building process (outlined [here](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases)) can be quite time and resource intensive.

    The Kraken2 database also comes with a few pre-built Bracken databases (filenames look like this: `database100mers.kmer_distrib`). Note that the "100mers" doesn't refer to k-mers (oligonucleotide frequency), but it is the length of the read the database was built for. In this workshop, we've used the 100bp database to estimate abundances for 126bp reads. [According to the developers, this is acceptable](https://github.com/jenniferlu717/Bracken/issues/260). As long as the read length of the built database is $\le$ than the length of reads to be classified, Bracken will provide good enough estimates.

    If you want to build your own Bracken database based on your library minimum read lengths, you can follow the process [here](https://github.com/jenniferlu717/Bracken#step-0-build-a-kraken-1krakenuniqkraken2-database). Take note that you will require the initial sequence and taxonomy files used for building the `Kraken2` database (not part of the files in the pre-built databases) and these (especially for NCBI nt and bacterial databses) will take a long time to download and lots of memory to build.

## Reconstruct rRNA using `PhyloFlash`-`EMIRGE`

Another method that we can use to obtain taxonomic composition of our metagenomic libraries is to extract/reconstruct the 16S ribosomal RNA gene. However, modern assemblers often struggle to produce contiguous 16S rRNA genes due to (1) taxa harbouring multiple copies of the gene and (2) repeat regions that short reads cannot span ([see here for details](https://github.com/ablab/spades/issues/803)). As such, software such as [EMIRGE](https://github.com/csmiller/EMIRGE) ([Miller et al., 2011](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r44)) and [PhyloFlash](https://github.com/HRGV/phyloFlash) ([Gruber-Vodicka et al., 2020](https://journals.asm.org/doi/10.1128/msystems.00920-20)) were developed to reconstruct 16S rRNA genes. These methods leverage the expansive catalogue of 16S rRNA genes available in databases such as SILVA in order to subset reads and then reconstruct the full-length gene.

For this workshop, we will use PhyloFlash to obtain sequences and abundances of 16S rRNA sequences. The reason being that PhyloFlash can also run EMIRGE as part of its routine, compare those reconstructions with sequences generated by PhyloFlash's map-assemble, and produce a set of 16S rRNA gene sequences as well as estimate relative abundances via coverage estimation.

!!! note "Limitations"

    Both EMIRGE and PhyloFlash are closed reference methods, meaning they cannot produce sequences from novel taxa (i.e., taxa that do not exist in 16S rRNA databases). In this regard, the other reliable option to obtain sequences from novel taxa is by amplicon sequencing. 

!!! warning "`-lib`"

    A quirk of PhyloFlash is that all outputs will be placed in the current running directory. The only punctuation allowed is "-" and "_". Make sure to name your "LIB" sensibly so you can track which library you're working with.

!!! terminal "code"

    ```bash linenums="1"
    nano phyloflash.sl
    ```

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      phyloflash
    #SBATCH --partition     milan
    #SBATCH --time          15:00
    #SBATCH --mem           16G
    #SBATCH --cpus-per-task 12
    #SBATCH --error         %x.%j.err
    #SBATCH --output        %x.%j.out

    module purge
    module load Miniconda3/23.10.0-1
    module load USEARCH/11.0.667-i86linux32

    # Set conda variables
    source $(conda info --base)/etc/profile.d/conda.sh
    export PYTHONNOUSERSITE=1

    # Activate environment
    export CONDA_ENVS_PATH=/nesi/project/nesi02659/MGSS_2024/resources/tools
    source activate phyloflash

    # Run PhyloFlash per sample pair
    for r1 in ../3.assembly/sample*_R1.fastq.gz; do
        phyloFlash.pl -lib $(basename $r1 _R1.fastq.gz) \
                      -read1 ${r1} -read2 ${r1/R1/R2} \
                      -everything \
                      -CPUs $SLURM_CPUS_PER_TASK
    done

    phyloFlash_compare.pl \
        --zip $(ls -1 sample*.phyloFlash.tar.gz | paste -sd ",") \
        --task heatmap,barplot,matrix,ntu_table \
        --out allsamples_compare
    
    conda deactivate
    ```

A successful run of PhyloFlash (lines 24-29) would have generated the following per-sample outputs:

* `tar.gz` archive containing relevant data files
* An `.html` page summary
* A timestamped `.log` detailing what was done

The reconstructed 16S rRNA sequences are stored in the per-sample archive. We will extract them to explore the data files inside.

!!! terminal "code"

    ```sh
    for pf in *phyloFlash.tar.gz; do
        outdir=$(basename ${pf} .tar.gz)
        mkdir -p ${outdir}
        tar -xvzf ${pf} -C ${outdir}
    done
    ```

The reconstructed sequences based on the `phyloFlash` and `EMIRGE` workflows are stored in `sample?.phyloFlash/sample?.all.final.fasta`.

In lines 32-34 of our slurm script, we also ran a cross-sample comparison script. This generated 4 files:

* `allsamples_compare.barplot.pdf` A relative abundance barplot of identified taxa
* `allsamples_compare.heatmap.pdf` A heatmap of identified taxa
* `allsamples_compare.matrix.tsv` A dissimilarity matrix (abundance-weighted, UniFrac-like)
* `allsamples_compare.ntu_table.tsv` (long/tidy-format table of relative abundance per sample per taxa)

The `allsamples_compare.ntu_table.tsv` is especially useful if you intend to perform specific downstream analyses in `R`.

## *(Optional)* Evaluating assemblies using `MetaQUAST`

For more genome-informed evaluation of the assembly, we can use the `MetaQUAST` tool to view our assembled metagenome. This is something of an optional step because, like `QUAST`, `MetaQUAST` aligns your assembly against a set of reference genomes. Under normal circumstances we wouldn't know the composition of the metagenome that led to our assembly. In this instance determining the optimal reference genomes for a `MetaQUAST` evaluation is a bit of a problem. For your own work, the following tools could be used to generate taxonomic summaries of your metagenomes to inform your reference selection:

!!! abstract ""
    1. [Kraken2](https://ccb.jhu.edu/software/kraken2/) (DNA based, *k*-mer classification)
    1. [CLARK](http://clark.cs.ucr.edu/) (DNA based. *k*-mer classification)
    1. [Kaiju](http://kaiju.binf.ku.dk/) (Protein based, BLAST classification)
    1. [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml) (DNA based, sequence alignment classification)
    1. [MeTaxa2](https://microbiology.se/software/metaxa2/) or [SingleM](https://github.com/wwood/singlem) (DNA based, 16S rRNA recovery and classification)
    1. [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) (DNA based, clade-specific marker gene classification)

A good summary and comparison of these tools (and more) was published by [Ye *et al.*](https://www.ncbi.nlm.nih.gov/pubmed/31398336).

However, since we **_do_** know the composition of the original communities used to build this mock metagenome, `MetaQUAST` will work very well for us today. In your `4.evaluation/` directory you will find a file called `ref_genomes.txt`. This file contains the names of the genomes used to build these mock metagenomes. We will provide these as the reference input for `MetaQUAST`.

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      metaquast
    #SBATCH --partition     milan
    #SBATCH --time          00:15:00
    #SBATCH --mem           10GB
    #SBATCH --cpus-per-task 10
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Load module
    module purge
    module load QUAST/5.0.2-gimkl-2018b
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/4.evaluation
    
    # Run metaquast    
    metaquast.py --references-list ref_genomes.txt --max-ref-number 21 \
                 -t $SLURM_CPUS_PER_TASK \
                 --labels SPAdes,SPAdes.m1000,IDBAUD,IDBAUD.m1000 \
                 --output-dir metaquast_results/ \
                 spades_assembly/spades_assembly.fna \
                 spades_assembly/spades_assembly.m1000.fna \
                 idbaud_assembly/idbaud_assembly.fna \
                 idbaud_assembly/idbaud_assembly.m1000.fna
    ```

By now, you should be getting familiar enough with the console to understand what most of the parameters here refer to. The one parameter that needs explanation is the `--max-ref-number` flag, which we have set to 21. This caps the maximum number of reference genomes to be downloaded from NCBI which we do in the interest of speed. Since there are 21 taxa in the file `ref_genomes.txt` (10 prokaryote species and 11 viruses), `MetaQUAST` will download one reference genome for each. If we increase  the `--max-ref-number` flag we will start to get multiple reference genomes per taxa provided which is usually desirable.

We will now look at a few interesting assembly comparisons.

!!! note "Viewing HTML in Jupyter Hub"

    The NeSI Jupyter hub does not currently support viewing HTML that require Javascript (even if the browser you are running it in does). To view a basic version of the report, download the report file by navigating to the `4.evaluation/quast_results/` folder, right-click `report.html/` and select download. The downloaded file will then open within a new tab in the browser. 
    
    !!! warning ""
    
        Rendering the full report requires the other folders from within `quast_results/` to also be downloaded and available in the same directory as `report.html`. Unfortunately, the Jupyter hub environment does not appear to currently support downloading entire folders using this method.

An example of the `MetaQUAST` output files are also available for download. You will need to download both [references](../resources/quast_references.zip) and [results](../resources/quast_results_sans_reference.zip). Unzip both within the same directory.

!!! success ""

    === "Brief summary<br>of assemblies"

        ![image](../figures/ex5_fig1_shortsummary_2022.png)

    === "Comparison of NGA50<br>between assemblies"

        ![image](../figures/ex5_fig2_nga50_2022.png)

    === "Comparison of<br>aligned contigs"

        ![image](../figures/ex5_fig3_contigsmatched_2022.png)

    === "Inspection of<br>unaligned contigs"

        ![image](../figures/ex5_fig4_contigsunmatched_2022.png)

---
