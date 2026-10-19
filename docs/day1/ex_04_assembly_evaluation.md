# Assembly evaluation

!!! info "Objectives"

    * [Evaluating the resource consumption of various assemblies](#evaluating-the-resource-consumption-of-various-assemblies)
    * [Evaluating the assemblies using `BBMap`](#evaluating-the-assemblies-using-bbmap)

---

![image](../theme_images/eval_assembly.png){.center width="450"}

## Evaluating the resource consumption of various assemblies

Check to see if your assembly jobs have completed. If you have multiple jobs running or queued, the easiest way to check this is to simply run the `squeue` command.

!!! terminal-2 "Check job progress"

    ```bash
    squeue --me
    ```

!!! circle-check "Terminal output"

    ```
         JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
      10014267     milan spades_a    drost  R      56:00      1 mc094
    ```

If there are no jobs listed, either everything running has completed or failed. To get a list of all jobs we have run in the last day, we can use the `sacct` command. By default this will report all jobs for the day but we can add a parameter to tell the command to report all jobs run since the date we are specifying.

!!! terminal-2 "Check progress for jobs started on specific date"

    ```bash
    sacct -S 2026-10-16
    ```

!!! circle-check "Terminal output"

    ```
    JobID           JobName  Partition    Account  AllocCPUS      State ExitCode 
    ------------ ---------- ---------- ---------- ---------- ---------- -------- 
    9996006      spades_as+      milan  nesi02659         24     FAILED    255:0 
    9996006.bat+      batch             nesi02659         24     FAILED    255:0 
    9996006.ext+     extern             nesi02659         24  COMPLETED      0:0 
    9996013      spades_as+      milan  nesi02659         24    TIMEOUT      0:0 
    9996013.bat+      batch             nesi02659         24  CANCELLED     0:15 
    9996013.ext+     extern             nesi02659         24  COMPLETED      0:0 
    9996014      spades_as+      milan  nesi02659         24 CANCELLED+      0:0 
    9996014.bat+      batch             nesi02659         24  CANCELLED     0:15 
    9996014.ext+     extern             nesi02659         24  COMPLETED      0:0 
    10014267     spades_as+      milan  nesi02659         32    RUNNING      0:0 
    10014267.ba+      batch             nesi02659         32    RUNNING      0:0 
    10014267.ex+     extern             nesi02659         32    RUNNING      0:0 
    ```

Each job has been broken up into several lines, but the main ones to keep an eye on are the base `JobID` values. 

We can see here the time elapsed for each job, and the number of CPU hours used during the run. If we want a more detailed breakdown of the job we can use the `nn_seff` command

!!! terminal-2 "Check job resource use"

    ```bash
    nn_seff 9996013
    ```

!!! circle-check "Terminal output"

    ```
    Job ID: 9996013
    State: TIMEOUT
    Cores: 12
    Job Wall-time:        101%  00:30:12 of 00:30:00 time limit
    Avg CPU Utilisation:   16%  00:59:33 of 06:02:24 core-walltime
    Peak Mem Utilisation:  38%  3.75 GB of 10.00 GB
    ```

Here we see some of the same information, but we also get some information regarding how well our job used the resources we allocated to it. You can see here that my CPU and memory usage was inefficient yet the job ran out of time. In the future, I can request more time and fewer CPUs (cores). As the job did not finish I still don't know the typical memory requirement of this task.

CPU efficiency is harder to interpret as it can be impacted by the behaviour of the program. For example, mapping tools like `bowtie` and `BBMap` can more or less use all of their threads, all of the time and achieve nearly 100% efficiency. More complicated processes, like those performed in `SPAdes` go through periods of multi-thread processing and periods of single-thread processing, drawing the average efficiency down.

---

## Evaluating the assemblies using `BBMap`

Evaluating the quality of a raw metagenomic assembly is quite a tricky process. Since, by definition, our community is a mixture of different organisms, the genomes from some of these organisms assemble better than those of others. It is possible to have an assembly that looks 'bad' by traditional metrics that still yields high-quality genomes from individual species, and the converse is also true.

A few quick checks I recommend are to see how many contigs or scaffolds your data were assembled into, and then see how many contigs or scaffolds you have above a certain minimum length threshold. We will use `seqmagick` for performing the length filtering, and then just count sequence numbers using `grep`.

These steps will take place in the `4.evaluation/` folder, which contains copies of our `SPAdes` assemblies, as well as assemblies created with `IDBA-UD`, which is covered in the [bonus material](./ex_05_assembly_idba-ud.md).

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal "code"

    ```bash
    # Load seqmagick
    module purge
    module load seqmagick/0.8.4-gimkl-2020a-Python-3.8.2

    # Navigate to working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/4.evaluation/

    # Filter assemblies and check number of contigs
    seqmagick convert --min-length 1000 \
                      spades_assembly/spades_assembly.fna \
                      spades_assembly/spades_assembly.m1000.fna
    grep -c '>' spades_assembly/spades_assembly.fna \
                spades_assembly/spades_assembly.m1000.fna

    seqmagick convert --min-length 1000 \
                      idbaud_assembly/idbaud_assembly.fna \
                      idbaud_assembly/idbaud_assembly.m1000.fna
    grep -c '>' idbaud_assembly/idbaud_assembly.fna \
                idbaud_assembly/idbaud_assembly.m1000.fna
    ```

!!! circle-check "Terminal output"

    === "`SPAdes`"

        ```
        spades_assembly/spades_assembly.fna:1327
        spades_assembly/spades_assembly.m1000.fna:934
        ```

    === "`IDBA-UD`"

        ```
        idbaud_assembly/idbaud_assembly.fna:5056
        idbaud_assembly/idbaud_assembly.m1000.fna:1996
        ```

If you have your own assemblies and you want to try inspect them in the same way, try that now. Note that the file names will be slightly different to the files provided above. If you followed the exact commands in the previous exercise, you can use the following commands.

!!! terminal "code"

    ```bash
    seqmagick convert \
        --min-length 1000 \
        ../3.assembly/my_spades_assembly/scaffolds.fasta my_spades_assembly.m1000.fna

    seqmagick convert \
        --min-length 1000 \
        ../3.assembly/my_idbaud_assembly/scaffold.fa my_idbaud_assembly.m1000.fna
    ```

!!! note "Choice of software: sequence file manipulation"

    The tool `seqtk` is also available on NeSI and performs many of the same functions as `seqmagick`. My choice of `seqmagick` is mostly cosmetic as the parameter names are more explicit so it's easier to understand what's happening in a command when I look back at my log files. Regardless of which tool you prefer, we strongly recommend getting familiar with either `seqtk` or `seqmagick` as both perform a lot of common FASTA and FASTQ file manipulations.

As we can see here, the `SPAdes` assembly has completed with fewer contigs assembled than the `IDBA-UD`, both in terms of total contigs assembled and contigs above the 1,000 bp size. This doesn't tell us a lot though - has `SPAdes` managed to assemble fewer reads, or has it managed to assemble the sequences into longer (and hence fewer) contigs? We can check this by looking at the N50/L50 (see more information about this statistic [here](https://www.molecularecologist.com/2017/03/29/whats-n50/)) of the assembly using the `stats.sh` tool that comes with `BBMap`.

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

    ```
    A       C       G       T       N       IUPAC   Other   GC      GC_stdev
    0.2541  0.2475  0.2453  0.2531  0.0018  0.0000  0.0000  0.4928  0.0958

    Main genome scaffold total:             934
    Main genome contig total:               2703
    Main genome scaffold sequence total:    34.293 MB
    Main genome contig sequence total:      34.231 MB       0.182% gap
    Main genome scaffold N/L50:             53/160.826 KB
    Main genome contig N/L50:               107/72.909 KB
    Main genome scaffold N/L90:             302/15.325 KB
    Main genome contig N/L90:               812/4.643 KB
    Max scaffold length:                    1.222 MB
    Max contig length:                      1.045 MB
    Number of scaffolds > 50 KB:            152
    % main genome in scaffolds > 50 KB:     77.10%


    Minimum         Number          Number          Total           Total           Scaffold
    Scaffold        of              of              Scaffold        Contig          Contig  
    Length          Scaffolds       Contigs         Length          Length          Coverage
    --------        --------------  --------------  --------------  --------------  --------
        All                    934           2,703      34,293,018      34,230,627    99.82%
       1 KB                    934           2,703      34,293,018      34,230,627    99.82%
     2.5 KB                    744           2,447      33,967,189      33,909,776    99.83%
       5 KB                    581           2,131      33,379,814      33,328,450    99.85%
      10 KB                    394           1,585      32,031,121      31,994,989    99.89%
      25 KB                    236             916      29,564,265      29,545,757    99.94%
      50 KB                    152             599      26,441,339      26,429,391    99.95%
     100 KB                     92             415      22,238,822      22,230,143    99.96%
     250 KB                     31             153      12,606,774      12,603,418    99.97%
     500 KB                      6              32       4,821,355       4,821,095    99.99%
       1 MB                      1               2       1,221,548       1,221,538   100.00%
    ```

!!! danger "N50 and L50 in `BBMap`"

    Unfortunately, the N50 and L50 values generated by `stats.sh` are switched. N50 should be a length and L50 should be a count. The results table below shows the corrected values based on `stats.sh` outputs.

But what we can highlight here is that the statistics for the `SPAdes` assembly, with short contigs removed, yielded an N50 of 72.5 kbp at the contig level. We will now compute those same statistics from the other assembly options.

!!! terminal "code"

    ```bash
    stats.sh in=spades_assembly/spades_assembly.m1000.fna
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
