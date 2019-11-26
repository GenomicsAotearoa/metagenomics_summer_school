# Evaluating the assemblies

### Objectives

* Evaluate the resource consumption of various assemblies
* Evaluate the assemblies
* Future considerations

---

### Evaluate the resource consumption of various assemblies

Check to see if your jobs from last night have completed. If you have multiple jobs running or queued, the easiest way to check this is to simply run the **squeue* command from yesterday.

```bash
squeue -u <user name>

#  JOBID     USER ACCOUNT            NAME  ST REASON    START_TIME                TIME TIME_LEFT NODES CPUS
```

Since there are no jobs listed, either everything running has completes or failed. To get a list of all jobs we have run in the last day, we can use the **sacct** command. By default this will report all jobs for the day but we can add a parameter to tell the command to report all jobs run since the date we are specifying.

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

Each job has been broken up into several lines, but the main ones to keep an eye are the base JobID values, and the values suffixed with *.0*. The first of these references the complete job. The later (and any subsequent suffixes like *.1*, *.2*) are the individual steps in the script that were called with the **srun** command.

We can see here the time ellapsed for each job, and the number of CPU hours used during the run. If we want a more detailed breakdown of the job we can use the **seff** command

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

Here we see some of the same information, but we also get some information regarding how well our job used the resources we allocated to it. You can see here that my CPU and memory usage was not particularly efficient, in hindsight I could have request a lot less RAM and still had the job run to completion. CPU efficiency is harder to be certain of as it impacted by the behaviour of the program. For example, mapping tools like **bowtie** and **BBMap** can more or less use all of their threads, all of the time and achieve nearly 100% efficiency. More complicated processes, like those performed in **SPAdes** go through periods of multi-thread processing and periods of single-thread processing, drawing the average efficiency down.

---

### Evaluate the assemblies

Evaluating the quality of a raw metagenomic assembly is quite a tricky process. Since, by definition, our community is amixture of different organisms, the genomes from some of these organisms assemble better than those of others. Is is possible to have an assembly that looks 'bad' by tranditional metrics that still yields high-quality genomes from individual species, and the converse is also true. A few quick checks I recommend are to see how many contigs or scaffolds your data were assembled into, and then see how many contigs or scaffolds you have above a certain minimum length threshold. I usually use **seqmagick** for performing the length filtering, and then just count sequence numbers using **grep**.

```bash
module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

seqmagick convert 

# TO COMPLETE
/nesi/nobackup/ga02676/Metagenomics_summerschool
spades_assembly/scaffolds.fasta
idbaud_assembly/scaffold.fa
```



---

### Future considerations

#### Co-assembly vs single assemblies
