# Introduction to binning

!!! info "Objectives"

    * [Remove short contigs from the data set](#remove-short-contigs-from-the-data-set)
    * [Obtain coverage profiles for assembled contigs via read mapping](#obtain-coverage-profiles-for-assembled-contigs-via-read-mapping)
    * [*Optional:* Read mapping using a slurm array](#optional-read-mapping-using-an-array)

---

### Remove short contigs from the data set

Ideally, we do not want to be creating bins from all of the assembled contigs, as there is often a long tail of contigs which are only several $k$-mers long. These have little biological meaning, as they are too short for robust gene annotation, and they can introduce a significant degree of noise in the clustering algorithms used for binning. We therefore identify a suitable threshold for a minimum length of contigs to be considered for binning.

We have already done this in the [previous exercise](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day1/ex5_evaluating_assemblies.md) so we could either use the existing filtering at 1,000 bp in length, or move to something stricter. Most binning tools have a default cut-off for minimum contig size - `MetaBAT` uses a default minimum of 2,500 bp, and recommends at least 1,500 bp. By contrast, `MaxBin` sets the minimum length at 1,000 bp.

---

### Obtain coverage profiles for assembled contigs via read mapping

Binning is done using a combination of information encoded in the *composition* and *coverage* of the assembled contigs. *Composition* refers to *k*-mer (usually tetranucleotide) frequency profiles of the contigs, which are generally conserved within a genome. By contrast, *coverage* is a reflection of the abundance of the contigs in the assembly. Organisms which are more abundant will contribute more genomic material to the metagenome, and hence their DNA will be, on average, more abundant in the sample. When binning, we can look for pieces of DNA which are not assembled together, but have similar *composition* and occur at approximately equal abundances in the sample to identify contigs which likely originate in the same genome.

The composition of the contigs is calculated by the binning tool at run time, but to obtain coverage information we must map our unassembled reads from each sample against the assembly to generate the differential abundance profiles for each contig. This is achieved using `bowtie2` to map the reads against the assembly, then `samtools` to sort and compress the resulting file.

#### Creating a mapping index

Before we can map reads, we need to create a `bowtie2` index file from the assembly, for use in read mapping. Navigate into the `5.binning/` folder to begin.

```bash
module load Bowtie2/2.4.5-GCC-11.3.0

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/

bowtie2-build spades_assembly/spades_assembly.m1000.fna spades_assembly/bw_spades
```

If you look inside the `spades_assembly/` folder you will now see the following:

```bash
ls spades_assembly/
# bw_spades.1.bt2  bw_spades.3.bt2  bw_spades.rev.1.bt2  spades_assembly.fna
# bw_spades.2.bt2  bw_spades.4.bt2  bw_spades.rev.2.bt2  spades_assembly.m1000.fna
```

These files ending in *.bt2* are the index files for `bowtie2`, and are specific to this tool. If you wish to map using an alternate tool (for example `bowtie` or `BBMap`) you will need to create index/database files using these programs.

Generally speaking, we don't need to know the names of the index files, as they are simply referred to be the output name we specified (*bw_spades*) when running `bowtie2`.

#### Mapping the reads

We will create a slurm script to perform the mapping steps, as these benefit greatly from the multithreaded capacity of NeSI and we will use a *for loop* to iterate over each set of reads to simplify our script.

The full script is provided here, and we will discuss it below.

Open a new script using nano:

```bash 
nano spades_mapping.sl
```

Paste in the following script. Remember to replace <YOUR FOLDER> with your own folder.

```bash
#!/bin/bash -e
#SBATCH --account nesi02659
#SBATCH --job-name spades_mapping
#SBATCH --res SummerSchool
#SBATCH --time 00:05:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --error spades_mapping.err
#SBATCH --output spades_mapping.out

module purge
module load Bowtie2/2.3.5-GCC-7.4.0 SAMtools/1.8-gimkl-2018b

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/

# Step 1
for i in sample1 sample2 sample3 sample4;
do

  # Step 2
  bowtie2 --minins 200 --maxins 800 --threads 10 --sensitive \
          -x spades_assembly/bw_spades \
          -1 ../3.assembly/${i}_R1.fastq.gz -2 ../3.assembly/${i}_R2.fastq.gz \
          -S ${i}.sam

  # Step 3
  samtools sort -@ 10 -o ${i}.bam ${i}.sam

done
```

Now run the script using `sbatch`

```bash
sbatch spades_mapping.sl
```

#### Step 1 - Loop through the sample files

Since we just want to perform the same set of operations on each file, we can use a *for loop* to repeat each operation on a set of files. The structure of the loop, and use of variables was covered on the first day.

For large sets of files, it can be beneficial to use a slurm *array* to send the jobs out to different nodes and distribute the process across many independent jobs. An example of how we could modify the above script is given at the bottom of this exercise, but is not necessary for the purposes of this workshop.

#### Step 2 - Map the reads using *bowtie2*

This is performed using the following parameters

|Parameter|Function|
|:---|:---|
|**--minins ...**|Minimum insert size, determines the minimum distance between the start of each read pair|
|**--maxins ...**|Maximum insert size, determines the maximum distance between the start of each read pair|
|**--threads ...**|Number of threads to use in read mapping|
|**--sensitive**|Specifies where we want to be positioned in the trade-off between speed and sensitivity. See the manual for more information|
|**-x ...**|The base name of our assembly index. Should be exactly the same as what was specified when running **bowtie2-build**|
|**-1 ...** / **-2 ...**|The forward and reverse read pairs to map to the assembly|
|**-S ...**|Name of the output file, to be written in *sam* format|

#### Step 3 - Sorting and compressing results

The default output format for most maping tools is the Sequence Alignment/Map (*sam*) format. This is a compact text representation of where each short read sits in the contigs. You can view this file using any text viewer, although owing to the file size `less` is a good idea.

Generally I wouldn't bother with this - there is a lot of information in here and unless you are looking to extract specific information from the alignment directly, this is just an intermediate file in our workflow. In order to save disk space, and prepare the file for downstream analysis we now perform two final steps:

1. Sort the mapping information
1. Compress the *sam* file into its binary equivalent, *bam*

Which is achieved with the following parameters

|Parameter|Function|
|:---|:---|
|**sort**|Subcommand for `samtools` to invoke the `sort` operation|
|**-@ ...**|Number of threads to use for sorting and compressing|
|**-o ...**|Output file name. When we specify the *bam* extension `samtools` automatically compresses the output|

Compressing the file to the *bam* format is an important step as when working with real data *sam* files can be massive and our storage capacity on NeSI is limited. It is also helpful to sort the mapping information so that reads mapped to a contig are listed in order of their start position. For example

```bash
# Unsorted reads
Ref: REFERENCECONTIG
Map: --------ECONT--
Map: REFE-----------
Map: --FERENCECO----
Map: -----------NTIG

# Sorted reads
Ref: REFERENCECONTIG
Map: REFE-----------
Map: --FERENCECO----
Map: --------ECONT--
Map: -----------NTIG
```

Reads will initially be mapped in an unsorted order, as they are added to the *sam* file in more or less the same order as they are encountered in the original *fastQ* files.

Sorting the mapping information is an important prerequisite for performing certain downstream processes. Not every tool we use requires reads to be sorted, but it can be frustrating having to debug the instances where read sorting matters, so we typically just get it done as soon as possible and then we don't have to worry about it again.

In newer versions of `samtools` we can perform the sorting and compressing in a single operation (as shown in the script above). For older versions of `samtools`, you may need to use a command of the following form.

```bash
samtools view -bS sample1.sam | samtools sort -o sample1.bam
```

---

### *Optional: Read mapping using an array*

If you have a large number of files to process, it might be worth using a slurm array to distribute you individual mapping jobs across many separate nodes. An example script for how to perform this is given below, although it will not be covered in this workshop.

Open a new script using nano:

```bash 
nano spades_mapping_array.sl
```

```bash
#!/bin/bash -e
#SBATCH --account nesi02659
#SBATCH --job-name spades_mapping_array
#SBATCH --res SummerSchool
#SBATCH --time 00:20:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --array 0-3
#SBATCH --cpus-per-task 10
#SBATCH --error spades_mapping_array.%j.err
#SBATCH --output spades_mapping_array.%j.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load Bowtie2/2.3.5-GCC-7.4.0 SAMtools/1.8-gimkl-2018b

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/

srun bowtie2-build spades_assembly/spades_assembly.m1000.fna spades_assembly/bw_spades

# Load the sample names into a bash array
samples=(sample1 sample2 sample3 sample4)

# Activate the srun command, using the SLURM_ARRAY_TASK_ID variable to
# identify which position in the `samples` array to use
srun bowtie2 --minins 200 --maxins 800 --threads 10 --sensitive -x spades_assembly/bw_spades \
             -1 ../3.assembly/${samples[ $SLURM_ARRAY_TASK_ID ]}_R1.fastq.gz \
             -2 ../3.assembly/${samples[ $SLURM_ARRAY_TASK_ID ]}_R2.fastq.gz \
             -S ${samples[ $SLURM_ARRAY_TASK_ID ]}.sam

srun samtools sort -o ${samples[ $SLURM_ARRAY_TASK_ID ]}.bam ${samples[ $SLURM_ARRAY_TASK_ID ]}.sam
```

submit the job to slurm

```bash
sbatch spades_mapping_array.sl
```
