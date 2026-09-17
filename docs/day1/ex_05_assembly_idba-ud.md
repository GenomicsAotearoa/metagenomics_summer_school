# BONUS: Assembling contigs with IBDA-UD





All work for this exercise will occur in the `3.assembly/` directory.

---

## The standard input files for `SPAdes` and `IDBA-UD`

Although they both make use of the same types of data, both `SPAdes` and `IDBA-UD` have their own preferences for how sequence data is provided to them. To begin, we will look at the types of data accepted by `SPAdes`:

!!! terminal-2 "Navigate to working directory"

    ```bash
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/3.assembly
    ```

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"




Awkwardly, while `SPAdes` accepts multiple input libraries (i.e. samples) in a single assembly, this behaviour does not work with the `-meta` flag enabled, which is needed in our example to activate metagenome assembly mode. We therefore need to concatenate our four individual samples together ready for sequencing.

!!! terminal-2 "Concatenate samples by read direction"

    ```bash
    cat sample1_R1.fastq.gz sample2_R1.fastq.gz sample3_R1.fastq.gz sample4_R1.fastq.gz > for_spades_R1.fq.gz
    cat sample1_R2.fastq.gz sample2_R2.fastq.gz sample3_R2.fastq.gz sample4_R2.fastq.gz > for_spades_R2.fq.gz
    ```

Note that these FASTQ files are compressed, yet we can concatenate them together with the `cat` command regardless. This is a nice feature of `.gz` files that is handy to remember.

By contrast, what does `IDBA-UD` accept?

!!! terminal-2 "Load `IDBA-UD` and check parameters"

    ```bash
    # Load module
    module purge
    module load IDBA-UD/1.1.3-GCC-11.3.0

    # Check parameters
    idba_ud --help
    ```

??? circle-check "`IDBA-UD` parameters"

    ```
      -o, --out arg (=out)                   output directory
      -r, --read arg                         fasta read file (<=128)
          --read_level_2 arg                 paired-end reads fasta for second level scaffolds
          --read_level_3 arg                 paired-end reads fasta for third level scaffolds
          --read_level_4 arg                 paired-end reads fasta for fourth level scaffolds
          --read_level_5 arg                 paired-end reads fasta for fifth level scaffolds
      -l, --long_read arg                    fasta long read file (>128)
          --mink arg (=20)                   minimum k value (<=124)
          --maxk arg (=100)                  maximum k value (<=124)
          --step arg (=20)                   increment of k-mer of each iteration
          --inner_mink arg (=10)             inner minimum k value
          --inner_step arg (=5)              inner increment of k-mer
          --prefix arg (=3)                  prefix length used to build sub k-mer table
          --min_count arg (=2)               minimum multiplicity for filtering k-mer when building the graph
          --min_support arg (=1)             minimum supoort in each iteration
          --num_threads arg (=0)             number of threads
          --seed_kmer arg (=30)              seed kmer size for alignment
          --min_contig arg (=200)            minimum size of contig
          --similar arg (=0.95)              similarity for alignment
          --max_mismatch arg (=3)            max mismatch of error correction
          --min_pairs arg (=3)               minimum number of pairs
          --no_bubble                        do not merge bubble
          --no_local                         do not use local assembly
          --no_coverage                      do not iterate on coverage
          --no_correct                       do not do correction
          --pre_correction                   perform pre-correction before assembly
    ```

'Short' or 'long' reads, and only a single file for each. This means that if we want to assemble our community data using `IDBA-UD` we will need to pool the paired-end data into a single, interleaved FASTA file. Interleaved means that instead of having a pair of files that contain the separate forward and reverse sequences, the read pairs are in a single file in alternating order. For example

```bash
# Paired-end file, forward
>read1_1
...
>read2_1
...

# Paired-end file, reverse
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

!!! terminal-2 "Decompress, interleave, and then concatenate sequence files"

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


## Setting the *k*-mer size


For `IDBA-UD`, we can select the *k*-mer size using

!!! terminal "code"

    ```bash
    idba_ud --mink 21 --maxk 121 --step 22
    ```

Unlike `SPAdes`, we do not have fine-scale control over the *k*-mer sizes used in the assembly. We instead provide `IDBA-UD` with the first and last *k*-mer size to use, then specify the increment to use between these. In either case, it is important that we are always assembling using a *k*-mer of odd lengths in order to avoid the creation of palindromic *k*-mers.

### Specifying the number of threads

This is simple in either assembler:

!!! terminal "code"

    ```bash

    idba_ud --num_threads 20 ...
    ```

The only thing to keep in mind is that these tools have different default behaviour. If no thread count is specified by the user, `SPAdes` will assemble with 16 threads. `IDBA-UD` will use all available threads, which can be problematic if you are using a shared compute environment that does not use a resource management system like slurm.

### Setting a memory limit

No such parameter exists in `IDBA-UD`, but it requires far less RAM than `SPAdes`, so you are less likely to need it.

---

   

### Submitting an `IDBA-UD` job to NeSI using slurm

!!! terminal-2 "Create a new slurm script using `nano` to run an equivalent assembly with `IDBA-UD`"

    ```bash
    nano idbaud_assembly.sl
    ```

Paste or type in the following:

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      idbaud_assembly
    #SBATCH --partition     milan
    #SBATCH --time          00:20:00
    #SBATCH --mem           4GB
    #SBATCH --cpus-per-task 12
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Prepare modules
    module purge
    module load IDBA-UD/1.1.3-GCC-11.3.0
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/3.assembly
    
    # Run IDBA-UD
    idba_ud --num_threads $SLURM_CPUS_PER_TASK --mink 33 --maxk 99 --step 22 \
            -r for_idba.fna -o idbaud_assembly/
    ```

!!! terminal-2 "Submit the script as a slurm job"

    ```bash
    sbatch idbaud_assembly.sl
    ```

When your job starts running, files with suffixes `.err` and `.out` will be created in the directory from where you submitted your job. These files will have have your job name and job identification number as file names.

---
