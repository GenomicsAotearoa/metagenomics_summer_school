# APPENDIX (ex9): Generating input files for *VizBin* from *DAS_Tool* curated bins

The final bins that we obtained in the previous step (output from `DAS_Tool`) have been copied into `6.bin_refinement/dastool_out/_DASTool_bins/`

**1. Generalise bin naming and add bin IDs to sequence headers**

We will first modify the names of our bins to be simply numbered 1 to n bins. We will use a loop to do this, using the wildcard ( * ) to loop over all files in the _DASTool_bins folder, copying to the new example_data_unchopped/ folder and renaming as bin_[1-n].fna. The sed command then adds the bin ID to the start of sequence headers in each of the new bin files (this will be handy information to have in the sequence headers for downstream processing).

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/6.bin_refinement/

# Make a new directory the renamed bins
mkdir example_data_unchopped/

# Copy and rename bins into generic <bin_[0-n].fna> filenames
i=0
for file in dastool_out/_DASTool_bins/*;
do
    # Copy and rename bin file
    cp ${file} example_data_unchopped/bin_${i}.fna
    # extract bin ID
    binID=$(basename example_data_unchopped/bin_${i}.fna .fna)
    # Add bin ID to sequence headers
    sed -i -e "s/>/>${binID}_/g" example_data_unchopped/bin_${i}.fna
    # Increment i
    ((i+=1))
done
```

**2. Fragment contigs**

Using the `cut_up_fasta.py` script that comes with the binning tool `CONCOCT`, cut contigs into 20k fragments to add better density to the cluster.

```bash
# Load CONCONCT module
module load CONCOCT/1.0.0-gimkl-2018b-Python-2.7.16

# Make directory to add chopped bin data into
mkdir example_data_20k/

# loop over .fna files to generate chopped (fragmented) files using CONCONT's cut_up_fasta.py
for bin_file in example_data_unchopped/*;
do
    bin_name=$(basename ${bin_file} .fna)
    cut_up_fasta.py -c 20000 -o 0 --merge_last ${bin_file} > example_data_20k/${bin_name}.chopped.fna
done
```

**3. Concatenate fragmented bins**

Concatenate chopped bins into a single *fastA* file.

Concatenate chopped bins into a single `all_bins.fna` file to use as input for both subcontig read mapping via `Bowtie2` and visualisation via `VizBin`.

```bash
cat example_data_20k/*.fna > all_bins.fna
```

**4. Read mapping of subcontigs (fragmented contigs based on 20k length)**

**4a. Build mapping index**

Build `Bowtie2` mapping index based on the concatenated chopped bins.

```bash
mkdir -p read_mapping/

module load Bowtie2/2.3.5-GCC-7.4.0

bowtie2-build all_bins.fna read_mapping/bw_bins
```

**4b. Map sample reads to index**

Map quality filtered reads to the index using `Bowtie2`.

Example slurm script:

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J 6.bin_refinement_mapping
#SBATCH --time 00:05:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH -e 6.bin_refinement_mapping.err
#SBATCH -o 6.bin_refinement_mapping.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load Bowtie2/2.3.5-GCC-7.4.0 SAMtools/1.8-gimkl-2018b

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/6.bin_refinement/

# Step 1
for i in sample1 sample2 sample3 sample4;
do

  # Step 2
  bowtie2 --minins 200 --maxins 800 --threads 10 --sensitive \
          -x read_mapping/bw_bins \
          -1 ../3.assembly/${i}_R1.fastq.gz -2 ../3.assembly/${i}_R2.fastq.gz \
          -S read_mapping/${i}.sam

  # Step 3
  samtools sort -@ 10 -o read_mapping/${i}.bam read_mapping/${i}.sam

done
```

*Note: These settings are appropriate for this workshop's mock data. Full data sets will likely require considerably greater memory and time allocations.*

**5. Generate coverage table of subcontigs (contigs fragmented based on 20k length)**

Use `MetaBAT`'s `jgi_summarise_bam_contig_depths` to generate a coverage table.

```bash
# Load module
module load MetaBAT/2.13-GCC-7.4.0

# calculate coverage table
jgi_summarize_bam_contig_depths --outputDepth example_data_20k_cov.txt read_mapping/sample*.bam
```

**6. Generate annotation table for `VizBin`**

Using the chopped bin files (`example_data_20k/`) and the coverage table generated above (`example_data_20k_cov.txt`), we can use the following script to generate an annotation table in the format that `VizBin` is expecting. Note that here we are including columns for per-sub-contig coverage based on *sample1* (see note at the start of this exercise), label (bin ID), and length values, and storing this in `all_bins.sample1.vizbin.ann`.

What this script is doing is taking each fasta file and picking out the names of the contigs found in that file (bin). It is then looking for any coverage information for sample1 which is associated with that contig in the `example_data_20k_cov.txt` file, and adding that as a new column to the file. (If you wished to instead look based on sample2, you would need to modify the `cut` command in the line `sample1_cov=$(grep -P "${contigID}\t" example_data_20k_cov.txt | cut -f4)` accordingly (e.g. `cut -f6`).

```bash
# Set up annotation file headers
echo "coverage,label,length" > all_bins.sample1.vizbin.ann

# loop through bin .fna files
for bin_file in example_data_20k/*.fna; do
    # extract bin ID
    binID=$(basename ${bin_file} .fna)
    # loop through each sequence header in bin_file
    for header in `grep ">" ${bin_file}`; do
        contigID=$(echo ${header} | sed 's/>//g')
        # identify this line from the coverage table (example_data_20k_cov.txt), and extract contigLen (column 2) and coverage for sample1.bam (column 4)
        contigLen=$(grep -P "${contigID}\t" example_data_20k_cov.txt | cut -f2)
        sample1_cov=$(grep -P "${contigID}\t" example_data_20k_cov.txt | cut -f4)
        # Add to vizbin.ann file
        echo "${sample1_cov},${binID},${contigLen}" >> all_bins.sample1.vizbin.ann
    done
done
```

We now have the `all_bins.fna` and `all_bins.sample1.vizbin.ann` files that were provided at the [start of the VizBin exercise](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day2/ex9_refining_bins.md#prepare-input-files-for-vizbin).

---
