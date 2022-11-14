# Gene annotation (part 2) and coverage calculation

!!! info "Objectives"

    * [Gene prediction and annotation with `DRAM`](#gene-prediction-and-annotation-with-dram-distilled-and-refined-annotation-of-metabolism)
    * [Annotation of the MAGs with `DRAM`](#annotation-of-the-mags-with-dram)
    * [Calculate per-sample coverage stats for prokaryotic bins and viral contigs](#calculate-per-sample-coverage-stats-of-the-filtered-prokaryote-bins-and-viral-contigs)
    * [Select initial goal](#select-initial-goal)

---

### Gene prediction and annotation with *DRAM* (Distilled and Refined Annotation of Metabolism) 

[DRAM](http://dx.doi.org/10.1093/nar/gkaa621) is a tool designed to profile microbial (meta)genomes for metabolisms known to impact ecosystem functions across biomes. `DRAM` annotates MAGs and viral contigs using KEGG (if provided by user), UniRef90, PFAM, CAZy, dbCAN, RefSeq viral, VOGDB (Virus Orthologous Groups), and the MEROPS peptidase database. It is also highly customizable to other custom user databases. 

`DRAM` only uses assembly-derived *fastA* files input by the user. These input files may come from unbinned data (metagenome contig or scaffold files) or genome-resolved data from one or many organisms (isolate genomes, single-amplified genome (SAGs), MAGs).

`DRAM` is run in two stages: annotation and distillation. 

![](https://github.com/mcastudillo/MAG-annotation-with-DRAM/blob/main/figures/DRAM_workflow.png)

#### Annotation

The first step in `DRAM` is to annotate genes by assigning database identifiers to genes. Short contigs (default < 2,500 bp) are initially removed. Then, `Prodigal` is used to detect open reading frames (ORFs) and to predict their amino acid sequences. Next, `DRAM` searches all amino acid sequences against multiple databases, providing a single *Raw* output. When gene annotation is complete, all results are merged in a single tab-delimited annotation table, including the best hit for each database for user comparison. 

#### Distillation 

After genome annotation, a distill step follows with the aim to curate these annotations into useful functional categories, creating genome statistics and metabolism summary files, which are stored in the *Distillate* output. The genome statistics provides most genome quality information required for [MIMAG](https://www.nature.com/articles/nbt.3893) standards, including `GTDB-tk` and `checkM` information if provided by the user. The summarised metabolism table includes the number of genes with specific metabolic function identifiers (KO, CAZY ID, etc) for each genome, with information obtained from multiple databases. The *Distillate* output is then further distilled into the *Product*, an html file displaying a heatmap, as well as the corresponding data table. We will investigate all these files later on.  

---

### Annotation of the MAGs with *DRAM*

Beyond annotation, `DRAM` aims to be a data compiler. For that reason, output files from both `CheckM` and `GTDB_tk` steps can be input to `DRAM` to provide both taxonomy and genome quality information of the MAGs. 

#### *DRAM* input files

For these exercises, we have copied the relevant input files into the folder `10.gene_annotation/DRAM_input_files/`. `gtdbtk.bac120.classification_pplacer.tsv` was taken from the earlier `8.coverage_and_taxonomy/gtdbtk_out/` outputs, and `checkm.txt` from the result of re-running `CheckM` on the final refined filtered bins in `6.bin_refinement/filtered_bins`.

Navigate to the `10.gene_annotation/` folder

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/
```

The `CheckM` output file (`checkm.txt`) can be input as it is. However, in order to use the file with the `gtdb_tk` taxonomy (`gtdbtk.bac120.classification_pplacer.tsv`) we should modify it first to include column headers 'bin_id' and 'classification'

First, take a quick look at the current format of the taxonomy file using `less`

```bash
less DRAM_input_files/gtdbtk.bac120.classification_pplacer.tsv

#bin_3.filtered  d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__
#bin_8.filtered  d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Synechococcales;f__Cyanobiaceae;g__Prochlorococcus_C;s__
#bin_2.filtered  d__Bacteria;p__Planctomycetota;c__Brocadiae;o__Brocadiales;f__Brocadiaceae;g__;s__
#bin_5.filtered  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__
#bin_9.filtered  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__
#bin_4.filtered  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Nitrosomonadaceae;g__Nitrosomonas;s__
#bin_7.filtered  d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Nitrobacter;s__
#bin_0.filtered  d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Arcobacter;s__
#bin_1.filtered  d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Nautiliales;f__Nautiliaceae;g__;s__
#bin_6.filtered  d__Bacteria;p__Desulfobacterota_A;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio;s__
```

Now, use `sed` to add the headers as a new line. The `^` character indicates to make the additions at the start of the line (in this case, the first line). `\t` and `\n` represent a tab space and a newline character, respectively.

```bash
sed -i '1s/^/bin_id\tclassification\n/' DRAM_input_files/gtdbtk.bac120.classification_pplacer.tsv


less DRAM_input_files/gtdbtk.bac120.classification_pplacer.tsv

#bin_id  classification
#bin_3.filtered  d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__
#bin_8.filtered  d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Synechococcales;f__Cyanobiaceae;g__Prochlorococcus_C;s__
#bin_2.filtered  d__Bacteria;p__Planctomycetota;c__Brocadiae;o__Brocadiales;f__Brocadiaceae;g__;s__
#bin_5.filtered  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__
#bin_9.filtered  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__
#bin_4.filtered  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Nitrosomonadaceae;g__Nitrosomonas;s__
#bin_7.filtered  d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Nitrobacter;s__
#bin_0.filtered  d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Arcobacter;s__
#bin_1.filtered  d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Nautiliales;f__Nautiliaceae;g__;s__
#bin_6.filtered  d__Bacteria;p__Desulfobacterota_A;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio;s__

```

#### *DRAM* annotation

In default annotation mode, `DRAM` only requires as input the directory containing all the bins we would like to annotate in *fastA* format (either .fa or .fna). There are few parameters that can be modified if not using the default mode. Once the annotation step is complete, the mode `distill` is used to summarise the obtained results. 

*NOTE: due to the increased memory requirements, UniRef90 database is not default and the flag `–use_uniref` should be specified in order to search amino acid sequences against UniRef90. In this exercise, due to memory and time constraints, we won't be using the UniRef90 database.*

We will start by making sure `DRAM` is loaded properly.

```sh
module purge
module load DRAM/1.3.5-Miniconda3

DRAM.py --help

# usage: DRAM.py [-h] {annotate,annotate_genes,distill,strainer,neighborhoods,merge_annotations} ...
# 
# positional arguments:
#   {annotate,annotate_genes,distill,strainer,neighborhoods,merge_annotations}
#     annotate            Annotate genomes/contigs/bins/MAGs
#     annotate_genes      Annotate already called genes, limited functionality compared to annotate
#     distill             Summarize metabolic content of annotated genomes
#     strainer            Strain annotations down to genes of interest
#     neighborhoods       Find neighborhoods around genes of interest
#     merge_annotations   Merge multiple annotations to one larger set
# 
# options:
#   -h, --help            show this help message and exit
```

#### Submitting *DRAM* annotation as a slurm job

To run this exercise we first need to set up a slurm job. We will use the results for tomorrow's distillation step. 

Create a new script

```bash
nano dram_annnotation.sl
```

Paste in the script (update all of the cases of `<YOUR FOLDER>`)

```sh
#!/bin/bash -e

#SBATCH --account       nesi02659
#SBATCH --job-name      DRAM_annotation
#SBATCH --res           SummerSchool
#SBATCH --time          5:00:00
#SBATCH --mem           20Gb
#SBATCH --cpus-per-task 24
#SBATCH --error         slurm-DRAM_annot.%A-%a.err 
#SBATCH --output        slurm-DRAM_annot.%A-%a.out 


module purge
module load DRAM/1.3.5-Miniconda3

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/

DRAM.py annotate -i 'predictions/*.filtered.fna' \
--checkm_quality DRAM_input_files/checkm.txt \
--gtdb_taxonomy DRAM_input_files/gtdbtk.bac120.classification_pplacer.tsv \
-o annotation_dram --threads $SLURM_CPUS_PER_TASK
```

Submit the job

```bash
sbatch dram_annnotation.sl
```

The program will take 4-4.5 hours to run, so we will submit the jobs and inspect the results tomorrow morning.

---
### Calculate per-sample coverage stats of the filtered prokaryote bins and viral contigs

One of the first questions we often ask when studying the ecology of a system is: What are the pattens of abundance and distribution of taxa across the different samples? With bins of metagenome-assembled genome (MAG) data, we can investigate this by mapping the quality-filtered unassembled reads back to the refined bins to then generate coverage profiles. Genomes in higher abundance in a sample will contribute more genomic sequence to the metagenome, and so the average depth of sequencing coverage for each of the different genomes provides a proxy for abundance in each sample. 

As per the preparation step at the start of the binning process, we can do this using read mapping tools such as `Bowtie`, `Bowtie2`, and `BBMap`. Here we will follow the same steps as before using `Bowtie2`, `samtools`, and `MetaBAT`'s `jgi_summarize_bam_contig_depths`, but this time inputting our refined filtered bins. 

These exercises will take place in the `8.coverage_and_taxonomy/` folder. Our final filtered refined bins from the previous bin refinement exercise have been copied to the `8.coverage_and_taxonomy/filtered_bins/` folder.

First, concatenate the bin data into a single file to then use to generate an index for the read mapper.

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/

cat filtered_bins/*.fna > filtered_bins.fna
```

Now build the index for `Bowtie2` using the concatenated bin data. We will also make a new directory `bin_coverage/` to store the index and read mapping output into.

```bash
mkdir -p bin_coverage/

# Load Bowtie2
module purge
module load Bowtie2/2.4.5-GCC-11.3.0

# Build Bowtie2 index
bowtie2-build filtered_bins.fna bin_coverage/bw_bins
```

Map the quality-filtered reads (from `../3.assembly/`) to the index using `Bowtie2`, and sort and convert to `.bam` format via `samtools`.

Create a new script

```bash
nano mapping_bins.sl
```

Paste in the script (replacing `<YOUR FOLDER>`)

```bash 
#!/bin/bash -e

#SBATCH --account       nesi02659
#SBATCH --res           SummerSchool
#SBATCH --job-name      mapping_filtered_bins
#SBATCH --time          00:05:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 10
#SBATCH --error         mapping_filtered_bins.err
#SBATCH --output        mapping_filtered_bins.out

module purge
module load Bowtie2/2.4.5-GCC-11.3.0 SAMtools/1.15.1-GCC-11.3.0

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/

# Step 1
for i in sample1 sample2 sample3 sample4;
do

  # Step 2
  bowtie2 --minins 200 --maxins 800 --threads 10 --sensitive \
          -x bin_coverage/bw_bins \
          -1 ../3.assembly/${i}_R1.fastq.gz -2 ../3.assembly/${i}_R2.fastq.gz \
          -S bin_coverage/${i}.sam

  # Step 3
  samtools sort -@ 10 -o bin_coverage/${i}.bam bin_coverage/${i}.sam

done
```

Submit the script

```bash
sbatch mapping_bins.sl
```

Finally, generate the per-sample coverage table for each contig in each bin via `MetaBAT`'s `jgi_summarize_bam_contig_depths`.

```bash 
# Load MetaBAT
module load MetaBAT/2.15-GCC-11.3.0

# calculate coverage table
jgi_summarize_bam_contig_depths --outputDepth bins_cov_table.txt bin_coverage/sample*.bam
```

The coverage table will be generated as `bins_cov_table.txt`. As before, the key columns of interest are the `contigName`, and each `sample[1-n].bam` column.

!!! note "Note"
    Here we are generating a per-sample table of coverage values for **each contig** within each bin. To get per-sample coverage of **each bin** as a whole, we will need to generate average coverage values based on all contigs contained within each bin. We will do this in `R` during our data visualisation exercises on day 4 of the workshop, leveraging the fact that we added bin IDs to the sequence headers.*

### Calculate per-sample coverage stats of viral contigs

Here we can follow the same steps as outlined above for the bin data, but with a concatenated *fastA* file of viral contigs. 

To quickly recap: 

* In previous exercises, we first used `VIBRANT` to identify viral contigs from the assembled reads, generating a new fasta file of viral contigs: `spades_assembly.m1000.phages_combined.fna` 
* We then processed this file using `CheckV` to generate quality information for each contig, and to further trim any retained (prokaryote) sequence on the ends of prophage contigs. 

The resultant *fasta* files generated by `CheckV` (`proviruses.fna` and `viruses.fna`) have been copied to to the `8.coverage_and_taxonomy/checkv` folder for use in this exercise.

!!! note "Note"
    Due to the rapid mutation rates of viruses, with full data sets it will likely be preferable to first further reduce viral contigs down based on a percentage-identity threshold using a tool such as `BBMap`'s `dedupe.sh`. This would be a necessary step in cases where you had opted for generating multiple individual assemblies or mini-co-assemblies (and would be comparable to the use of a tool like `dRep` for prokaryote data), but may still be useful even in the case of single co-assemblies incorporating all samples.*

We will first need to concatenate these files together.

```bash
cat checkv/proviruses.fna checkv/viruses.fna > checkv_combined.fna
```

Now build the index for `Bowtie2` using the concatenated viral contig data. We will also make a new directory `viruses_coverage/` to store the index and read mapping output into.

```bash
mkdir -p viruses_coverage/

# Load Bowtie2
module load Bowtie2/2.4.5-GCC-11.3.0

# Build Bowtie2 index
bowtie2-build checkv_combined.fna viruses_coverage/bw_viruses
```

Map the quality-filtered reads (from `../3.assembly/`) to the index using `Bowtie2`, and sort and convert to `.bam` format via `samtools`.

Create a new script

```bash
nano mapping_viruses.sl
```

Paste in the script (replacing `<YOUR FOLDER>`)

```bash 
#!/bin/bash -e

#SBATCH --account       nesi02659
#SBATCH --res           SummerSchool
#SBATCH --job-name      mapping_filtered_viruses
#SBATCH --time          00:05:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 10
#SBATCH --error         mapping_filtered_viruses.err
#SBATCH --output        mapping_filtered_viruses.out

module purge
module load Bowtie2/2.4.5-GCC-11.3.0 SAMtools/1.15.1-GCC-11.3.0

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/

# Step 1
for i in sample1 sample2 sample3 sample4;
do

  # Step 2
  bowtie2 --minins 200 --maxins 800 --threads 10 --sensitive \
          -x viruses_coverage/bw_viruses \
          -1 ../3.assembly/${i}_R1.fastq.gz -2 ../3.assembly/${i}_R2.fastq.gz \
          -S viruses_coverage/${i}.sam

  # Step 3
  samtools sort -@ 10 -o viruses_coverage/${i}.bam viruses_coverage/${i}.sam

done
```

Run the script

```bash
sbatch mapping_viruses.sl
```

Finally, generate the per-sample coverage table for each viral contig via `MetaBAT`'s `jgi_summarize_bam_contig_depths`.

```bash 
# Load MetaBAT
module load MetaBAT/2.15-GCC-11.3.0

# calculate coverage table
jgi_summarize_bam_contig_depths --outputDepth viruses_cov_table.txt viruses_coverage/sample*.bam
```

The coverage table will be generated as `viruses_cov_table.txt`. As before, the key columns of interest are the `contigName`, and each `sample[1-n].bam` column.

!!! note "Note"
    Unlike the prokaryote data, we have not used a binning process on the viral contigs (since many of the binning tools use hallmark characteristics of prokaryotes in the binning process). Here, `viruses_cov_table.txt` is the final coverage table. This can be combined with `CheckV` quality and completeness metrics to, for example, examine the coverage profiles of only those viral contigs considered to be "High-quality" or "Complete".* 

#### Normalising coverage values

Having generated per-sample coverage values, it is usually necessary to also normalise these values across samples of differing sequencing depth. In this case, the mock metagenome data we have been working with are already of equal depth, and so this is an unnecessary step for the purposes of this workshop. 

For an example of one way in which the `cov_table.txt` output generated by `jgi_summarize_bam_contig_depths` above could then be normalised based on average library size, see the [Normalise per-sample coverage Appendix](https://genomicsaotearoa.github.io/metagenomics_summer_school/resources/3_APPENDIX_ex11_Normalise_coverage_example/).

---

### Select initial goal

It is now time to select the goals to investigate the genomes you have been working with. We ask you to select one of the following goals:

!!! quote ""

    1. Denitrification (Nitrate or nitrite to nitrogen)
    2. Ammonia oxidation (Ammonia to nitrite or nitrate)
    3. Anammox (Ammonia and nitrite to nitrogen)
    4. Sulfur oxidation (SOX pathway, thiosulfate to sulfate)
    5. Sulfur reduction (DSR pathway, sulfate to sulfide)
    6. Photosynthetic carbon fixation
    7. Non-photosynthetic carbon fixation (Reverse TCA or Wood-Ljundahl)
    8. Non-polar flagella expression due to a chromosomal deletion
    9. Plasmid-encoded antibiotic resistance
    10. Aerobic (versus anaerobic) metabolism

Depending on what you are looking for, you will either be trying to find gene(s) of relevance to a particular functional pathway, or the omission of genes that might be critical in function. In either case, make sure to use the taxonomy of each MAG to determine whether it is likely to be a worthwhile candidate for exploration, as some of these traits are quite restricted in terms of which organisms carry them.

To conduct this exersise, you should use the information generated with ```DRAM``` as well as the annotation files we created previously that will be available in the directory ```10.gene_annotation/gene_annotations```. 

Please note that we have also provided further annotation files within the directory ```10.gene_annotation/example_annotation_tables``` that contain information obtained after annotating the MAGs against additional databases (UniProt, UniRef100, KEGG, PFAM and TIGRfam). These example files can also be downloaded from [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/example_annotation_tables.zip). These files were created by using an in-house python script designed to aggregate different annotations and as part of the environmental metagenomics worflow followed in Handley's lab. Information about using this script as well as the script is available [here](https://github.com/GenomicsAotearoa/environmental_metagenomics/blob/master/metagenomic_annotation/3.aggregation.md)  
