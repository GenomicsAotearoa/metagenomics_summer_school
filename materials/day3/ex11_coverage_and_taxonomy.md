# Per-sample coverage stats and assigning taxonomy

### Objectives

* Calculate per-sample coverage stats of the filtered prokaryote bins
* Calculate per-sample coverage stats of the viral contigs output by `VIBRANT`
* Assign taxonomy to the refined bins
* Overview of using `vContact2` to predict taxonomy of viral contigs

---

### Calculate per-sample coverage stats of the filtered prokaryote bins

*WIP (MH)*

---

### Calculate per-sample coverage stats of viral contigs

*WIP (MH)*

---

### Assign taxonomy to the refined bins

It is always valuable to know the taxonomy of our binned MAGs, so that we can link them to the wider scientific literature. In order to do this, there are a few different options available to us:

1. Extract 16S rRNA gene sequences from the MAGs and classify them
1. Annotate each gene in the MAG and take the consensus taxonomy
1. Use a profiling tool like `Kraken`, which matches pieces of DNA to a reference database using *k*-mer searches
1. Identify a core set of genes in the MAG, and use these to compute a species phylogeny

For this exercise, we will use the last option in the list, making use of the `GTDB-TK` software (available on [github](https://github.com/Ecogenomics/GTDBTk)) to automatically identify a set of highly conserved, single copy marker genes which are diagnostic of the bacterial (120 markers) and archaeal (122 markers) lineages. Briefly, `GTDB-TK` will perform the following steps on a set of bins.

1. Attempt to identify a set of 120 bacterial marker genes, and 122 archaeal marker genes in each MAG.
1. Based on the recovered numbers, identify which domain is a more likely assignment for each MAG
1. Create a concatenated alignment of the domain-specific marker genes, spanning approximately 41,000 amino acid positions
1. Filter the alignment down to approximately 5,000 informative sites
1. Insert each MAG into a reference tree create from type material and published MAGs
1. Scale the branch lengths of the resulting tree, as described in [Parks et al.](https://www.ncbi.nlm.nih.gov/pubmed/30148503), to identify an appropriate rank to each branch event in the tree
1. Calculate ANI and AAI statistics between each MAG and its nearest neighbours in the tree
1. Report the resulting taxonomic assignment, and gene alignment

This can all be achieved in a single command, although it must be performed through a slurm script due to the high memory requirements of the process.

```bash
#!/bin/bash
#SBATCH -A nesi02659
#SBATCH -J gtdbtk_test
#SBATCH --partition ga_bigmem
#SBATCH --res SummerSchool
#SBATCH --time 00:30:00
#SBATCH --mem 140GB
#SBATCH --cpus-per-task 10
#SBATCH -e gtdbtk_test.err
#SBATCH -o gtdbtk_test.out

module load GTDB-Tk/0.2.2-gimkl-2018b-Python-2.7.16

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/6.bin_refinement/

gtdbtk classify_wf -x fna --cpus 10 --genome_dir filtered_bins/ --out_dir gtdbtk_out/
```

As usual, lets look at the parameters here

|Parameter|Function|
|:---|:---|
|**classify_wf**|Specifies the sub-workflow from `GTDB-TK` that we wish to use|
|**-x ...**|Specify the file extension for MAGs within our input directory.<br>Default is *.fna*, but it's always good practice to specify it anyway|
|**--cpus ...**|Number of CPUs to use when finding marker genes, and performing tree insertion operations|
|**--genome_dir ...**|Input directory containing MAGs as individual *fastA* files|
|**--out_dir ...**|Output directory to write the final set of files|

Before submitting your job, think carefully about which set of MAGs you want to classify. You could either use the raw `DAS_Tool` outputs in the `dastool_out/_DASTool_bins/` folder, the renamed set of bins in the `example_data_unchopped/` folder, the set of curated bins in the `filtered_bins/` folder, or your own set of refined bins. Whichever set you choose, make sure you select the correct input folder and extension setting as it may differ from the example here.

When the task completes, you will have a number of output files provided. The main ones to look for are `gtdbtk.bac120.summary.tsv` and `gtdbtk.arch122.summary.tsv` which report the taoxnomies for your MAGs, split at the domain level. These file are only written if MAGs that fall into the domain were found in your data set, so for this exercise we do not expect to see the `gtdbtk.arch122.summary.tsv` file.

If you are interested in performing more detailed phylogenetic analysis of the data, the filtered multiple sequence alignment (MSA) for the data are provided in the `gtdbtk.bac120.msa.fasta` and `gtdbtk.arch122.msa.fasta` files.

Have a look at your resulting taxonomy. The classification of your MAGs will be informative when addressing your research goal for this workshop.

---

### Overview of using `vContact2` to predict taxonomy of viral contigs

Even more so than prokaryote taxonomy, establishing a coherent system for viral taxonomy is complex and continues to evolve. Just in the last year, the International Committee on Taxonomy of Viruses ([ICTV](https://talk.ictvonline.org/)) overhauled the classification code into [15 hierarchical ranks](https://www.nature.com/articles/s41564-020-0709-x). Furthermore, the knowledge gap in databases of known and taxonomically assigned viruses remains substantial, and so identifying the putative taxonomy of viral contigs from environmental metagenomics data remains challenging.

There are a number of approaches that can be used to attempt to predict the taxonomy of the set of putative viral contigs output by programs such as `VIBRANT`, `VirSorter`, and `VirFinder`. For example, [vContact2](https://www.nature.com/articles/s41587-019-0100-8) is one such method that uses 'guilt-by-contig-association' to predict the potential taxonomy of viral genomic sequence data based on relatedness to known viruses within a reference database (such as viral RefSeq). The principle is that, to the extent that the 'unknown' viral contigs cluster closely with known viral genomes, we can then expect that they are closely related enough to be able to predict a shared taxonomic rank. Anecdotally however, in my own experience with this processes I have been able to predict the taxonomy of *very* few of the viral contigs ouput by `VIBRANT`, `VirSorter`, or `VirFinder` from an environmental metagenomic data set (due to not clustering closely enough with known viruses in the reference database).

Running `vContact2` can require a reasonable amount of computational resources, and so we won't be running this in the workshop today. The required process is outlined below for reference, should you wish to experiment with this on your own data in the future.

Further information for installation and running of vContact2 can also be found [here](https://bitbucket.org/MAVERICLab/vcontact2/src/master/).

### Overview of viral taxonomy prediction via `vContact2`

NOTE: these steps are based on having `vContact2` set up as a `conda` environment. This documentation will be updated should `vContact2` become available as a NeSI module.

**1. Predict genes via `prodigal`**

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy
```

Example slurm script:

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J prodigal
#SBATCH --time 00:05:00
#SBATCH --mem=1GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -e prodigal.err
#SBATCH -o prodigal.out

# Load dependencies
module load prodigal/2.6.3-GCCcore-7.4.0

# Set up working directories
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy
mkdir -p viral_taxonomy

# Run main analyses 
srun prodigal -p meta -q \
-i ../7.viruses/checkv/checkv_combined.fna \
-a viral_taxonomy/checkv_combined.faa 
```

**2. Generate required mapping file for `vContact2`**

Use `vContact2`'s `vcontact2_gene2genome` script to generate the required mapping file from the output of `prodigal`.

*NOTE: update `/path/to/conda/envs/vContact2/bin` in the below script to the appropraite path.*

```bash
# activate vcontact2 conda environment
module purge
module load Miniconda3
source activate vContact2

# Load dependencies
export PATH="/path/to/conda/envs/vContact2/bin:$PATH"
module load DIAMOND/0.9.32-GCC-9.2.0
module load MCL/14.137-gimkl-2020a

# run vcontact2_gene2genome
vcontact2_gene2genome -p viral_taxonomy/checkv_combined.faa -o viral_taxonomy/viral_genomes_g2g.csv -s 'Prodigal-FAA'

# deactivate conda environment
conda deactivate
```

**3. Run `vContact2`**

Example slurm script:

*NOTE: update `/path/to/conda/envs/vContact2/bin` and `/path/to/conda/envs/vContact2/bin/cluster_one-1.0.jar` in the below script to the appropraite paths.*

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J vcontact2
#SBATCH --time 02:00:00
#SBATCH --mem=20GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -e vcontact2.err
#SBATCH -o vcontact2.out

# Set up working directories
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/viral_taxonomy/

# activate vcontact2 conda environment
module purge
module load Miniconda3
source activate vContact2

# Load dependencies
export PATH="/path/to/conda/envs/vContact2/bin:$PATH"
module load DIAMOND/0.9.32-GCC-9.2.0
module load MCL/14.137-gimkl-2020a

# Run vcontact2
srun vcontact2 \
-t 20 \
--raw-proteins checkv_combined.faa \
--rel-mode Diamond \
--proteins-fp viral_genomes_g2g.csv \
--db 'ProkaryoticViralRefSeq201-Merged' \
--c1-bin /path/to/conda/envs/vContact2/bin/cluster_one-1.0.jar \
--output-dir vConTACT2_Results

# deactivate conda environment
conda deactivate
```

**4. Predict taxonomy of viral contigs based on ouput of `vContact2`**

`vContact2` doesn't actually *assign* taxonomy to your input viral contigs. It instead provides an output outlining which reference viral genomes your viral contigs clustered with (if they clustered with any at all). Based on how closely they clustered with any reference genome(s), you can then use this to *predict* the likely taxonomy of the contig. 

Note from the `vContact2` online docs:

> One important note is that the taxonomic information is not included for user sequences. This means that each user will need to find their genome(s) of interest and check to see if reference genomes are located in the same VC. If the user genome is within the same VC subcluster as a reference genome, then there's a very high probability that the user genome is part of the same genus. If the user genome is in the same VC but not the same subcluster as a reference, then it's highly likely the two genomes are related at roughly genus-subfamily level. If there are no reference genomes in the same VC or VC subcluster, then it's likely that they are not related at the genus level at all.

The summary output of `vContact2` is the file `vConTACT2_Results/genome_by_genome_overview.csv`. As the comment above notes, one approach would be to search this file for particular contigs of interest, and see if any reference genomes fall into the same viral cluster (VC), using this reference to predict the taxonomy of the contig of interest.

The following `python` script is effectively an automated version of this for all input contigs (*Note: this script has not been widely tested, and should be used with some degree of caution*). This script groups contigs (and reference genomes) that fall into each VC together, and then for each included contig and genome, outputs a list of all taxonomies (at the ranks of 'Order', 'Family', and 'Genus' separately) that were found in that cluster. The predictions (i.e. the list of all taxonomies found in the same VC) for each is output to the table `tax_predict_table.txt`. 

*NOTE: The taxonomies are deliberately enclosed in square brackets (`[ ]`) to highlight the fact that these are *predictions*, rather than definitive taxonomy *assignments*.

For future reference, a copy of this script is available for download [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/scripts/vcontact2_tax_predict.py)

```bash
module load Python/3.8.2-gimkl-2020a

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/

./vcontact2_tax_predict.py \
-i viral_taxonomy/vConTACT2_Results/genome_by_genome_overview.csv \
-o viral_taxonomy/
```

---
