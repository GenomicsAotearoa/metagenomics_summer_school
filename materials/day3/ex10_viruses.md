# Identifying viral contigs in metagenomics data

### Objectives

* Identifying viral contigs
* Examine prophage identified by `VIBRANT`
* Examine viral metabolism and auxiliary metabolic genes (AMGs) outputs from `VIBRANT`
* Check quality and estimate completeness of the viral contigs via `CheckV`

---

### Identifying viral contigs

Viral metagenomics is a rapidly progressing field, and new software are constantly being developed and released each year that aim to better identify and characterise viral genomic sequences from assembled metagenomic sequence reads. 

Currently, the most commonly used methods are [VIBRANT](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0), [VirSorter](https://peerj.com/articles/985/), and [VirFinder](https://link.springer.com/epdf/10.1186/s40168-017-0283-5?author_access_token=YQgkTWibFIFPtRICkTjZF2_BpE1tBhCbnbw3BuzI2RMCpVMGldKV8DA9scozc7Z-db3ufPFz9-pswHsYVHyEsCrziBuECllLPOgZ6ANHsMeKF5KejrdDKdeASyDkxB5wfFDq523QSd01cnqxCLqCiQ%3D%3D) (or the machine learning implementation of this, [DeepVirFinder](https://github.com/jessieren/DeepVirFinder)). A number of recent studies use either one of these tools, or a combination of several at once.

[VIBRANT](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0) uses a machine learning approach based on protein similarity (non-reference-based similarity searches with multiple HMM sets), and is in principle applicable to bacterial and archaeal DNA and RNA viruses, integrated proviruses (which are excised from contigs by `VIBRANT`), and eukaryotic viruses. 

[VirSorter](https://peerj.com/articles/985/) uses a predicted protein homology reference database-based approach, together with searching for a number of pre-defined metrics based on known viral genomic features. The authors note that `VirSorter` is currently best applied to DNA viruses (including prophage, which are also excised from contigs by `VirSorter`), but is likely poor with RNA viruses (from metatranscriptome data) and is also poor with eukaryotic viruses (as the database currently lacks eukaryotic viruses, and the genomic features incorporated were developed based on viruses of prokaryotes). 

* *Hot off the press: an early test version of `VirSorter2` was recently released, which expands `VirSorter`'s target viruses to now include dsDNAphage, ssDNA and RNA viruses, and the viral groups Nucleocytoviricota and lavidaviridae.* 

[DeepVirFinder](https://github.com/jessieren/DeepVirFinder) uses a machine learning based approach based on *k*-mer frequencies. Having developed a database of the differences in *k*-mer frequencies between prokaryote and viral genomes, `VirFinder` examines assembled contigs and identifies whether their *k*-mer frequencies are comparable to known viruses in the database, using this to predict viral genomic sequence. This method has some limitation based on the viruses that were included when building the database (bacterial DNA viruses, but very few archaeal viruses, and, at least in some versions of the software, no eukaryotic viruses). However, tools are also provided to build your own database should you wish to develop an expanded one. Due to its distinctive *k*-mer frequency-based approach, `VirFinder` may also have the capability of identifying some novel viruses overlooked by tools such as `VIBRANT` or `VirSorter`.

**Further documentation**

* [VIBRANT](https://github.com/AnantharamanLab/VIBRANT)
* [VirSorter](https://github.com/simroux/VirSorter)
* [DeepVirFinder](https://github.com/jessieren/DeepVirFinder)

---

### Identifying viral contigs using `VIBRANT`

For this exercise, we will use `VIBRANT` to identify viral contigs from our assembled contigs. As an added bonus, `VIBRANT` also provides outputs for gene and protein predictions, identified integrated prophage, and metabolism (including viral metabolism and auxiliary-metabolic genes (AMGs) - prokaryote host-derived genes that have been integrated in, and are encoded by, viral genomes). 

These exercises will take place in the `7.viruses/` folder.

#### Run *Vibrant*

For this we will input the assembled contigs from the `SPAdes` assembly we performed earlier. These assembly files are available at `7.viruses/spades_assembly/`

*NOTE: The path to the required databases for the `VIBRANT` NeSI module have been linked to the variable `$DB_PATH`. This needs to be included in the command (using the `-d` flag) for `VIBRANT` to be able to locate them.*

Example slurm script:

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH --res SummerSchool
#SBATCH -J vibrant
#SBATCH --time 00:30:00
#SBATCH --mem=4GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -e vibrant.err
#SBATCH -o vibrant.out
#SBATCH --profile=task

# Load VIBRANT module
module load VIBRANT/1.2.1-gimkl-2020a

# Set up working directories
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/7.viruses/
mkdir -p vibrant

# Run VIBRANT
srun VIBRANT_run.py -t 16 \
-i spades_assembly/spades_assembly.m1000.fna \
-d $DB_PATH \
-folder vibrant/
```

#### Examine outputs of *VIBRANT*

Exercise: `VIBRANT` provides a number of different outputs. Explore through the various folders within the `vibrant/` folder and identify some that might be of particular interest. Open some of these files to see if you can find the following information:

* How many viral contigs did `VIBRANT` identify?
* Of these, how many did `VIBRANT` identify as prophage?
* What are some annotations of interest within the output annotations file? 
  * *NOTE: the `VIBRANT` annotations file includes multiple columns for both **prokaryote** and **viral** protein predictions. Be careful as to which column you are looking at (as well as its associated confidence score) when assessing viral annotations vs. AMGs*.
* Among these annotations, how many were flagged as AMGs?
* What broad metabolic categories did the AMGs fall into? 
  * *NOTE: as well as providing this information in table format, `VIBRANT` also generates a number of summary figures. These figures are not viewable within NeSI, but can be downloaded and opened on your local machine (e.g. via `scp ...`) if you wish to look more closely at these*.
* Discussion point: How might we investigate whether identified putative AMGs are actually *within* the viral genomes, rather than residual contaminating host genomic sequence attached to the end of integrated prophage (but incompletely trimmed off in the excision process)?

---

### Check quality and estimate completeness of the viral contigs via *CheckV*

[CheckV](https://www.biorxiv.org/content/10.1101/2020.05.06.081778v1.abstract) was recently developed as an analogue to `CheckM`. `CheckV` first performs a 'contaminating sequence' trim, removing any retained (prokaryote) host sequence on the end of contigs with integrated prophage, and then assesses the quality and completeness of the assembled viral contigs. The quality of the contigs are also categoriesed based on the recently developed [Minimum Information about an Unclutivated Virus Genome](https://www.nature.com/articles/nbt.4306) (MIUViG) standards for reporting sequences of unclutivated virus geneomes (such as those recovered from metagenomic sequencing data). The MIUViG were developed as an extension of the [Minimum Information about any (x) Sequence](https://www.nature.com/articles/nbt.1823) ([MIxS](https://gensc.org/mixs/)) standards, which include, among others, standards for Metagenome-Assembled Genomes (MIMAG).

Installation and further instructions for `CheckV` can be found [here](https://bitbucket.org/berkeleylab/checkv/src/master/).

#### Run *CheckV*

Run `CheckV` providing the *fastA* file of combined (virus and prophage) viral contigs output by `VIBRANT` as input (`spades_assembly.m1000.phages_combined.fna`).

Example slurm script:

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH --res SummerSchool
#SBATCH -J checkv
#SBATCH --time 00:10:00
#SBATCH --mem=3GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -e checkv.err
#SBATCH -o checkv.out
#SBATCH --profile=task

# Load the module
module load CheckV/0.7.0-gimkl-2020a-Python-3.8.2

# Set up working directories
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/7.viruses/
mkdir -p checkv_out

# Run main analyses 
checkv_in="vibrant/VIBRANT_spades_assembly.m1000/VIBRANT_phages_spades_assembly.m1000/spades_assembly.m1000.phages_combined.fna"

srun checkv end_to_end ${checkv_in} checkv_out -t 10 --quiet
```

#### Examine outputs of *CheckV*

`CheckV` provides summary outputs for contamination, completeness, repeats, and an overall quality summary. Have a brief look at some examples of the information you can draw from each of these `CheckV` outputs. 

Exercise: Examining `checkv_out/quality_summary.tsv`

* How many viral contigs meet the "High-quality" (MIUViG) standard?
* How many might we consider "complete" genomes based on `CheckV`'s completeness estimation?
* Were any terminal repeat regions identified for any of the contigs?
* How many prophage were identified by `CheckV`? Why might this differ to `VIBRANT`'s identification of prophage above?
* Are there any suspicious contigs that you might want to flag for closer examination (and/or careful consideration in downstream analyses)?

---
