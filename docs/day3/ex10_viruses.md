# Identifying viral contigs in metagenomic data

!!! info "Objectives"

    * [Identifying viral contigs using `VIBRANT`](#identifying-viral-contigs-using-vibrant)
    * [Examine prophage identified by `VIBRANT`](#examine-outputs-of-vibrant)
    * [Examine viral metabolism and auxiliary metabolic genes (AMGs) outputs from `VIBRANT`](#examine-outputs-of-vibrant)
    * [Check quality and estimate completeness of the viral contigs via `CheckV`](#check-quality-and-estimate-completeness-of-the-viral-contigs-via-checkv)
    * [Introduction to *vContact2* for predicting taxonomy of viral contigs](#introduction-to-vcontact2-for-predicting-taxonomy-of-viral-contigs)
    
---

### Identifying viral contigs

Viral metagenomics is a rapidly progressing field, and new software are constantly being developed and released each year that aim to better identify and characterise viral genomic sequences from assembled metagenomic sequence reads. 

Currently, the most commonly used methods are [VIBRANT](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0), [VirSorter](https://peerj.com/articles/985/), and [VirFinder](https://link.springer.com/epdf/10.1186/s40168-017-0283-5?author_access_token=YQgkTWibFIFPtRICkTjZF2_BpE1tBhCbnbw3BuzI2RMCpVMGldKV8DA9scozc7Z-db3ufPFz9-pswHsYVHyEsCrziBuECllLPOgZ6ANHsMeKF5KejrdDKdeASyDkxB5wfFDq523QSd01cnqxCLqCiQ%3D%3D) (or the machine learning implementation of this, [DeepVirFinder](https://github.com/jessieren/DeepVirFinder)). A number of recent studies use either one of these tools, or a combination of several at once.

!!! info ""

    === "[VIBRANT](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0)"

        Uses a machine learning approach based on protein similarity (non-reference-based similarity searches with multiple HMM sets), and is in principle applicable to bacterial and archaeal DNA and RNA viruses, integrated proviruses (which are excised from contigs by `VIBRANT`), and eukaryotic viruses. 

        **More info** [github.com/AnantharamanLab/VIBRANT](https://github.com/AnantharamanLab/VIBRANT)

    === "[VirSorter](https://peerj.com/articles/985/)"

        Uses a predicted protein homology reference database-based approach, together with searching for a number of pre-defined metrics based on known viral genomic features. The authors note that `VirSorter` is currently best applied to DNA viruses (including prophage, which are also excised from contigs by `VirSorter`), but is likely poor with RNA viruses (from metatranscriptome data) and is also poor with eukaryotic viruses (as the database currently lacks eukaryotic viruses, and the genomic features incorporated were developed based on viruses of prokaryotes). 
        

        * *Hot off the press: an early test version of `VirSorter2` was recently released, which expands `VirSorter`'s target viruses to now include dsDNAphage, ssDNA and RNA viruses, and the viral groups Nucleocytoviricota and lavidaviridae.* 

        **More info** [github.com/simroux/VirSorter](https://github.com/simroux/VirSorter)

    === "[DeepVirFinder](https://github.com/jessieren/DeepVirFinder)"

        Uses a machine learning based approach based on *k*-mer frequencies. Having developed a database of the differences in *k*-mer frequencies between prokaryote and viral genomes, `VirFinder` examines assembled contigs and identifies whether their *k*-mer frequencies are comparable to known viruses in the database, using this to predict viral genomic sequence. This method has some limitation based on the viruses that were included when building the database (bacterial DNA viruses, but very few archaeal viruses, and, at least in some versions of the software, no eukaryotic viruses). However, tools are also provided to build your own database should you wish to develop an expanded one. Due to its distinctive *k*-mer frequency-based approach, `VirFinder` may also have the capability of identifying some novel viruses overlooked by tools such as `VIBRANT` or `VirSorter`.

        **More info** [github.com/jessieren/DeepVirFinder](https://github.com/jessieren/DeepVirFinder)


---

### Identifying viral contigs using `VIBRANT`

For this exercise, we will use `VIBRANT` to identify viral contigs from our assembled contigs. As an added bonus, `VIBRANT` also provides outputs for gene and protein predictions, identified integrated prophage, and metabolism (including viral metabolism and auxiliary-metabolic genes (AMGs) - prokaryote host-derived genes that have been integrated in, and are encoded by, viral genomes). 

These exercises will take place in the `7.viruses/` folder.

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/7.viruses/
```

#### Run *Vibrant*

For this we will input the assembled contigs from the `SPAdes` assembly we performed earlier. These assembly files are available at `7.viruses/spades_assembly/`

*NOTE: The path to the required databases for the `VIBRANT` NeSI module have been linked to the variable `$DB_PATH`. This needs to be included in the command (using the `-d` flag) for `VIBRANT` to be able to locate them.*

Create a new script

```bash
nano vibrant.sl
```
!!! warning "Warning"

    Paste in the following (updating `<YOUR FOLDER>`)

!!! terminal "code"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      vibrant
    #SBATCH --time          00:40:00
    #SBATCH --mem           4GB
    #SBATCH --cpus-per-task 16
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out

    # Load modules
    module purge
    module load VIBRANT/1.2.1-gimkl-2020a

    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/7.viruses/

    # Output directory
    mkdir -p vibrant

    # Run VIBRANT
    VIBRANT_run.py -t $SLURM_CPUS_PER_TASK \
                   -i spades_assembly/spades_assembly.m1000.fna \
                   -d $DB_PATH \
                   -folder vibrant/
    ```

Submit the script as a slurm job

```bash
sbatch vibrant.sl
```

#### Examine outputs of *VIBRANT*

Exercise: `VIBRANT` provides a number of different outputs. Explore through the various folders within the `vibrant/` folder and identify some that might be of particular interest. Open some of these files to see if you can find the following information:

!!! quote ""

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

Create a new script

```bash
nano checkv.sl
```
!!! warning "Warning"

    Paste in the following (updating `<YOUR FOLDER>`)

!!! terminal "code"

    ```bash
    #!/bin/bash -e

    #SBATCH --account       nesi02659
    #SBATCH --job-name      CheckV
    #SBATCH --time          00:10:00
    #SBATCH --mem           3GB
    #SBATCH --cpus-per-task 10
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out

    # Load modules
    module purge
    module load CheckV/1.0.1-gimkl-2022a-Python-3.10.5

    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/7.viruses/

    # Output directory
    mkdir -p checkv_out

    # Run CheckV
    checkv_in="vibrant/VIBRANT_spades_assembly.m1000/VIBRANT_phages_spades_assembly.m1000/spades_assembly.m1000.phages_combined.fna"

    checkv end_to_end $checkv_in checkv_out/ -t $SLURM_CPUS_PER_TASK
    ```

Submit the script as a slurm job

```bash
sbatch checkv.sl
```

#### Examine outputs of *CheckV*

`CheckV` provides summary outputs for contamination, completeness, repeats, and an overall quality summary. Have a brief look at some examples of the information you can draw from each of these `CheckV` outputs. 

!!! question "Exercise: Examining `checkv_out/quality_summary.tsv`"

    * How many viral contigs meet the "High-quality" (MIUViG) standard?
    * How many might we consider "complete" genomes based on `CheckV`'s completeness estimation?
    * Were any terminal repeat regions identified for any of the contigs?
    * How many prophage were identified by `CheckV`? Why might this differ to `VIBRANT`'s identification of prophage above?
    * Are there any suspicious contigs that you might want to flag for closer examination (and/or careful consideration in downstream analyses)?

---

### Introduction to *vContact2* for predicting taxonomy of viral contigs

Even more so than prokaryote taxonomy, establishing a coherent system for viral taxonomy is complex and continues to evolve. In 2020, the International Committee on Taxonomy of Viruses ([ICTV](https://talk.ictvonline.org/)) overhauled the classification code into [15 hierarchical ranks](https://www.nature.com/articles/s41564-020-0709-x). Furthermore, the knowledge gap in databases of known and taxonomically assigned viruses remains substantial, and so identifying the putative taxonomy of viral contigs from environmental metagenomics data remains challenging.

There are a number of approaches that can be used to attempt to predict the taxonomy of the set of putative viral contigs output by programs such as `VIBRANT`, `VirSorter`, and `VirFinder`. [vContact2](https://www.nature.com/articles/s41587-019-0100-8) is one such method that uses 'guilt-by-contig-association' to predict the potential taxonomy of viral genomic sequence data based on relatedness to known viruses within a reference database (such as viral RefSeq). The principle is that, to the extent that the 'unknown' viral contigs cluster closely with known viral genomes, we can then expect that they are closely related enough to be able to predict a shared taxonomic rank. 

!!! note "Note"
    Anecdotally, however, in my own experience with this processes I have unfortunately been unable to predict the taxonomy of the vast majority of the viral contigs ouput by `VIBRANT`, `VirSorter`, or `VirFinder` from an environmental metagenomic data set (due to not clustering closely enough with known viruses in the reference database).

Running `vContact2` can require a considerable amount of computational resources, and so we won't be running this in the workshop today. The required process is outlined for reference in an [Appendix for this exercise](../resources/4_APPENDIX_ex11_viral_taxonomy_prediction_via_vContact2.md), should you wish to experiment with this on your own data in the future. 

For today, we have provided the final two output files from this process when applied to our mock metagenome data. These can be viewed in the folder `8.coverage_and_taxonomy/vConTACT2_Results/` via `head` or `less`.

```bash
less vConTACT2_Results/genome_by_genome_overview.csv
```

```bash
less vConTACT2_Results/tax_predict_table.txt
```

A few notes to consider: 

!!! quote ""

    * You will see that the `genome_by_genome_overview.csv` file contains entries for the full reference database used as well as the input viral contigs (contigs starting with `NODE`). 
    * You can use a command such as `grep "NODE" vConTACT2_Results/genome_by_genome_overview.csv | less` to view only the lines for the input contigs of interest. 
      * Note also that these lines however will *not* contain taxonomy information. 
      * See the notes in the [Appendix](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/APPENDIX_ex11_viral_taxonomy_prediction_via_vContact2.md) for further information about why this might be.
    * As per the notes in the [Appendix](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/APPENDIX_ex11_viral_taxonomy_prediction_via_vContact2.md), the `tax_predict_table.txt` file contains *predictions* of potential taxonomy (and or taxonom*ies*) of the input viral contigs for order, family, and genus, based on whether they clustered with any viruses in the reference database.
      * Bear in mind that these may be lists of *multiple* potential taxonomies, in the cases where viral contigs clustered with multiple reference viruses representing more than one taxonomy at the given rank.
      * *NOTE: The taxonomies are deliberately enclosed in square brackets (`[ ]`) to highlight the fact that these are **predictions**, rather than definitive taxonomy **assignments**.*
    
---
