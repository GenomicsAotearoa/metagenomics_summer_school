### Gene prediction and annotation with *DRAM* (Distilled and Refined Annotation of Metabolism) 

[DRAM](http://dx.doi.org/10.1093/nar/gkaa621) is a tool designed to profile microbial (meta)genomes for metabolisms known to impact ecosystem functions across biomes. `DRAM` annotates MAGs and viral contigs using KEGG (if provided by user), UniRef90, PFAM, CAZy, dbCAN, RefSeq viral, VOGDB (Virus Orthologous Groups) and the MEROPS peptidase database. It is also highly customizable to other custom user databases. 

`DRAM` only uses assembly-derived *fastA* files input by the user. These input files may come from unbinned data (metagenome contig or scaffold files) or genome-resolved data form one or many organisms (isolate genomes, single-amplified genome (SAGs), MAGs).

`DRAM` is run in two stages: annotation and distillation. 

![](https://github.com/mcastudillo/MAG-annotation-with-DRAM/blob/main/figures/DRAM_workflow.png)

#### Annotation

The first step in `DRAM` is to annotate genes by assigning database identifiers to genes. Short contigs (default < 2,500 bp) are initially removed. Then, `Prodigal` is used to detect open reading frames (ORFs) and to predict their amino acid sequences. Next, `DRAM` searches all amino acid sequences against multiple databases, providing a single *Raw* output. When gene annotation is complete, all results are merged in a single tab-delimited annotation table, including best hit for each database for user comparison. 

#### Distillation 

After genome annotation, a distill step follows with the aim to curate these annotations into useful functional categories, creating genome statistics and metabolism summary files, and stored in the *Distillate* output. The genome statistics provides most genome quality information required for [MIMAG](https://www.nature.com/articles/nbt.3893), including `GTDB-tk` and `checkM` information if provided by the user. Summarised metabolism table include the number of genes with specific metabolic function identifiers (KO, CAZY ID, etc) fore each genome, with information obtained from multiple databases. The *Distillate* output is then further distilled into the *Product*, an html file displaying a heatmap, as well as the corresponding data table. We will investigate all these files later on.  

---

### Annotation of the MAGs with *DRAM*

Beyond annotation, `DRAM` aims to be a data compiler. For that reason, output files from both `CheckM` and `GTDB_tk` steps can be input to `DRAM` to provide both taxonomy and genome quality information of the MAGs. 

For these exercises, we have copied the relevant input files into the folder `10.gene_annotation/DRAM_input_files/`. `gtdbtk.bac120.classification_pplacer.tsv` was taken from the earlier `8.coverage_and_taxonomy/gtdbtk_out/` outputs, and `checkm.txt` from the result of re-running `CheckM` on the final refined filtered bins in `6.bin_refinement/filtered_bins`.

The `CheckM` output file (`checkm.txt`) can be input as it is. However, in order to use the file with the `gtdb_tk` taxonomy (`gtdbtk.bac120.classification_pplacer.tsv`) we should modify it first to include column headers 'bin_id' and 'classification'

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


sed -i '1s/^/bin_id\tclassification\n/' DRAM_input_files/gtdbtk.bac120.classification_pplacer2.tsv


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

In default annotation mode, `DRAM` takes as only input the directory containing all the bins we would like to annotate in *fastA* format (either .fa or .fna). There are few parameters that can be modified if not using the default mode. Once the annotation step is done, the mode `distill` is used to summarise the obtained results. 

*NOTE: due to the increased memory requirements, UniRef90 database is not default and the flag `–use_uniref` should be specified in order to search amino acid sequences against UniRef90. In this exercise, due to memory and time constraints, we won't be using UniRef90 database.*

```bash

module purge
module load Miniconda3/4.8.3
module load gimkl/2020a

export CONDA_PKGS_DIRS=/nesi/project/nesi02659/.conda/pkgs
export CONDA_ENVS_PATH=/nesi/project/nesi02659/.conda/envs

source activate DRAM

DRAM.py --help

# usage: DRAM.py [-h] {annotate,annotate_genes,distill,strainer,neighborhoods} ...

# positional arguments:
#   {annotate,annotate_genes,distill,strainer,neighborhoods}
#    annotate            Annotate genomes/contigs/bins/MAGs
#    annotate_genes      Annotate already called genes, limited functionality compared to annotate
#    distill             Summarize metabolic content of annotated genomes
#    strainer            Strain annotations down to genes of interest
#    neighborhoods       Find neighborhoods around genes of interest

#optional arguments:
#  -h, --help            show this help message and exit

```

To run this exercise we first need to set up a slurm job. We will use the results for tomorrow's distillation step. 

*NOTE: Currently DRAM has to be run from the directory where* ```DRAM-setup.py``` *was ran in order to work, that is why we start the slurm script with* ```cd /nesi/project/nesi02659/.conda/dramdbsetup```


```bash
#!/bin/bash -e
#SBATCH -J DRAM_annotation
#SBATCH -A nesi02659
#SBATCH --res SummerSchool
#SBATCH --time=5:00:00
#SBATCH --mem=20Gb
#SBATCH -e slurm-DRAM_annot.%A-%a.err 
#SBATCH -o slurm-DRAM_annot.%A-%a.out 

cd /nesi/project/nesi02659/.conda/dramdbsetup

module purge
module load Miniconda3/4.8.3
module load gimkl/2020a

export CONDA_PKGS_DIRS=/nesi/project/nesi02659/.conda/pkgs
export CONDA_ENVS_PATH=/nesi/project/nesi02659/.conda/envs

source activate DRAM

DRAM.py annotate -i '/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/bins_for_DRAM/*.fna' \
--checkm_quality /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/DRAM_input_files/checkm.txt \
--gtdb_taxonomy /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/DRAM_input_files/gtdbtk.bac120.classification_pplacer.tsv \
-o /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/annotation_mgss

```

The program will take 4-4.5 hours to run, so we will submit the jobs and inspect the results tomorrow morning. 
