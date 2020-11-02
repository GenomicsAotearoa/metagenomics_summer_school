# Gene annotation

### Objectives

* **BLAST**-like gene annotations and domain annotations
* Overview of annotation databases
* Evaluating the quality of gene assignment
* *Discussion: Differences in taxonomies (GTDB, NCBI etc)*
* Annotating MAGs with **DRAM**

---

### *BLAST*-like gene annotations and domain annotations

Broadly speaking, there are two ways we perform gene annotations with protein sequences. Both compare our sequences of interest against a curated set of protein sequences for which function is known, or is strongly suspected. In each case, there are particular strenths to the approach and for particular research questions, one option may be favoured over another.

#### BLAST-like annotation

The first of these is the `BLAST` algorithm for sequence alignment. This approach performs pairwise alignment between the gene of interest (query sequence) and the sequences in the database (target sequence). `BLAST` searches each potential target sequence for *k*-mers identified in the query sequence. Where these *k*-mers are found in targets, the ends are extended out to try to create longer regions of highly similar sequence spans. Across this span, the tool identifies the longest span of characters (nucleotide or amino acid) that match within a scoring framework to return the length of the region (coverage) and the sequence identity over the span (identity).

The original tool for performing this kind of analysis was the `BLAST` tool. While `BLAST` and its variants are still excellent tools for performing this kind of sequence annotation, they suffer from a slow runtime speed due to the need to test each query sequence against every target sequence in the database. For this reason, several tools have been published which take the basic approach of `BLAST`, but augment it with methods to reduce the number of pairwise comparisons needed to identify targets with high sequence similarity to the query. Two popular pieces of software are the tools `usearch` [here](http://www.drive5.com/usearch/) and `diamond` [here](https://github.com/bbuchfink/diamond).

#### HMM-profiling of domains

An alternate method for attributing function to query sequences it to consider them as a collection of independently functioning protein folding domains. This is the approach used in the `HMMer` (http://hmmer.org/) software, and the *Pfam*, *TIGRfam*, and *PANTHER* databases. In these analyses, the database consists not of individual sequences, but of Hidden Markov models built from a collection of proteins that share a common domain. These profiles build out a statistical map of the amino acid transitions (from position to position), variations (differences at a position), and insertions/deletions between positions in the domain across the different observations in the training database and apply these maps to the query data.

---

### Annotating MAGs with against the *UniProt* database with *diamond*

For this exercise we are going to use diamond for performing our annotation. We have chosen to use this tool because it is faster than BLAST, and usearch comes with licencing restrictions that make it hard to work with in a shared computing environment like NeSI.

For this exercise we have created a diamond-compatible database from the 2018 release of the UniProt database

In general, diamond takes a simple pair of input files - the protein coding sequences we wish to annotate and the database we will use for this purpose. There are a few parameters that need to be tweaked for obtaining a useful output file, however.

```bash
module load DIAMOND/0.9.25-gimkl-2018b

diamond help
# diamond v0.9.25.126 | by Benjamin Buchfink <buchfink@gmail.com>
# Licensed under the GNU GPL <https://www.gnu.org/licenses/gpl.txt>
# Check http://github.com/bbuchfink/diamond for updates.

# Syntax: diamond COMMAND [OPTIONS]

# Commands:
# ...
# blastp  Align amino acid query sequences against a protein reference database
# ...

# General options:
# --threads (-p)         number of CPU threads
# --db (-d)              database file
# --out (-o)             output file
# --outfmt (-f)          output format
# ...
```

There are two output formats we can chose from which are useful for our analysis. We will obtain our output in the BLAST tabular format, which provides the annotation information in a simple-to-parse text file that can be viewed in any text or spreadsheet viewing tool. This will allow us to investigate and evaluate the quality of our annotations.  For now, just annotate a single bin:

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.gene_annotation/

diamond blastp -p 1 --db /nesi/project/nesi02659/mg_workshop/uniprot_nr_200213.diamond \
               --max-target-seqs 5 --evalue 0.001 \
               -q example_data/bin_0.filtered.genes.no_metadata.faa \
               --outfmt 6 -o bin_0.diamond.txt
```

Awkwardly, `diamond` does not provide the headers for what the columns in the output table mean. [This table](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) is a handy reference for how to interpret the output.

From here we can view important stastics for each query/target pairing such as the number of identify residues between sequences and the aligned length between query and target.

Before we proceed with this exercise, lets set up a slurm job to annotate each of our MAGs. We will use this tomorrow.

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J annotate_uniprot
#SBATCH --partition ga_bigmem
#SBATCH --res SummerSchool
#SBATCH --time 02:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e annotate_uniprot_dmnd.err
#SBATCH -o annotate_uniprot_dmnd.out

module load DIAMOND/0.9.25-gimkl-2018b

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.gene_annotation/

mkdir -p gene_annotations/

for prot_file in example_data/*.genes.no_metadata.faa;
do
  out_file=$(basename ${prot_file} .faa)

  diamond blastp -p 20 --max-target-seqs 5 --evalue 0.001 \
                 --db /nesi/nobackup/nesi02659/MGSS_resources/uniprot_nr_200213.diamond \ ## MODIFY PATHS ONCE THE DIRECTORY IS SET UP
                 -q ${prot_file} --outfmt 6 -o gene_annotations/${out_file}.uniprot.txt
done
```
---

## Annotating MAGs against the Pfam database with hmm

The standard software for performing this kind of annotation is [hmmer](http://hmmer.org/). Compared to BLAST, FASTA, and other sequence alignment and database search tools based on older scoring methodology, HMMER aims to be significantly more accurate and more able to detect remote homologs because of the strength of its underlying mathematical models. In the past, this strength came at significant computational expense, but in the new HMMER3 project, HMMER is now essentially as fast as BLAST. First, let's have a look at hmmer options.


```
module load HMMER/3.1b2-gimkl-2017a

hmmsearch -h 

# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: hmmsearch [options] <hmmfile> <seqdb>

Basic options:
  -h : show brief help on version and usage

Options directing output:
  -o <f>           : direct output to file <f>, not stdout
  -A <f>           : save multiple alignment of all hits to file <f>
  --tblout <f>     : save parseable table of per-sequence hits to file <f>
  --domtblout <f>  : save parseable table of per-domain hits to file <f>
  --pfamtblout <f> : save table of hits and domains to file, in Pfam format <f>
  --acc            : prefer accessions over names in output
  --noali          : don't output alignments, so output is smaller
  --notextw        : unlimit ASCII text output line width
  --textw <n>      : set max width of ASCII text output lines  [120]  (n>=120)

Options controlling reporting thresholds:
  -E <x>     : report sequences <= this E-value threshold in output  [10.0]  (x>0)
  -T <x>     : report sequences >= this score threshold in output
  --domE <x> : report domains <= this E-value threshold in output  [10.0]  (x>0)
  --domT <x> : report domains >= this score cutoff in output

Options controlling inclusion (significance) thresholds:
  --incE <x>    : consider sequences <= this E-value threshold as significant
  --incT <x>    : consider sequences >= this score threshold as significant
  --incdomE <x> : consider domains <= this E-value threshold as significant
  --incdomT <x> : consider domains >= this score threshold as significant

Options controlling model-specific thresholding:
  --cut_ga : use profile's GA gathering cutoffs to set all thresholding
  --cut_nc : use profile's NC noise cutoffs to set all thresholding
  --cut_tc : use profile's TC trusted cutoffs to set all thresholding

Options controlling acceleration heuristics:
  --max    : Turn all heuristic filters off (less speed, more power)
  --F1 <x> : Stage 1 (MSV) threshold: promote hits w/ P <= F1  [0.02]
  --F2 <x> : Stage 2 (Vit) threshold: promote hits w/ P <= F2  [1e-3]
  --F3 <x> : Stage 3 (Fwd) threshold: promote hits w/ P <= F3  [1e-5]
  --nobias : turn off composition bias filter

Other expert options:
  --nonull2     : turn off biased composition score corrections
  -Z <x>        : set # of comparisons done, for E-value calculation
  --domZ <x>    : set # of significant seqs, for domain E-value calculation
  --seed <n>    : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
  --tformat <s> : assert target <seqfile> is in format <s>: no autodetection
  --cpu <n>     : number of parallel CPU workers to use for multithreads

```

We are now going to submit another slurm job to annotate our MAGs using the [Pfam database](https://pfam.xfam.org/). Matching sequences to a Pfam entry allows us to transfer the functional information from an experimentally characterised sequence to uncharacterised sequences in the same entry. Pfam then provides comprehensive annotation for each entry.

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J annotate_pfam
#SBATCH --partition ga_bigmem
#SBATCH --res SummerSchool
#SBATCH --time 02:00:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e annotate_pfam_hmm.err
#SBATCH -o annotate_pfam_hmm.out

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.gene_annotation/


for prot_file in example_data/*.genes.no_metadata.faa;
do
  out_file=$(basename ${prot_file} .faa)

  hmmsearch --tblout ${out_file}.pfam.txt -E 1e-3 --cpu 10 ../databases/Pfam-A.hmm ${prot_file}
  
done

```
---

### Evaluating the quality of gene assignment

Determining how trustworthy a gene annotation is can be a very tricky process. How similar to do protein sequences need to be to perform the same function? The answer is surprisingly low - a bioinformatic analysis performed in 1999 identified that proteins with as little as 20 - 35% sequence identity can still share the same function ([Rost, 1999](https://doi.org/10.1093/protein/12.2.85)) but this is not a universal occurrence. When evaluating your annotations, we consider the following questions:

1. What is the amino acid identity along the aligned region?
1. What is the amino acid *similarity* between the aligned region?
1. What is the coverage as a percentage or the query and target genes?
1. If we infer a phylogeny of this query gene with references from the target family, is a stable tree resolved?
1. Does the inclusion of this gene function make sense in the context of the organism's taxonomy?
1. Does the gene sit on a long contig that is core to the MAG, or is it a short contig carrying only a single gene?
1. If we are uncertain of a particular annotation, does the predicted gene occur in an operon? If so, are the other genes present in the annotation?

We must also remain aware of the potential for incorrectly annotated genes in the annotation database and that proteins can perform multiple functions (and may therefore be attributed multiple, inconsistent annotations). Furthermore, it is also important to consider exactly which part of the target gene the alignment is happening across. There are several catalytic centers of enzymes, such as the Fe-S cofactor, which are shared across many different proteins and if your annotation is only spanning one of these regions then it may simply be the case that you are identifying a generic electron accepting or donating domain.

---

### Differences in taxonomies

Another way to determine if a annotation 'belongs' in the MAG of interest is to consider the predicted taxonomy of the query gene with that of the MAG itself. For example, if you detect a *Desulfovibrio*-like *dsrA* seuqence in a bin that has been classified as belonging to the genus *Desulfovibrio* then it is probably a safe bet that the annotation is correct.

However, when comparing taxonomic assignments, it is important to be aware of the differing taxonomic schemes that are circulating in the microbiological and bioinformatic literature and to know how to reconcile their differences. Similar to how the 16S rRNA gene taxonomies provided by *SILVA*, *Greengenes*, and *RDP* taxonomies all differ in some aspects there are multiple competing taxonomies in protein databases.

#### Examples of various genome-level and protein taxonomies

* [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy)
* [Genome Taxonomy Database](https://gtdb.ecogenomic.org/)

This problem exists because despite the existance of a formal [Code](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.000778) for the naming of bacteria and archaea, because 

1. There are no rules governing how we define the grouping of these names together, other than for type species
1. Defunct synonyms and basonyms are not correctly purged from taxonomy lists (this is quite noticable with the NCBI taxonomy)
1. Valid names cannot be assigned for uncultivate organisms, meaning there are many informal placeholder names in the literature. For example, clades like WPS-2, SAR324, and SAUL are widely cited in the literature despite having no official standing

It is therefore important to periodically sanity check your taxonomic annotations in order to avoid splitting taxa based on spelling differences or the use of historic names that have since been reclassified.

---

## Gene prediction and annotation with DRAM (Distilled and Refined Annotation of Metabolism) 

DRAM is a tool designed to profile microbial (meta)genomes for metabolisms known to impact ecosystem functions across biomes. DRAM annotates MAGs and viral contigs using KEGG (if provided by user), UniRef90, PFAM, CAZy, dbCAN, RefSeq viral, VOGDB (Virus Orthologous Groups) and the MEROPS peptidase database. It is also highly customizable to other custom user databases. 
DRAM only uses assembly-derived FASTA files input by the user. These input files may come from unbinned data (metagenome contig or scaffols files) or genome-resolved data form one or many organisms (isolate genomes, single-amplified genome (SAGs), MAGs).
DRAM is run in two stages: annotation and distillation. 

![](https://github.com/mcastudillo/MAG-annotation-with-DRAM/blob/main/figures/DRAM_workflow.png)

### Annotation
The first step in DRAM is to annotate genes by assigning database identifiers to genes. Short contigs (default < 2,500 bp) are initially removed. Then, Prodigal is used to detect open reading frames (ORFs) and to predict their amino acid sequences. Next, DRAM searches all amino acid sequences against multiple databases, providing a single *Raw* output. When gene annotation is complete, all results are merged in a single tab-delimited annotation table, including best hit for each database for user comparison. 


### Distillation 
After genome annotation, a distill step follows with the aim to curate these annotations into useful functional categories, creating genome statistics and metabolism summary files, and stored in the *Distillate* output. The genome statistics provides most genome quality information required for [MIMAG](https://www.nature.com/articles/nbt.3893), including GTDB-tk and checkM information if provided by user. Summarised metabolism table include the number of genes with specific metabolic function identifiers (KO, CAZY ID, etc) fore each genome, with information obtained from multiple databases. The *Distillate* output is then further distilled into the *Product*, an html file displaying a heatmap, as well as the corresponding data table. We will investigate all these files later on.  



## Annotation of the MAGs with DRAM 
Beyond annotation, DRAM aims to be a data compiler. For that reason, output files from both CheckM and GTDB_tk steps can be input to DRAM to provide both taxonomy and genome quality information of the MAGs. CheckM output file (`checkm.txt`) can be input as it is. However, in order to use the file with the gtdb_tk taxonomy (`gtdbtk.bac120.classification_pplacer.tsv`) we should modify it first to include column headers 'bin_id' and 'classification'

```
nano gtdbtk.bac120.classification_pplacer.tsv

# bin_1   d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__
# bin_0   d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Synechococcales;f__Cyanobiaceae;g__Prochlorococcus_C;s__
# bin_9   d__Bacteria;p__Planctomycetota;c__Brocadiae;o__Brocadiales;f__Brocadiaceae;g__;s__
# bin_2   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__
# bin_3   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__
# bin_4   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Nitrosomonadaceae;g__Nitrosomonas;s__
# bin_5   d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Nitrobacter;s__
# bin_7   d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Arcobacter;s__
# bin_6   d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Nautiliales;f__Nautiliaceae;g__;s__
# bin_8   d__Bacteria;p__Desulfobacterota_A;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio;s__


# bin_id  classification
# bin_1   d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__
# bin_0   d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Synechococcales;f__Cyanobiaceae;g__Prochlorococcus_C;s__
# bin_9   d__Bacteria;p__Planctomycetota;c__Brocadiae;o__Brocadiales;f__Brocadiaceae;g__;s__
# bin_2   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__
# bin_3   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__
# bin_4   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Nitrosomonadaceae;g__Nitrosomonas;s__
# bin_5   d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Nitrobacter;s__
# bin_7   d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Arcobacter;s__
# bin_6   d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Nautiliales;f__Nautiliaceae;g__;s__
# bin_8   d__Bacteria;p__Desulfobacterota_A;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio;s__

```
In default annotation mode, `DRAM` takes as only input the directory containing all the bins we would like to annotate in fasta format (either .fa or .fna). There are few parameter that can be modified if not using the default mode. Once the annotation step is done, the mode `distill` is used to summarise the obtained results. **Note:** due to the increased memory requirements, UniRef90 database is not default and the flag `–use_uniref` should be specified in order to search amino acid sequences against UniRef90. In this exercise, due to memory and time constrains, we won't be using UniRef90 database.
```
module load Python
module load DRAM

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

```
#!/bin/bash -e
#SBATCH --job-name=DRAM_annotation
#SBATCH --account=ga02676
#SBATCH --time=6:00:00
#SBATCH --mem=25Gb
#SBATCH -e slurm-DRAM_annot.%A-%a.err #Standard error
#SBATCH -o slurm-DRAM_annot.%A-%a.out #Standard output
#SBATCH --mail-type ALL
#SBATCH --mail-user <your_email>


module load Miniconda3/4.7.10
source activate /nesi/nobackup/ga02676/Metagenomics_summerschool/carmen/00.DRAM/DRAM_env #Modify it when Dini gets DRAM set-up
module load Python

DRAM.py annotate -i '/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.annotation/example_data/*.fna' --checkm_quality /path/.../checkm.txt --gtdb_taxonomy /path/.../gtdbtk.bac120.classification_pplacer.tsv  -o annotation 

## End of slurm script

```

The program will take 4-4.5 hours to run, so we will submit the jobs and inspect the results tomorrow morning. 



