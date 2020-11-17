# Gene annotation (part 1)

### Objectives

* [**BLAST**-like gene annotations and **domain** annotations](#blast-like-gene-annotations-and-domain-annotations)
* [Overview of annotation databases](#blast-like-gene-annotations-and-domain-annotations)
* [Evaluating the quality of gene assignment](#evaluating-the-quality-of-gene-assignment)
* [*Discussion: Differences in taxonomies (GTDB, NCBI etc)*](#differences-in-taxonomies)

---

### *BLAST*-like gene annotations and domain annotations

Broadly speaking, there are two ways we perform gene annotations with protein sequences. Both compare our sequences of interest against a curated set of protein sequences for which function is known, or is strongly suspected. In each case, there are particular strenths to the approach and for particular research questions, one option may be favoured over another.

#### BLAST-like annotation

The first of these is the `BLAST` algorithm for sequence alignment. This approach performs pairwise alignment between the gene of interest (query sequence) and the sequences in the database (target sequence). `BLAST` searches each potential target sequence for *k*-mers identified in the query sequence. Where these *k*-mers are found in targets, the ends are extended out to try to create longer regions of highly similar sequence spans. Across this span, the tool identifies the longest span of characters (nucleotide or amino acid) that match within a scoring framework to return the length of the region (coverage) and the sequence identity over the span (identity).

The original tool for performing this kind of analysis was the `BLAST` tool. While `BLAST` and its variants are still excellent tools for performing this kind of sequence annotation, they suffer from a slow runtime speed due to the need to test each query sequence against every target sequence in the database. For this reason, several tools have been published which take the basic approach of `BLAST`, but augment it with methods to reduce the number of pairwise comparisons needed to identify targets with high sequence similarity to the query. Two popular pieces of software are the tools `usearch` [here](http://www.drive5.com/usearch/) and `diamond` [here](https://github.com/bbuchfink/diamond).

#### HMM-profiling of domains

An alternate method for attributing function to query sequences it to consider them as a collection of independently functioning protein folding domains. This is the approach used in the [HMMer](http://hmmer.org/) software, and the *Pfam*, *TIGRfam*, and *PANTHER* databases. In these analyses, the database consists not of individual sequences, but of Hidden Markov models built from a collection of proteins that share a common domain. These profiles build out a statistical map of the amino acid transitions (from position to position), variations (differences at a position), and insertions/deletions between positions in the domain across the different observations in the training database and apply these maps to the query data.

These exercises will take place in the `10.gene_annotation/` folder.

---

### Annotating MAGs against the *UniProt* database with *diamond*

For this exercise we are going to use diamond for performing our annotation. We have chosen to use this tool because it is faster than BLAST, and usearch comes with licencing restrictions that make it hard to work with in a shared computing environment like NeSI.

For this exercise we have created a diamond-compatible database from the 2018 release of the UniProt database.

For input files, the `predictions/` results from the previous gene prediction exercise have been copied over to `10.gene_annotation/predictions/`.

In general, diamond takes a simple pair of input files - the protein coding sequences we wish to annotate and the database we will use for this purpose. There are a few parameters that need to be tweaked for obtaining a useful output file, however.

```bash
module purge
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

There are two output formats we can chose from which are useful for our analysis. We will obtain our output in the BLAST tabular format, which provides the annotation information in a simple-to-parse text file that can be viewed in any text or spreadsheet viewing tool. This will allow us to investigate and evaluate the quality of our annotations. 

Awkwardly, `diamond` does not provide the headers for what the columns in the output table mean. [This table](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) is a handy reference for how to interpret the output.

From here we can view important stastics for each query/target pairing such as the number of identify residues between sequences and the aligned length between query and target.

Lets set up a slurm job to annotate each of our MAGs. 

Create a new script

```bash
nano annotate_uniprot.sl
```

Paste in the script (update `<YOUR FOLDER>`)

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J annotate_uniprot
#SBATCH --res SummerSchool
#SBATCH --time 02:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e annotate_uniprot_dmnd.err
#SBATCH -o annotate_uniprot_dmnd.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load DIAMOND/0.9.25-gimkl-2018b

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/

mkdir -p gene_annotations/

for prot_file in predictions/*.genes.no_metadata.faa;
do
  out_file=$(basename ${prot_file} .faa)

  diamond blastp -p 20 --max-target-seqs 5 --evalue 0.001 \
        --db /nesi/nobackup/nesi02659/MGSS_resources_2020/databases/uniprot_nr_200213.diamond \
        -q ${prot_file} --outfmt 6 -o gene_annotations/${out_file}.uniprot.txt
done

```

Submit the script

```bash
sbatch annotate_uniprot.sl
```

---

### Annotating MAGs against the *Pfam* database with *hmmer*

The standard software for performing HMM-profiling annotation is [hmmer](http://hmmer.org/). Compared to `BLAST`, `FASTA`, and other sequence alignment and database search tools based on older scoring methodology, `HMMER` aims to be significantly more accurate and more able to detect remote homologs because of the strength of its underlying mathematical models. In the past, this strength came at significant computational expense, but in the new `HMMER3` project, `HMMER` is now essentially as fast as `BLAST`. 

`HMMER` will search one or more profiles against a sequence database for sequence hommologs, and for making sequence alignments, implementing profile hidden Markov models. In this exercise, we will perform a search using `hmmsearch`. For each profile in *hmmfile*, `HMMER` uses that query profile to search the target database of sequences indicated in *seqdb*, and output ranked lists of the sequences with the most significant matches to the profile. `hmmsearch` accepts any *fastA* file as target database input. It also accepts EMBL/UniProtKB text format, and Genbank format. It will automatically determine what format your file is in so you don’t have to specify it. 

As we did with `diamond`, we will also have to modify some parameters to get the desired ouotput. 

```bash
module load HMMER/3.1b2-gimkl-2017a

hmmsearch -h 

# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Usage: hmmsearch [options] <hmmfile> <seqdb>

# Basic options:
#  -h : show brief help on version and usage

# Options directing output:
# ...
# --tblout <f>     : save parseable table of per-sequence hits to file <f>
# ....

# Options controlling reporting thresholds:
# ...
# -E <x>     : report sequences <= this E-value threshold in output  [10.0]  (x>0)
# ...

# Other expert options:
# ...
# --cpu <n>     : number of parallel CPU workers to use for multithreads
# ...
```

We are now going to submit another slurm job to annotate our MAGs using the [Pfam database](https://pfam.xfam.org/). Matching sequences to a `Pfam` entry allows us to transfer the functional information from an experimentally characterised sequence to uncharacterised sequences in the same entry. `Pfam` then provides comprehensive annotation for each entry.

Create a new script

```bash
nano annotate_pfam.sl
```

Paste in the script (update `<YOUR FOLDER>`)

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J annotate_pfam
#SBATCH --res SummerSchool
#SBATCH --time 02:00:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH -e annotate_pfam_hmm.err
#SBATCH -o annotate_pfam_hmm.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load HMMER/3.1b2-gimkl-2017a

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation/

for prot_file in predictions/*.genes.no_metadata.faa;
do
  out_file=$(basename ${prot_file} .faa)

  hmmsearch --tblout gene_annotations/${out_file}.pfam.txt -E 1e-3 --cpu 10 /nesi/nobackup/nesi02659/MGSS_resources_2020/databases/Pfam-A.hmm ${prot_file}
  
done

```

Submit the script

```bash
sbatch annotate_pfam.sl
```

---

### Evaluating the quality of gene assignment

Determining how trustworthy a gene annotation is can be a very tricky process. How similar do protein sequences need to be to perform the same function? The answer is surprisingly low. A bioinformatic analysis performed in 1999 identified that proteins with as little as 20 - 35% sequence identity can still share the same function ([Rost, 1999](https://doi.org/10.1093/protein/12.2.85)), but this is not a universal occurrence. When evaluating annotations, consider the following questions:

1. What is the amino acid identity along the aligned region?
1. What is the amino acid *similarity* between the aligned region?
1. What is the coverage as a percentage or the query and target genes?
1. If we infer a phylogeny of this query gene with references from the target family, is a stable tree resolved?
1. Does the inclusion of this gene function make sense in the context of the organism's taxonomy?
1. Does the gene sit on a long contig that is core to the MAG, or is it a short contig carrying only a single gene?
1. If we are uncertain of a particular annotation, does the predicted gene occur in an operon? If so, are the other genes present in the annotation?

We must also remain aware of the potential for incorrectly annotated genes in the annotation database and that proteins can perform multiple functions (and may therefore be attributed multiple, inconsistent annotations). Furthermore, it is also important to consider exactly which part of the target gene the alignment is happening across. There are several catalytic centers of enzymes, such as the Fe-S cofactor, which are shared across many different proteins, and if your annotation is only spanning one of these regions then it may simply be the case that you are identifying a generic electron accepting or donating domain.

---

### Differences in taxonomies

Another way to determine if a annotation 'belongs' in the MAG of interest is to consider the predicted taxonomy of the query gene with that of the MAG itself. For example, if you detect a *Desulfovibrio*-like *dsrA* seuqence in a bin that has been classified as belonging to the genus *Desulfovibrio* then it is probably a safe bet that the annotation is correct.

However, when comparing taxonomic assignments, it is important to be aware of the differing taxonomic schemas that are circulating in the microbiological and bioinformatic literature and to know how to reconcile their differences. Similar to how the 16S rRNA gene taxonomies provided by *SILVA*, *Greengenes*, and *RDP* taxonomies all differ in some aspects, there are multiple competing taxonomies in protein databases.

#### Examples of various genome-level and protein taxonomies

* [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy)
* [Genome Taxonomy Database](https://gtdb.ecogenomic.org/)

This problem exists despite the existance of a formal [Code](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.000778) for the naming of bacteria and archaea, because 

1. There are no rules governing how we define the grouping of these names together, other than for type species
1. Defunct synonyms and basonyms are not correctly purged from taxonomy lists (this is quite noticable with the NCBI taxonomy)
1. Valid names cannot be assigned for uncultivate organisms, meaning there are many informal placeholder names in the literature. For example, clades like WPS-2, SAR324, and SAUL are widely cited in the literature despite having no official standing

It is therefore important to periodically sanity check your taxonomic annotations in order to avoid splitting taxa based on spelling differences or the use of historic names that have since been reclassified.

---




