# Gene annotation

### Objectives

* **BLAST**-like gene annotations and domain annotations
* Overview of annotation databases
* Evaluating the quality of gene assignment
* *Discussion: Differences in taxonomies (GTDB, NCBI etc)*

---

### *BLAST*-like gene annotations and domain annotations

Broadly speaking, there are two ways we perform gene annotations with protein sequences. Both compare our sequences of interest against a curated set of protein sequences for which function is known, or is strongly suspected. In each case, there are particular strenths to the approach and for particular research questions, one option may be favoured over another.

#### BLAST-like annotation

The first of these is the `BLAST` algorithm for sequence alignment. This approach performs pairwise alignment between the gene of interest (query sequence) and the sequences in the database (target sequence). `BLAST` searches each potential target sequence for *k*-mers identified in the query sequence. Where these *k*-mers are found in targets, the ends are extended out to try to create longer regions of highly similar sequence spans. Across this span, the tool identifies the longest span of characters (nucleotide or amino acid) that match within a scoring framework to return the length of the region (coverage) and the sequence identity over the span (identity).

The original tool for performing this kind of analysis was the `BLAST` tool. While `BLAST` and its variants are still excellent tools for performing this kind of sequence annotation, they suffer from a slow runtime speed due to the need to test each query sequence against every target sequence in the database. For this reason, several tools have been published which take the basic approach of `BLAST`, but augment it with methods to reduce the number of pairwise comparisons needed to identify targets with high sequence similarity to the query. Two popular pieces of software are the tools `usearch` [here](http://www.drive5.com/usearch/) and `diamond` [here](https://github.com/bbuchfink/diamond).

#### HMM-profiling of domains

An alternate method for attributing function to query sequences it to consider them as a collection of independently functioning protein folding domains. This is the approach used in the `HMMer` (http://hmmer.org/) software, and the *Pfam*, *TIGRfam*, and *PANTHER* databases. In these analyses, the database consists not of individual sequences, but of Hidden Markov models built from a collection of proteins that share a common domain. These profiles build out a statistical map of the amino acid transitions (from position to position), variations (differences at a position), and insertions/deletions between positions in the domain across the different observations in the training database and apply these maps to the query data.

---

### Annotating MAGs with against the *NCBI NR* database with *diamond*

For this exercise we are only going to use a single tool for performing our annotation. We have chosen to use `diamond` because it is faster than `BLAST`, and `usearch` comes with licencing restrictions that make it hard to work with in a shared computing environment like NeSI.

For this exercise we have created a `diamond`-compatible database from a 2016 release of the NCBI non-redundant protein sequence database. The reasons for using this particular database will become apparent in a subsequent exercise.

In generaly, `diamond` takes a simple pair of input files - the protein coding sequences we wish to annotate and the database we will use for this purpose. There are a few parameters that need to be tweaked for obtaining a useful output file, however.

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

There are two output formats we can chose from which are useful for our analysis. We will obtain our output in the *BLAST tabular* format, which provides the annotation information in a simple-to-parse text file that can be viewed in any text or spreadsheet viewing tool. This will allow us to evaluate the quality of our annotations. For now, just annotate a single bin:

```bash
cd 8.gene_annotation/

diamond blastp -p 2 --db /nesi/project/nesi02659/mg_workshop//NCBI_nr_2016.dmnd \
               --max-target-seqs 5 --evalue 0.001 \
               -q example_data/bin_0.genes.no_metadata.faa \
               --outfmt 6 -o bin_0.diamond.txt
```

Awkwardly, `diamond` does not provide the headers for what the columns in the output table mean. [This table](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) is a handy reference for how to interpret the output.

From here we can view important stastics for each query/target pairing such as the number of identify residues between sequences and the aligned length between query and target.

Before we proceed with this exercise, lets set up a slurm job to annotate each of our MAGs using the *BLAST XML* output format. We will need this for tomorrow.

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J annotate_diamond
#SBATCH --partition ga_bigmem
#SBATCH --time 02:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e annotate_diamond.err
#SBATCH -o annotate_diamond.out

module load DIAMOND/0.9.25-gimkl-2018b

cd 8.gene_annotation/

mkdir -p gene_annotations/

for prot_file in example_data/*.genes.no_metadata.faa;
do
  out_file=$(basename ${prot_file} .faa)

  diamond blastp -p 20 --max-target-seqs 5 --evalue 0.001 \
                 --db /nesi/project/nesi02659/mg_workshop//NCBI_nr_2016.dmnd \
                 -q ${prot_file} --outfmt 5 -o gene_annotations/${out_file}.nr.xml
done
```

**Make sure you have changed the *outfmt* value to 5, not 6!**

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
