# Gene prediction

### Objectives

* Overview/refresher of **prodigal**
* Predicting protein coding sequences in MAGs
* Predicting protein coding sequences in unassembled reads
* Predicting RNA features and non-coding regions

---

### Overview/refresher of *prodigal*

At this stage we have recovered a numer of high quality genomes or population genomes. While there are interesting biological questions we can ask of the genomes at the DNA/organisational level, it is more likely that we are interested in the genes present in the organism.

How we predict genes in the metagenomic data varies depending on what features we are trying to detect. Most often, we are interested in putatively protein coding regions and open reading frames. For features that are functional but not not translated, such as ribosomal RNA and tRNA sequences we need to use alternate tools. When considering protein coding sequences, we avoid the use of the term 'open reading frame' (ORF). The nature of a fragmented assembly is that you may encounter a partial gene on the start or end of a contig that is a function gene, but lacks the start or stop codon due to issues with assembly or sequencing depth.

There are many software tools to predict gene sequences and in this workshop we will start with the tool **Prodigal** (PROkaryotic Dynamic Programming Genefinding ALgorithm). **Prodigal** has gone on to become one of the most popular microbial gene prediction algorithms as in incorporates modeling algorithms to profile the coding sequences within your genome and better identify the more cryptic (or partial) genes.

**Prodigal** is execellent for the following use cases:

1. Predicting protein-coding genes in draft genomes and metagenomes  
1. Quick and unsupervised execution, with minimal resource requirements
1. Ability to handle gaps, scaffolds, and partial genes 
1. Identification of translation initiation sites 
1. Multiple output formats, including either straight *fastA* files or the DNA sequence and protein translation for genes, as well as detailed summary statistics for each gene (e.g. contig length, gene length, GC content, GC skew, RBS motifs used, and start and stop codon usage)

**Prodigal** is not the best tool to use for the following cases:

1. Predicting RNA genes 
1. Handling genes with introns 
1. Deal with frame shifts

It is also not advised to use **prodigal** when making predictions through your unassembled reads. If you are working with unassembled data, **FragGeneScan** is a better tool, as it is more sensitive for partial genes and does not assume every piece of DNA in the input *fastA* file must be coding.

---

### Predicting protein coding sequences



---

### Predicting RNA features and non-coding regions

