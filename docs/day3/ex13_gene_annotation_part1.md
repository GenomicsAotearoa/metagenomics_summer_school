# Gene annotation I: BLAST-like and HMM

!!! info "Objectives"

    - [Annotation methods](#annotation-methods)
    - [Annotating MAGs against the *UniProt* database with `DIAMOND`](#annotating-mags-against-the-uniprot-database-with-diamond)
    - [Annotating MAGs against the *Pfam* database with `HMMER`](#annotating-mags-against-the-pfam-database-with-hmmer)
    - [Annotating signal peptides with `SignalP6`](#annotating-signal-peptides-of-predicted-genes-using-signalp6)
    - [Evaluating the quality of gene assignment](#evaluating-the-quality-of-gene-assignment)
    - [Differences in taxonomies](#differences-in-taxonomies)
    
---

## Annotation methods

Broadly speaking, there are two ways we perform gene annotations with protein sequences. Both compare our sequences of interest against a curated set of protein sequences for which function is known, or is strongly suspected. In each case, there are particular strengths to the approach and for particular research questions, one option may be favoured over another.

### `BLAST`-like annotation

The first of these is the `BLAST` algorithm for sequence alignment. This approach performs pairwise alignment between the gene of interest (query sequence) and the sequences in the database (target sequence). `BLAST` searches each potential target sequence for *k*-mers identified in the query sequence. Where these *k*-mers are found in targets, the ends are extended out to try to create longer regions of highly similar sequence spans. Across this span, the tool identifies the longest span of characters (nucleotide or amino acid) that match within a scoring framework to return the length of the region (coverage) and the sequence identity over the span (identity).

The original tool for performing this kind of analysis was the `BLAST` tool. While `BLAST` and its variants are still excellent tools for performing this kind of sequence annotation, they suffer from a slow runtime speed due to the need to test each query sequence against every target sequence in the database. For this reason, several tools have been published which take the basic approach of `BLAST`, but augment it with methods to reduce the number of pairwise comparisons needed to identify targets with high sequence similarity to the query. Two popular pieces of software are the tools [`USEARCH`](http://www.drive5.com/usearch/) and [`DIAMOND`](https://github.com/bbuchfink/diamond).

### HMM-profiling of domains

An alternate method for attributing function to query sequences is to consider them as a collection of independently functioning protein folding domains. This is the approach used in the [HMMER](http://hmmer.org/) software, and the *Pfam*, *TIGRfam*, and *PANTHER* databases. In these analyses, the database consists not of individual sequences, but of Hidden Markov models built from a collection of proteins that share a common domain. These profiles build out a statistical map of the amino acid transitions (from position to position), variations (differences at a position), and insertions/deletions between positions in the domain across the different observations in the training database and apply these maps to the query data.

These exercises will take place in the `10.gene_annotation_and_coverage/` folder. 

!!! terminal-2 "Navigate to working directory"

    ```
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation_and_coverage
    ```

---

## Annotating MAGs against the *UniProt* database with `DIAMOND`

For this exercise we are going to use `DIAMOND` for performing our annotation. We have chosen to use this tool because it is faster than `BLAST`, and `USEARCH` comes with licensing restrictions that make it hard to work with in a shared computing environment like NeSI.

For this exercise we have created a diamond-compatible database from the 2024 release of the *UniProt* database.

For input files, the `predictions/` results from the previous gene prediction exercise have been copied over to `10.gene_annotation_and_coverage/predictions/`.

In general, `diamond` takes a pair of input files - the protein coding sequences we wish to annotate and the database we will use for this purpose. There are a few parameters that need to be tweaked for obtaining a useful output file, however.

!!! terminal "code"

    ```bash
    module purge
    module load DIAMOND/2.1.6-GCC-11.3.0

    diamond help
    ```

!!! circle-check "Terminal output"

    ```
    diamond v2.0.15.153 (C) Max Planck Society for the Advancement of Science
    Documentation, support and updates available at http://www.diamondsearch.org
    Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)

    Syntax: diamond COMMAND [OPTIONS]

    Commands:
    ...
    blastp  Align amino acid query sequences against a protein reference database
    ...

    General options:
    --threads (-p)         number of CPU threads
    --db (-d)              database file
    --out (-o)             output file
    --outfmt (-f)          output format
    ...
    ```

There are two output formats we can chose from which are useful for our analysis. We will obtain our output in the BLAST tabular format, which provides the annotation information in a simple-to-parse text file that can be viewed in any text or spreadsheet viewing tool. This will allow us to investigate and evaluate the quality of our annotations. 

Awkwardly, `DIAMOND` does not provide the headers for what the columns in the output table mean. [This table](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) is a handy reference for how to interpret the output.

From here we can view important statistics for each query/target pairing such as the number of identical residues between sequences and the aligned length between query and target.

Before we begin, we need to create an directory for outputs.

!!! terminal "code"

    ```bash
    mkdir -p gene_annotations
    ```

Now, lets set up a slurm job to annotate each of our MAGs. 

Create a new script

!!! terminal-2 "Create script named `annotate_uniprot.sl`"

    ```bash
    nano annotate_uniprot.sl
    ```

!!! warning "Remember to update <YOUR FOLDER> to your own folder"

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e

    #SBATCH --account       nesi02659
    #SBATCH --job-name      annotate_uniprot
    #SBATCH --partition     milan
    #SBATCH --time          01:00:00
    #SBATCH --mem           20GB
    #SBATCH --cpus-per-task 20
    #SBATCH --array         0-9
    #SBATCH --error         %x_%A_%a.err
    #SBATCH --output        %x_%A_%a.out

    # Load modules
    module purge
    module load DIAMOND/2.1.6-GCC-11.3.0

    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation_and_coverage

    # Variables
    prot_file=predictions/bin_${SLURM_ARRAY_TASK_ID}.genes.no_metadata.faa
    out_file=$(basename ${prot_file} .faa)
    db=/nesi/nobackup/nesi02659/MGSS_2024/resources/databases/uniprot.20240724.dmnd

    # Run DIAMOND
    diamond blastp --threads $SLURM_CPUS_PER_TASK --max-target-seqs 5 --evalue 0.001 \
                   --db $db --query ${prot_file} --outfmt 6 \
                   --out gene_annotations/${out_file}.uniprot.txt
    ```

!!! terminal-2 "Submit the script"

    ```bash
    sbatch annotate_uniprot.sl
    ```

---

## Annotating MAGs against the *Pfam* database with `HMMER`

The standard software for performing HMM-profiling annotation is [HMMER](http://hmmer.org/). Compared to `BLAST`, `FASTA`, and other sequence alignment and database search tools based on older scoring methodology, `HMMER` aims to be significantly more accurate and more able to detect remote homologs because of the strength of its underlying mathematical models. In the past, this strength came at significant computational expense, but in the new `HMMER3` project, `HMMER` is now essentially as fast as `BLAST`. 

`HMMER` will search one or more profiles against a sequence database for sequence hommologs, and for making sequence alignments, implementing profile hidden Markov models. In this exercise, we will perform a search using `hmmsearch`. For each profile in *hmmfile*, `HMMER` uses that query profile to search the target database of sequences indicated in *seqdb*, and output ranked lists of the sequences with the most significant matches to the profile. `hmmsearch` accepts any *fastA* file as target database input. It also accepts EMBL/UniProtKB text format, and Genbank format. It will automatically determine what format your file is in so you don’t have to specify it. 

As we did with `diamond`, we will also have to modify some parameters to get the desired ouotput. 

!!! terminal "code"

    ```bash
    module load HMMER/3.3.2-GCC-11.3.0

    hmmsearch -h
    ```

??? abstract "`hmmsearch -h`" 

    ```
    hmmsearch :: search profile(s) against a sequence database
    HMMER 3.3.2 (Nov 2020); http://hmmer.org/
    Copyright (C) 2020 Howard Hughes Medical Institute.
    Freely distributed under the BSD open source license.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Usage: hmmsearch [options] <hmmfile> <seqdb>

    Basic options:
     -h : show brief help on version and usage

    Options directing output:
    ...
    --tblout <f>     : save parseable table of per-sequence hits to file <f>
    ....

    Options controlling reporting thresholds:
    ...
    -E <x>     : report sequences <= this E-value threshold in output  [10.0]  (x>0)
    ...

    Other expert options:
    ...
    --cpu <n>     : number of parallel CPU workers to use for multithreads
    ...
    ```

We are now going to submit another slurm job to annotate our MAGs using the [Pfam database](https://www.ebi.ac.uk/interpro/download/Pfam/). Pfam used to have a standalone website, but it has recently been integrated into InterPro maintained by the European Bioinformatics Institute (see [announcement](https://xfam.wordpress.com/2022/08/04/pfam-website-decommission/)). Matching sequences to a `Pfam` entry allows us to transfer the functional information from an experimentally characterised sequence to uncharacterised sequences in the same entry. `Pfam` then provides comprehensive annotation for each entry.

!!! terminal-2 "Create script named `annotate_pfam.sl`"

    ```bash
    nano annotate_pfam.sl
    ```

!!! warning "Remember to update <YOUR FOLDER> to your own folder"

!!! terminal "code"

    ```bash linenums="1"
    #!/bin/bash -e
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      annotate_pfam
    #SBATCH --partition     milan
    #SBATCH --time          00:20:00
    #SBATCH --mem           5GB
    #SBATCH --cpus-per-task 10
    #SBATCH --array         0-9
    #SBATCH --error         %x_%A_%a.err
    #SBATCH --output        %x_%A_%a.out
    
    # Load modules
    module purge
    module load HMMER/3.3.2-GCC-11.3.0
    
    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation_and_coverage
    
    # Variables
    prot_file=predictions/bin_${SLURM_ARRAY_TASK_ID}.genes.no_metadata.faa
    out_file=$(basename ${prot_file} .faa)
    db=/nesi/nobackup/nesi02659/MGSS_resources_2022/databases/Pfam-A.hmm
    
    # Run HMMER
    hmmsearch --tblout gene_annotations/${out_file}.pfam.txt -E 0.001 \
              --cpu $SLURM_CPUS_PER_TASK \
              ${db} ${prot_file}
    ```

!!! terminal-2 "Submit the script"

    ```bash
    sbatch annotate_pfam.sl
    ```

---

## Annotating signal peptides of predicted genes using `SignalP6`

In addition to putative function, we can also predict whether a protein is secreted or not using `SignalP`. In the 5<sup>th</sup> version, users needed to know the Gram status of the organism for the software to predict signal peptides. In 2022, the developers ([Teufel et al., 2022](https://www.nature.com/articles/s41587-021-01156-3)) produced a 6<sup>th</sup> version that uses protein language models (like ChatGPT, for proteins!) to predict the presence of signal peptides across the domains of life.

`SignalP6` has several models depending on required accuracy and hardware availability. Here, we are running the fast model on CPU. There are models built for GPU (they are interconvertible and the software package comes with utilities for easy conversion) that is blazingly fast (most language models require GPUs for accelerations). However, GPUs on NeSI are a coveted resource. For our small dataset, the fast model on CPU is sufficient as an example.

!!! note "`SignalP6` is licensed software"

    In terms of computational capability to detect signal peptides, there is nothing like SignalP6. It is free for academic use and other models can be requested from the developers. If you are not running this on NeSI, you will need to acquire your own license.

!!! terminal-code "code"

    ```bash
    nano annotate_signal_peptide.sl
    ```

!!! terminal-code "code"

    ```bash linenums="1"
    #!/bin/bash -e

    #SBATCH --account       nesi02659
    #SBATCH --job-name      signal_peptide
    #SBATCH --partition     milan
    #SBATCH --time          00:15:00
    #SBATCH --mem           10GB
    #SBATCH --cpus-per-task 24
    #SBATCH --error         %x_%A_%a.err
    #SBATCH --output        %x_%A_%a.out
    #SBATCH --array         0-9

    # Load modules
    module purge
    module load SignalP/6.0g-gimkl-2022a-Python-3.10.5

    # Working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation_and_coverage

    # Array variables
    files=(predictions/*.faa)
    bin=${files[$SLURM_ARRAY_TASK_ID]}
    prefix=$(basename ${bin} genes.no_metadata.faa)

    # Run SignalP
    signalp6 --fastafile ${bin} --output_dir signal_peptides/${prefix} \
             --format none --organism other --bsize 50 \
             --write_procs $SLURM_CPUS_PER_TASK
    ```

!!! warning "`SignalP` will probably `TIMEOUT`"

    The above script will almost, inevitably, time out before completing every MAG. However, that is not critical for our purposes today. We simply need at least one of the array jobs to complete (several will) for us to explore the output.

The output of the script above are directories per bin that contain 5 files:

* `predicted_results.txt` the main output with gene ID and the predicted signal peptide ("Other" means none predicted), the prediction probability (i.e., confidence), and the cleavage site
* `output.gff3` A GFF3-format output of sequences with start and end positions of signal peptides relative to each amino acid sequence
* `output.json` Similar to above, with more information about run parameters
* `processed_entries.fasta` Amino acid sequences with the signal peptide cleaved off (i.e., predicted mature proteins)
* `region_output.gff3` Signal peptides have specific regions characteristic of each type of signal peptide. These usually inform us about cleavage sites for those proteins. This output provides this information. 

!!! warning "`--format none`"

    If you choose to run the above script with other formats (e.g., txt, png, all), per-sequence data outputs and prediction plots (if requested) will be produced. For large numbers of amino acid sequences, this can be problematic as it uses up a lot of space. Unless you're needing to query specific parts of a few proteins, it is advised to stick to "none" to output summaries only.


## Evaluating the quality of gene assignment

Determining how trustworthy a gene annotation is can be a very tricky process. How similar do protein sequences need to be to perform the same function? The answer is surprisingly low. A bioinformatic analysis performed in 1999 identified that proteins with as little as 20 - 35% sequence identity can still share the same function ([Rost, 1999](https://doi.org/10.1093/protein/12.2.85)), but this is not a universal occurrence. When evaluating annotations, consider the following questions:

!!! quote ""
    1. What is the amino acid identity along the aligned region?
    1. What is the amino acid *similarity* between the aligned region?
    1. What is the coverage as a percentage of the query and target genes?
    1. If we infer a phylogeny of this query gene with references from the target family, is a stable tree resolved?
    1. Does the inclusion of this gene function make sense in the context of the organism's taxonomy?
    1. Does the gene sit on a long contig that is core to the MAG, or is it a short contig carrying only a single gene?
    1. If we are uncertain of a particular annotation, does the predicted gene occur in an operon? If so, are the other genes present in the annotation?

We must also remain aware of the potential for incorrectly annotated genes in the annotation database and that proteins can perform multiple functions (and may therefore be attributed multiple, inconsistent annotations). Furthermore, it is also important to consider exactly which part of the target gene the alignment is happening across. There are several catalytic centers of enzymes, such as the Fe-S cofactor, which are shared across many different proteins, and if your annotation is only spanning one of these regions then it may simply be the case that you are identifying a generic electron accepting or donating domain.

---

## Differences in taxonomies

Another way to determine if an annotation 'belongs' in the MAG of interest is to consider the predicted taxonomy of the query gene with that of the MAG itself. For example, if you detect a *Desulfovibrio*-like *dsrA* sequence in a bin that has been classified as belonging to the genus *Desulfovibrio* then it is probably a safe bet that the annotation is correct.

However, when comparing taxonomic assignments, it is important to be aware of the differing taxonomic schemas that are circulating in the microbiological and bioinformatic literature and to know how to reconcile their differences. Similar to how the 16S rRNA gene taxonomies provided by *SILVA*, *Greengenes*, and *RDP* taxonomies all differ in some aspects, there are multiple competing taxonomies in protein databases.

This problem exists despite the existence of a formal [Code](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.000778) for the naming of bacteria and archaea, because 

!!! quote ""
    1. There are no rules governing how we define the grouping of these names together, other than for type species
    1. Defunct synonyms and basonyms are not correctly purged from taxonomy lists (this is quite noticeable with the NCBI taxonomy)
    1. Valid names cannot be assigned for uncultivated organisms, meaning there are many informal placeholder names in the literature. For example, clades like WPS-2, SAR324, and SAUL are widely cited in the literature despite having no official standing

It is therefore important to periodically sanity check your taxonomic annotations in order to avoid splitting taxa based on spelling differences or the use of historic names that have since been reclassified.

---




