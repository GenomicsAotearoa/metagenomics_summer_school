# Manually refining bins

### Objectives

* Prepare input files for `VizBin`
* Project a *t-SNE* and examine bin clusters
* Refining bins by removing incorrectly assigned contigs
* *Optional:* Creating new `VizBin` profiles with different fragment lengths
* *Optional:* Scripts for processing data with `ESOMana`
* Assigning taxonomy to the refined bins

---

### Prepare input files for *VizBin*

[**VizBin**](http://claczny.github.io/VizBin/) is a handy, GUI-based tool for creating ordinations of our binning data using the [t-Distributed Stochastic Neighbor Embedding (t-SNE)](https://lvdmaaten.github.io/tsne/) algorithm to project high-dimensional data down into a 2D plot that preserves clustering information. There's a really good video on [YouTube](https://www.youtube.com/watch?v=NEaUSP4YerM) that explains how the algorithm works in high-level terms, but for our purposes you can really consider it as a similar approach to a PCA or NMDS.

On its own, `VizBin` takes a set of contigs and performs the *t-SNE* projection using compositional data. We can optionally provide it files that annotate contigs as belonging to particular bins and a file that adds coverage data to be considered when clustering. Unfortuantely, at this stage `VizBin` only allows a single coverage value per contig which is not ideal for our purposes. There is an optional exercise down below that shows you how to create your own *t-SNE* projection using multiple coverage values.

We will use a small `python` script to create these input files, as colouring contigs by bin is a really effective way to spot areas that might need refinement.

*Note: We are using `python` here for the sake of keeping the workshop easy to follow. It is possible to perform the equivalent steps in `bash` using a loop similar to how we performed the `DAS_Tool` preparation.*

```bash
module load Python/3.7.3-gimkl-2018b

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/6.bin_refinment/

# Step 1
cut -f1,4 example_data_20k.txt > sample1.txt

# Step 2
python build_vizbin_inputs.py -o vb_sample1 -c sample1.txt example_data_20k/*.fna

# Step 3
# Make a few different versions of the .ann file with various columns removed
cut -f2 -d ',' vb_sample1.vizbin.ann > vb_sample1.vizbin.bin_only.ann
cut -f1,2 -d ',' vb_sample1.vizbin.ann > vb_sample1.vizbin.no_length.ann
```

There are a few things going on in this set of commands.

##### Step 1

Here we are taking the file `example_data.txt`, which contains the per-sample coverage scores for each sub-contig in the data set. Because `VizBin` only uses coverage as a means to modify the visualisation, not the ordination itself, we can only pass the coverage values from a single sample at a time.

##### Step 2

The `build_vizbin_inputs.py` script is what produces our input files for `VizBin`. It the folder of bin files as the input, and optionally a table of sub-contig coverage values, which is defined by the `-c` parameter.

What this script is doing is taking each fasta file and picking out the names of the contigs found in that file (bin). It is then looking for any coverage information which is associated with that contig in the `sample1.txt` file, and adding that as a new column to the file. The resulting information is written to a file that carries the name from the `-o` parameter. This file is a comma-delimited table (*csv file*) that presents the information in the way that `VizBin` expects it.

```bash
head -n5 vb_sample1.vizbin.ann
# coverage,label,length
# 4.14625,bin_0,20000
# 3.5722400000000003,bin_0,20000
# 3.8721900000000002,bin_0,20000
# 3.80796,bin_0,20000
```

The order of rows in this file corresponds to the order of contigs in the *fastA* file that is produced by the script - `vb_sample1.vizbin.fna`.

##### Step 3

Here we are just producing a few variations of the *.ann* file, with various columns removed.

---

### Project a *t-SNE* and examine bin clusters

We can now use these files in `VizBin` to curate the contigs in our bins. We will load and view the data in a few different steps.

#### Loading the input files

To get started, simply click the 'Choose...' button then navigate to the *fastA* file.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_load_fasta.PNG)

Once this is imported, using the 'Show additional options' button to expose the advanced options, and add your 'bin only' *.ann* file into the 'Annotation file (optional)' field.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_load_ann.PNG)

#### Executing the *t-SNE*

For now leave all other parameters as default. Click the 'Start' button to begin building the ordination. When it completes, you should see an output similar to the following:

##### Contigs coloured by bin

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_bin_only.PNG)

##### Contigs coloured by bin, sized by length, shaded by coverage

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_all_ann.PNG)

Similar to any other projection technique, we interpret the closeness of points as a proxy for how similar they are, and because of our *.ann* file we can see which contigs belong to the same bin.

---

### Picking refined bins

We can use the interactive GUI to pick the boundaries of new bins, or to identify contigs which we do not believe should be retained in the data. Have a play around with the interface, testing out the following commands:

1. Left-click and drag: Highlight an area of the ordination to zoom into
1. Right-click, 'Zoom Out', 'Both Axes': Rest the view of view
1. Left-click several points: Create a selection of contigs to extract from the data
1. Right-click, 'Selection', 'Export': Save the selected contigs into a new file
1. Right-click, 'Selection', 'Clear selection': Clear the current selection

How you proceed in this stage is up to you. You can either select bins based on their boundary, and call these the refined bins. Alternatively, you could select outlier contigs and examine these in more detail to determine whether or not they were correctly placed into the bin. For example:

##### Highlight an area with problematic contigs

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_select_to_zoom.PNG)

##### Select the contigs to examine

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_select_outlier.PNG)

##### Export the contigs

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex10_export.PNG)

Which way you proceed really depends on how well the ordination resolves your bins, and it might be that both approaches are needed.

Use a combination of `VizBin` and `seqmagick` to remove contigs from bins when you do not trust the placement of the contig. We are aiming to reduce each bin to a trusted set of contigs.

---

### *Optional:* Creating new *VizBin* profiles with different fragment lengths

The data you have been working with was created using the `cut_up_fasta.py` script that comes with the binning tool `CONCOCT`. It was run to cut contigs into 20k fragments, to better add density to the cluster. If you would like to visualise the data using different contig fragment sizes, you can create these using the following commands:

```bash
module load CONCOCT/1.0.0-gimkl-2018b-Python-2.7.16

mkdir custom_chop/

for bin_file in example_data_unchopped/*;
do
    bin_name=$(basename ${bin_file} .fna)
    cut_up_fasta.py -c YOUR_CONTIG_SIZE -o 0 --merge_last ${bin_file} > custom_chop/${bin_name}.chopped.fna
done
```

You can then create an *.ann* file using the same `python` script as above, but we will not be able to add in the coverage information for this run. Because the `build_vizbin_inputs.py` script is written in version 3 of `python`, but the `CONCOCT` module loads version 2, we need to reload `python 3.7` before running the script again.

```bash
module load Python/3.7.3-gimkl-2018b

python build_vizbin_inputs.py -o vb_sample1 custom_chop/*
```

---

### *Optional:* Scripts for processing data with *ESOMana*

A suite of tools for creating input files for `ESOMana` can be found on github [here](https://github.com/tetramerFreqs/Binning).

The tool `ESOMana` can be downloaded from [SourceForge](http://databionic-esom.sourceforge.net/).

---

### Assigning taxonomy to the refined bins

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
#SBATCH --time 2:00:00
#SBATCH --mem 120GB
#SBATCH --cpus-per-task 10
#SBATCH -e gtdbtk_test.err
#SBATCH -o gtdbtk_test.out

module load GTDB-Tk/0.2.2-gimkl-2018b-Python-2.7.16

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/6.bin_refinment/

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
