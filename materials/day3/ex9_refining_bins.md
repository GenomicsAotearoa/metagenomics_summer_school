# Manually refining bins

### Objectives

* Prepare input files for `VizBin`
* Project a *t-SNE* and examine bin clusters
* Picking refined bins

---

### Prepare input files for *VizBin*

[**VizBin**](http://claczny.github.io/VizBin/) is a handy, GUI-based tool for creating ordinations of our binning data using the [t-Distributed Stochastic Neighbor Embedding (t-SNE)](https://lvdmaaten.github.io/tsne/) algorithm to project high-dimensional data down into a 2D plot that preserves clustering information. There's a really good video on [YouTube](https://www.youtube.com/watch?v=NEaUSP4YerM) that explains how the algorithm works in high-level terms, but for our purposes you can really consider it as a similar approach to a PCA or NMDS.

On its own, `VizBin` takes a set of contigs and performs the *t-SNE* projection using compositional data. We can optionally provide it files that annotate contigs as belonging to particular bins and a file that adds coverage data to be considered when clustering. Unfortuantely, at this stage `VizBin` only allows a single coverage value per contig which is not ideal for our purposes. There is an optional exercise down below that shows you how to create your own *t-SNE* projection using multiple coverage values.

We will use a small `python` script to create these input files, as colouring contigs by bin is a really effective way to spot areas that might need refinement.

*Note: We are using `python` here for the sake of keeping the workshop easy to follow. It is possible to perform the equivalent steps in `bash` using a loop similar to how we performed the `DAS_Tool` preparation.*

```bash
module load Python/3.7.3-gimkl-2018b

cd 6.bin_refinment/

cut -f1,2 example_data.txt > sample1.txt

python build_vizbin_inputs.py -o vb_sample1 -c sample1.txt example_data/*

# Make a few different versions of the .ann file with various columns removed
cut -f2 -d ',' vb_sample1.vizbin.ann > vb_sample1.vizbin.bin_only.ann
cut -f1,2 -d ',' vb_sample1.vizbin.ann > vb_sample1.vizbin.no_length.ann
```

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

---
