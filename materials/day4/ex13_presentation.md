# Presentation of data

### Objectives

* Overview and installing `R` packages
* Build a basic heatmap from `BLAST` data using `R`
* Build a gene alignment figure from annotation and `prodigal` data using `R`

---

### Overview and installing *R* packages

There are a number of powerful packages within the `R` software environment which can be used to create high quality images of genomic and metagenomic data. While each of these packages comes with its own documentation, these documents and tutorials usually assume that your data is already in some already-prepared format. Our data will almost never be in this format, though, so these exercises show two brief examples of how we can scrape data from our existing files to create useful figures. As such, these examples are more complicated than what you would get reading the tutorials and manuals of the plotting tools, but will be transferable to your own work.

Both exercises will take place with files in the `9.data_presentation/` folder.

Before we get started, you will need to load a module for the `R` software, and then install a few packages for later use i the tutorial. There are several versions of `R` installed on NeSI, so we will use the newest version available, 3.6.1.

```
module load R/3.6.1-gimkl-2018b
```

We can then start `R` in interactive mode by simply typing the command `R` into the terminal:

```
R
```

Once loaded we need to install three packages. This is done with the `install.packages()` command:

```R
install.packages('genoPlotR')
install.packages('tidyverse')
install.packages('pheatmap')
```

When you install the first package, you might be prompted to select a download mirror. Just choose option 1 to proceed. If you receive any prompts or errors regarding the library/location to install packages, just accept all prompts to create a local R library in your home directory.

---

###  Build a basic heatmap from annotation data using *R*

For this exercise, navigate into the `cazy_heatmap/` directory. These annotation files we will be using have been pre-computed by annotating the `prodigal` gene predictions against the [CAZy](http://www.cazy.org/) database using the [dbCAN](http://bcb.unl.edu/dbCAN2/) resource. Briefly, each annotation was made by:

1. Annotating each `prodigal` output against the **dbCAN** database using `hmmer`
1. Converting the raw `hmmer` output into a table using the `hmmscan-parser.py` script that bundles with **dbCAN**

#### Import the data into an R data.frame

The first thing we need to do is import these annotations into `R`. We will do this using the following workflow

1. Import `tidyverse` libraries for manipulating data tables
1. Create an empty data.frame (master data.frame)
1. Importing each of the text file into an `R`, then 
   - Append a new column to the text file, consisting of the name of the bin
   - Select only the bin name and gene annotation from the imported table
   - Append the text table to the master data.frame

First, we import our `R` libraries with the `library()` command. For this workflow, we need three libraries from the `tidyverse` package:

```R
library(dplyr)
library(tidyr)
library(tibble)
```

We can then import our data using the `list.files()` function to loop over each text file, and the `mutate`, `select`, and pipe (`%>%`) functions from the `dplyr` library.

```R
cazy_files = list.files('.')

# For each file, import it, drop unneeded columns, and add a column recording the bin name
cazy_df = data.frame()

for( cazy_file in cazy_files ) {
    df = read.table(cazy_file, sep='\t', stringsAsFactors=F, header=F) %>% 
         mutate('Bin' = cazy_file) %>%
         select(CAZy=V3, Bin)
    cazy_df = rbind(cazy_df, df)
}
```

We can inspect the final data.frame using the `head` command:

```R
head(cazy_df)
#         CAZy       Bin
# 1    AA3.hmm bin_0.txt
# 2 GH13_9.hmm bin_0.txt
# 3  GH104.hmm bin_0.txt
# 4    GT5.hmm bin_0.txt
# 5  GH116.hmm bin_0.txt
# 6   GT41.hmm bin_0.txt
```

We can also confirm that we have imported all of the text files by looking at the unique entries in the `Bin` column:

```R
unique(cazy_df$Bin)
# [1] "bin_0.txt" "bin_1.txt" "bin_2.txt" "bin_3.txt" "bin_4.txt" "bin_5.txt"
# [7] "bin_6.txt" "bin_7.txt" "bin_8.txt" "bin_9.txt"
```

We will now perform a summarising step, aggregating instances of multiple genes with the same annotation into a single count for each genome. We do this by

- For each bin in the data.frame
  - For each annotation in the bin
    - Count the number of times the annotation is observed

For the majority of cases this will probably be one, but there will be a few cases where multiple annotations have been seen.

This process is done using the `group_by` and `tally` functions from `dplyr`, again using pipes to pass the data between functions.

```R
cazy_counts = cazy_df %>% group_by(Bin, CAZy) %>% tally()

cazy_counts
# A tibble: 253 x 3
# Groups:   Bin [10]
#   Bin       CAZy            n
#   <chr>     <chr>       <int>
# 1 bin_0.txt AA3.hmm         1
# 2 bin_0.txt AA7.hmm         1
# 3 bin_0.txt CE11.hmm        1
# 4 bin_0.txt CE9.hmm         1
# 5 bin_0.txt GH100.hmm       1
# 6 bin_0.txt GH104.hmm       1
# 7 bin_0.txt GH116.hmm       1
# 8 bin_0.txt GH13_11.hmm     1
# 9 bin_0.txt GH13_20.hmm     1
#10 bin_0.txt GH13_9.hmm      1
# … with 243 more rows
```

We now have a data.frame-like object (a [tibble](https://tibble.tidyverse.org/)) with three columns. We can convert this into a gene matrix using the `pivot_wider` function from the `tidyr` library to create a genome x gene matrix in the following form:

|Bin|CAZy_1|CAZy_2|...|CAZy_n|
|:---:|:---:|:---:|:---:|:---:|
|bin_0|N. of genes|N. of genes|...|N. of genes|
|bin_1|N. of genes|...|...|...|
|...|...|...|...|...|
|bin_9|N. of genes|...|...|...|

```R
cazy_matrix = cazy_counts %>% pivot_wider(id_cols=Bin, names_from=CAZy, values_from=n, values_fill=list(n = 0))
```

Finally, we create the actual plot by passing this matrix into the `pheatmap` library. Before doing this, we need to take the text column `Bin` from the matrix and move it into the rownames, as this is how `pheatmap` infers the names of our samples. Also, if we left text data in the numeric input for a heatmap, the function would crash. We can quickly transfer the `Bin` column to the rownames using the `column_to_rownames` function from the `tibble` library.

```R
library(pheatmap)

cazy_matrix %>% column_to_rownames('Bin') %>% as.matrix(.) %>% pheatmap(.)
```

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex14_heatmap.PNG)

And there we go. This is a pretty basic heatmap, so there are a number of presentation issues with it. If you have time, try to do the following fixes to the heatmap by exploring the manual for `pheatmap` or other `tidyverse` and `R` functions.

1. Replace the column-wise clustering with alphabetic arrangement of the gene columns
1. Change the colour palette of the heatmap
1. Reduce the font size in the column (gene) labels
1. Remove the `.hmm` extensions from the gene names, and the `.txt` extensions from the bin names
1. Add colouring to the row (bins), marking bins belonging to the same phyla

---

###  Build a gene alignment figure from annotation and *prodigal* data using *R*

When investigating the evolution of genomes, we sometimes want to consider not only the presence/absence of genes in a genome, but also how they are arranged in an operon. For this next exercise, we are going to visualise several flagella genes from bin_4 and bin_5, comparing their arrangement between the two organisms.

For this exercise, navigate to the folder `flagella_synteny/`. Here, you have been provided with a copy of the `prodigal` gene predictions for each of the bins (`.faa` files) and a small table of the annotation of some key genes of interest (`.txt` files). This later files were created by manually searching the annotations obtained in the previous exercises.

#### Part 1 - Parsing files in *bash*

We will be performing this exercise in two stages. Firstly, in `bash`, we will use `cut` and `tail` to pull out the genes of interest from the `prodigal` files. The gene names will then be used to create a table of gene coordinates from the `prodigal` output using `grep`, `cut`, and `sed`.

```bash
cut -f3 bin_4.txt | tail -n+2 > bin_4.genes
cut -f3 bin_5.txt | tail -n+2  > bin_5.genes
```

We now have two new files, ending with the `.genes` suffix which are simply a list of the genes that we wish to extract from the `prodigal` files. We can then use each of these files as input in for `grep` to pull out the *fastA* entries that correspond to these genes.

```bash
grep -f bin_4.genes bin_4.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_4.coords
grep -f bin_5.genes bin_5.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_5.coords
```

As previously, we will quickly go through the steps of this command:

```bash
grep -f bin_4.genes bin_4.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_4.coords
|                                     |                      |              |                |
|                                     |                      |              |                Write the output
|                                     |                      |              |
|                                     |                      |              Replace each '#' character with tab
|                                     |                      |
|                                     |                      For each line, replace the '>' with empty text
|                                     |
|                                     Split the results into columns delimited by the # character,
|                                     then take columns 1 - 4.
|
Select the lines of the bin_4.genes.faa file that contain entries found in the bin_4.genes file.
```

Check the content of the `.coords` files now. You should see something like the following:

```bash
cat bin_4.coords
# 4f94130388f7df8dfba7383a12b9de4d_313     332390          333061          1
# 4f94130388f7df8dfba7383a12b9de4d_314     333076          334308          1
# 4f94130388f7df8dfba7383a12b9de4d_315     334348          335091          1
# 4f94130388f7df8dfba7383a12b9de4d_316     335159          335941          1
# 4f94130388f7df8dfba7383a12b9de4d_317     335979          336773          1
# 4f94130388f7df8dfba7383a12b9de4d_318     336795          337907          1
```

If you recall from the previous exercise on gene prediction, we have taken the first four entries from each line of the `prodigal` output, which consists of:

1. The gene name, written as [CONTIG]\_[GENE]
1. The start position of the gene
1. The stop position of the gene
1. The orienation of the gene

We will now use these tables, together with the annotation tables to create the gene synteny view in `R`.

#### Part 2 - Producing the figure in *R*

There are two libaries we need to import for this exercise,

```R
library(dplyr)
library(genoPlotR)
```

##### Part 2.1 - Creating dna_seq tables

We can now begin importing the data. First, we will import the `.coords` files, and set column names to the files.

```R
bin_4_coords = read.table('bin_4.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_4_coords) = c('name', 'start', 'end', 'strand')

bin_5_coords = read.table('bin_5.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_5_coords) = c('name', 'start', 'end', 'strand')
```

Take a look at the content of each of these data.frames, by entering their names into the terminal. You should notice that the coordinates occur at quite different positions between the genomes. If we were looking at compelte genomes, this would indicate their position relative to the origin of replication, but as these are unclosed genomes obtained from MAGs, they reflect the coordinates upon their particular contig.

For the purposes of our analysis, the large difference between coordinate positions is not important to us as we are not looking at genome-level organisation, only the relative positions of these genes within their local context. We will therefore reduce some of the leading digits off these numbers, so that they occur in the same scale.

```R
bin_4_coords$start = bin_4_coords$start - 330000
bin_4_coords$end = bin_4_coords$end - 330000

bin_5_coords$start = bin_5_coords$start - 1200000
bin_5_coords$end = bin_5_coords$end - 1200000
```

We now parse these data.frames into the *dna_seq* data class, which is defined by the `genoPlotR` library.

```R
bin_4_ds = dna_seg(bin_4_coords)
bin_5_ds = dna_seg(bin_5_coords)
```

##### Part 2.2 - Creating a comparison table

In order to plot the data, we must now create a table that links the start and stop positions between the two *dna_seq* objects. To do this, we will import the `.txt` files from before, then use several functions from the `tidyverse` libraries to create a single matched table.

Import the data, and append the annotation and coordinate data. This is done using two functions from base `R` - `read.table` and `cbind` (column bind, to merge two data.frames). We then pass the resulting data.frame through the `select` function to extract only the columns we are interested in.

```R
bin_4_ann = read.table('bin_4.txt', header=T, sep='\t', stringsAsFactors=F) %>%
            cbind(., bin_4_coords) %>%
            select(Annotation, start1=start, end1=end, strand1=strand)

bin_5_ann = read.table('bin_5.txt', header=T, sep='\t', stringsAsFactors=F) %>%
            cbind(., bin_5_coords) %>%
            select(Annotation, start2=start, end2=end, strand2=strand)
```

Here, both resulting data.frames have a common column `Annotation`, and then three data.frame-specific columns. The `start1`, `end1`, `start2`, and `end2` names are used directly by `genoPlotR`. The `strand1` and `strand2` we will use to create a new column.

We now want to merge these two data.frames according to the `Annotation` column. This can be done using an `inner_join` function, which takes two data.frames and combines them according to a shared column. Table join operations are a common function in database and table manipulation languages, and there are [several different types of operation](https://dplyr.tidyverse.org/reference/join.html). We wish to use `iner_join` specifically as it will only return instances where a value was found between both tables, and discard rows from either table that lack a partner.

```R
ann_merged = inner_join(bin_4_ann, bin_5_ann, by='Annotation')
```

We will then create a flag for where genes are in the same (forward) orientation, or are reversed. This is done by multiplying the two `strand` columns together, using the following bit-wise logic:

|Gene 1 (orientation)|Gene 1 (value)|Gene 2 (orientation)|Gene 2 (value)|Product|Interpretation|
|:---:|:---:|:---:|:---:|:---:|:---:|
|Forward|1|Forward|1|1|Genes are in same orientation|
|Forward|1|Reverse|-1|-1|Genes are reversed|
|Reverse|-1|Forward|1|-1|Genes are reversed|
|Reverse|-1|Reverse|-1|1|Genes are in same orientation|

Genes in the same orientation will be matched with gray lines, reversed genes with purple. This is achiveved with the `mutate` function and the `case_when` operator.

```R
ann_compute_orienation = ann_merged %>%
                         mutate( direction=case_when( strand1 != strand2 ~ -1, TRUE ~ 1) ) %>%
                         mutate( col=case_when( direction == 1 ~ 'grey', TRUE ~ 'purple') )
```

Finally, we want to strip away the unwanted columns, and cast the resulting table into a *comparison* object.

```R
bin4_bin5_comparison = ann_compute_orienation %>%
                       select(-strand1, -strand2, -Annotation) %>%
                       as.comparison(.)
```

##### Part 2.3 - Creating an annotation table

Before we create the plot, we need to create one more table with descriptive information. This is an *annotation* object, which contains the name of each gene in the figure along with the coordinates to write the name. We can use the `middle` function to automatically calculate the mid-point of each gene in bin_4 for our plotting positions, and the gene annotations from bin_4 as our text. I an chosing bin_4 here as there is a gene loss in bin_5 meaning it does not contain the full complement of genes in this figure.

```R
mid_pos = middle(bin_4_ds)

bin4_bin5_annot = annotation(x1=mid_pos, text=bin_4_ann$Annotation, rot=45)
```

##### Part 2.4 - Creating the plot

Now we can pass these data into the plotting function, to get our final figure.

```R
plot_gene_map(dna_segs=list(bin_4_ds, bin_5_ds),
              comparisons=list(bin4_bin5_comparison),
              annotations=bin4_bin5_annot )
```

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex14_gene_view.PNG)

All done! We can see here that compared to bin_4 the following differences are apparent in the gene operon for bin_5:

1. There has been a loss of *flgD* and *flgF*
1. There has been a duplication of *flgG*
1. *flgE* has been translocated upstream relative to the other genes
1. *flgF*, *flgH*, and *flgI* have been reversed relative to *flgE*

However, it is important to note that these do not all constitute independent events. Careful analysis would be needed to determine whether this rearrangement occurred in a single evolutionary event, or multiple successive ones.

---
