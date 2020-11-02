# Presentation of data

### Objectives

* Overview, `RStudio`, and using `R` in the `Jupyter` environment
* Data visualisation and accessibility
* Build a heatmap of average coverage per sample using `R`
* Build a sulfur assimilation gene alignment figure to investigate gene synteny using `R`
* Build a KEGG pathway map using `R`
* *Optional*: Build a basic heatmap from `BLAST` data using `R`
* *Optional*: Build nMDS plots using `R` to investigate beta-diversity. 

---

### Overview, *R/RStudio*, and using *R* in the *Jupyter* environment

There are a number of powerful packages within the `R` software environment which can be used to create high quality images of genomic and metagenomic data. While each of these packages comes with its own documentation, these documents and tutorials usually assume that your data is already in some already-prepared format. Our data will almost never be in this format, though, so these exercises show two brief examples of how we can scrape data from our existing files to create useful figures. As such, these examples are more complicated than what you would get reading the tutorials and manuals of the plotting tools, but will be transferable to your own work.

In your own work, it may be preferable to download the relevant files from NeSI (e.g. via `scp ...`) and work with them on a locally installed version of `RStudio` on your own machine. For today, to be able to run these `R` exercises in a stable environment within the NeSI platform, we will be running an `R` [kernel](https://en.wikipedia.org/wiki/Kernel_(operating_system)) from within a [Jupyter](https://jupyter.org/) environment. 

`Jupyter Notebooks` provide an interactive space that allows for mixing multiple languages within a single document, including [Markdown](https://en.wikipedia.org/wiki/Markdown), `Python`, and `R` (by default, `Markdown` and one coding language such as `R` can be used within one document, but there are add-ons available to expand this, should you wish to). `Jupyter Notebooks` can be extremely useful as a workspace that is the equivalent of an electronic "lab book". Today, we will be using it as an interactive space to run `R`. Note that, while the layout will be slightly different to `RStudio`, the commands we will be working through will work the same in both environments.

These exercises will take place with files in the `9.data_presentation/` folder.

---

### Data visualisation and accessibility

In this section, we will work through a number of example exercises for visualising various aspects of metagenomic data. 

As the fundamental point of data visualisation is *communication*, when building figures it is important to be mindful of aspects of your figures that might affect the accessibility of what you're trying to communicate (i.e. to maximise the number of people you will be communicating effectively with). A considerable number of your intended audience will be affected by one of the forms of colour vision deficiency (colour blindness). There are a number of useful resources available online for both selecting and testing the appropriateness of your colour selection. Some include:

* [ColorBrewer2](https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3) (select 'colorblindsafe')
* [chroma.js](https://gka.github.io/palettes/#/7|d|6e5300,7c8c00,00a63e|ffffe0,ff005e,93003a|1|1)
* Selecting and checking your colour choice using [Viz Palette](https://projects.susielu.com/viz-palette?colors=[%22#ffd700%22,%22#ffb14e%22,%22#fa8775%22,%22#ea5f94%22,%22#cd34b5%22,%22#9d02d7%22,%22#0000ff%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22achromatopsia%22)
* [Blog post](https://bconnelly.net/posts/creating_colorblind-friendly_figures/) by Brian Connelly.
* Several useful colour palettes designed by Martin Krzywinski are available [here](http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container)

We have been mindful to make appropriate colour selections throughout these examples, but please do let us know if you spot any we've overlooked.

---

### Getting started: logging in to the NeSI *Jupyter hub*

To get started, log in to NeSI's `Jupyter hub` space via a standard browser at [jupyter.nesi.org.nz/hub/login](https://jupyter.nesi.org.nz/hub/login). You will need to enter your NeSI user name, NeSI password, and two-factor authentication.

Select `Genomics Aotearoa Virtual Lab Training Access (nesi02659)` as the project, and select walltime of 8 hours. For now, leave the rest of the options at the default settings.

Within the `Jupyter` launcher, open a `Notebook` running the `R 4.0.1` module as the kernel. 

All of the required packages for these exercises are already installed within the `R 4.0.1 (gimkl-2020a)` kernel we are running. If you need to install these on your local `R` or `RStudio`, this can be done via the `install.packages()` command within `R` or `RStudio`:

```R
install.packages('ade4')
install.packages('genoPlotR')
install.packages('pheatmap')
install.packages('gplots')
install.packages('vegan')

# Install 'tidyverse' (includes: ggplot2, tibble, tidyr, readr, purrr, dplyr, string4, forcats)
install.packages('tidyverse')

# Install 'pathview' package (part of Bioconductor)
if (!require(BiocManager)) {
  install.packages("BiocManager")
  BiocManager::install("pathview", update = FALSE)
}

# Install 'KEGGREST' package (part of Bioconductor)
if (!require(BiocManager)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST", update = FALSE)
}
```

Spend a few minutes familiarising yourself with the `Jupyter Notebook` workspace. You'll see the coding language kernel running in the background of this `Notebook` in the top right of the pane. The main document works in *blocks*; click the `+` button to add additional blocks. The *Code* drop-down menu allows you to select between whether the current block is a `Markdown` or `code` (in this case, `R`) block. 

To start, create a title and description for this Notebook. 

* Create a title
  * Click on the first code block
  * use the drop-down menu to convert the block to `Markdown`
  * enter a title in the block preceeded by one or more `#` symbols (e.g. `# MGSS2020: Data presentation`). In `Markdown`, the `#` symbol denotes this to be rendered as a title.
  * Click on the block and then press `<shift> + <enter>` to 'run' the block (in this case, to render the `Markdown` text).
* Create a comment or sub-title
  * Create a second block
  * convert to `Markdown` and enter a short description of the workshop (e.g. `Data presentation exercises for the Genomics Aotearoa Metagenomics Workshop, 2020`). 
  * as above, press <shift> + <enter> to render the text
   
*NOTE: for `Markdown` blocks, simply double click on a block to be able to re-edit the rendered section.* 

*NOTE: To run a `code` block, use the same method of clicking on the block and pressing `<shift> + <enter>` to interactively run the code for that block (in this case, to run our `R` code).*

---

### Build a heatmap of average coverage per sample using *R*

One of the first questions we often ask when studying the ecology of a system is: What are the pattens of abundance and distribution of taxa across the different samples? In the previous [coverage and taxonomy](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md) exercises we generated per-sample coverage tables by mapping the quality-filtered unassembled reads back to the refined bins and the viral contigs to then generate coverage profiles for each. 

As a reminder:

> Genomes in higher abundance in a sample will contribute more genomic sequence to the metagenome, and so the average depth of sequencing coverage for each of the different genomes provides a proxy for abundance in each sample.

A simple way to present this information is via a heatmap. In this exercise we will build a clustered heatmap of these coverage profiles in `R`. Since we also have tables of taxonomy assignments (via `gtdb-tk` for MAGs) and/or predictions (via `vContact2` for viral contigs), we will also use these to add taxonomy information to the plot.

The coverage and taxonomy tables generateed in earlier exercises have been copied to `9.data_presentation/coverage/` a for use in these exercises.

In addition to this, a simple mapping file has also been created (within `9.data_presentation/coverage/mapping_file.txt`). This is a tab-delimited file listing each sample ID in the first column, and the sample "Group" in the second column (*Group_A*, *Group_B*, and *Group_C*). This grouping might represent, for example, three different sampling sites that we want to compare between. If you had other data (such as environmental measures, community diversity measures, etc.) that you wish to incorporate in other downstream analyses (such an fitting environmental variables to an ordination, or correlation analyses) you could also add new columns to this file and load them at the same time.

*NOTE: It is usually necessary to normalise coverage values based on sample depth. For example, by normalising all coverages for each sample based on either minimum or average library size (the number of sequencing reads per sample). In this case, the mock metagenome data we have been working with were already of equal depth, and so we will omit any normalisation step here. But this would not be the case with your own sequencing data sets.*

### Part 1 - Building a heatmap of MAG coverage per sample.

#### 1.1 Set working directory, load *R* libraries, and import data

First, set the working directory and load the required libraries.

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/9.data_presentation/')

# tidyverse libraries 
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Other libraries
library(gplots)
library(vegan)
```

*NOTE: after copying this code into a code block in `Jupyter`, remember that, to run the code, press `<shift> + <enter>` with the code block selected.*

Import the coverage and taxonomy tables. When importing the files, we will select only the information of use here and will do a few basic first data wrangling steps. For the coverage table, we will select the `contigName` column and each of the columns of coverage values (columns `sample[1-4].bam`). For taxonomy, we will select the `user_genome` column (converted to `Bin` ID) and `classification` (coverated to `taxonomy`), and will also use `gsub` to extract just the taxonomic ranks of interest (in this case, we will extract the *class* to colour code MAGs in the heat map, and *genus* to add to MAG names in the plot), and to add `Unassigned` to any MAGs not assigned at these ranks. (*NOTE: to view a different taxonomic rank, you will need to change the two `mutate(taxonomy_class = gsub(...))` rows below accordingly*).

```R
# coverage table
cov_MAG <- read_tsv("coverage/bins_cov_table.txt") %>% 
  select(c('contigName', ends_with('.bam'))) 

# taxonomy table
tax_MAG <- read_tsv("coverage/gtdbtk.bac120.summary.tsv") %>% 
  mutate(Bin = gsub("(.*).filtered", "\\1", .$user_genome)) %>% 
  mutate(taxonomy = gsub(".*;c(.*);o.*", "class\\1", classification)) %>% 
  mutate(taxonomy = gsub("^class__$", "Unassigned", taxonomy)) %>% 
  mutate(taxonomy_genus = gsub(".*;g__(.*);.*", "\\1", classification)) %>% 
  mutate(taxonomy_genus = gsub("^$", "Unassigned", taxonomy_genus)) %>% 
  select(c('Bin', 'taxonomy', 'taxonomy_genus'))

# mapping file (import both columns as factors: col_types set to factor, factor)
map.df <- read_tsv("coverage/mapping_file.txt", col_types = "ff")
```

#### 1.2 wrangle data

As noted during the [coverage and taxonomy](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md) exercises, it is important to remember that we currently have a table of coverage values for all *contigs* contained within each MAG. Since we're aiming to present coverage for each *MAG*, we need to reduced these contig coverages into a single mean coverage value per MAG per sample. 

In the following code, we first strip the `.bam` extensions off of our sample names. We will then leverage the fact that we added bin IDs to each of the contig headers earlier to re-extract the bin ID for each using `gsub`, use the `group_by()` function to group by `Bin`, and the `summarise()` function to return the per-sample mean coverage of each set of contigs contained within each bin. 

```R
# Extract BinID, group by Bin, calculate mean coverages for each set of contigs per Bin
cov_MAG <- cov_MAG %>%
  rename_at(
    vars(contains("sample")),
    list(~ str_replace(., ".bam", ""))
  ) %>%
  mutate(Bin = gsub("(.*)_NODE.*", "\\1", .$contigName)) %>% 
  group_by(Bin) %>% 
  summarise(across(where(is.numeric), mean))
```

Finally, collate the coverage and taxonomy tables into a single table for downstream use, and rename bins to include the genus name.

```R
# Collate coverage and taxonomy
collate_data_MAG <- cov_MAG %>% 
  left_join(tax_MAG) %>% 
  mutate(Bin = paste(Bin, taxonomy_genus, sep="_")) %>%
  select(-taxonomy_genus) %>%
  mutate(taxonomy = replace_na(taxonomy, "Unassigned"))
```

In many real-life cases, to enhance the visualisation across a wide range of coverage values, you may wish to perform a transformation on your data. 

Perform a log(2)-transformation on the coverage data (those values in the columns that `matches("sample[1-9]")`). Note, here we are calculating log2(x + 1) to allow for any cases where coverage values are 0 in any of the samples (log2(0 + 1) = 0). 

```
# Log(2)-transform coverage data
collate_data_MAG_log2 <- collate_data_MAG
collate_data_MAG_log2[,names(select(collate_data_MAG_log2, contains("sample")))] <- log2(collate_data_MAG_log2[,names(select(collate_data_MAG_log2, contains("sample")))] + 1)
```

To have the data in a format that is ready for `heatmap2` we will perform a few final data wrangling steps. Here we are re-ordering the file by row sums, selecting only the coverage columns (`contains("sample")`) and `taxonomy` column, ensuring `taxonomy` is set as a factor, and setting the rownames based on the `Bin` column.

```R
coverage.heatmap.data.MAG <- collate_data_MAG_log2 %>%
  mutate(sum = rowSums(select(collate_data_MAG_log2, contains("sample")))) %>% 
  arrange(desc(sum)) %>%
  select("Bin", contains("sample"), "taxonomy") %>%
  mutate(taxonomy = factor(taxonomy)) %>% 
  droplevels() %>% 
  column_to_rownames(var = "Bin") %>%
  as.data.frame() 
```

As there are only 10 bins in our mock metagenome data, spread over 7 different *classes* (the taxonomic rank we've chosen), we can include colour codes for all classes here. If you had many more cases of the selected taxonomic rank, it would be necessary to select a subset to colour code (e.g. the top 10 taxa, and grouping all others into `Other`) so as not to have so many colours that things become unintelligible.

**We *don't* need to run this code block today** (in fact, it will return an error due to having too few bins compared to the number we're trying to subset). But if you wished to do this for your own work, you could run something like the following to identify the 10 most abundant taxa (at the selected rank), and change the taxonomy of all others to `Other`. 

> ```R
> ## Identify most abundant Phyla in MAG data
> # Aggregate rows based on taxonomic assignments, reorder by overall relative abundance
> cov_MAG.tax <- collate_data_MAG_log2 %>% 
>   replace_na(list(taxonomy = "Unassigned")) %>% 
>   mutate(sum = rowSums(select(collate_data_MAG_log2, contains("sample")))) %>% 
>   group_by(taxonomy) %>% 
>   summarise_if(is.numeric, list(sum = sum)) %>% 
>   arrange(desc(sum_sum)) %>% 
>   mutate(tax_summary = factor(c(taxonomy[1:10], rep("Other", length(taxonomy)-10))))
> 
> ## Add summaries of top 14 taxa to coverage.heatmap.data.MAG
> coverage.heatmap.data.MAG <- coverage.heatmap.data.MAG %>% 
>   rownames_to_column("Bin") %>% 
>   left_join(select(cov_MAG.tax, c("taxonomy", "tax_summary"))) %>% 
>   mutate(taxonomy = factor(tax_summary)) %>% 
>   column_to_rownames("Bin")
> ```

#### 1.3 Calculate hierarchical clustering of columns (samples) and rows (MAGs).

It can be useful to arrange our rows and/or columns by some form of clustering. The plotting function `heatmap2` can do this for us. However, here we will first perform this ourselves to be able to view the clustering output separately. We can then pass the same clustering to `heatmap2` using the `as.dendrogram()` function within the `heatmap2` command.

Some data sets will encounter an issue with the clustering calculation due to some variables having insufficient variance. Let's first perform a filter to remove any potential problem rows from the data.

```R
# Subset the relevant columns
coverage.heatmap.data.MAG.filt <- select(coverage.heatmap.data.MAG, contains("sample"), "taxonomy")

# Filter out zero-variance rows
coverage.heatmap.data.MAG.filt <- coverage.heatmap.data.MAG.filt[!apply(select(coverage.heatmap.data.MAG.filt, -"taxonomy")==select(coverage.heatmap.data.MAG.filt, -"taxonomy")[[1]],1,all),]
```

Run the hierarchical clustering calculations based on [Bray-Curtis dissimilarity](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity) by column and row.

```R
cov_clus.avg.col <- hclust(vegdist(t(select(coverage.heatmap.data.MAG.filt, contains("sample"))), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE), "aver")
cov_clus.avg.row <- hclust(vegdist(select(coverage.heatmap.data.MAG.filt, contains("sample")), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE), "aver")
```

Let's have a quick look at the outputs of each using the basic `plot` function.

```R
plot(cov_clus.avg.col, hang = -1, cex = 1.5)
plot(cov_clus.avg.row, hang = -1, cex = 1.5)
```

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_MAGs_BC_hclust_Samples.png)

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_MAGs_BC_hclust_MAGs.png)

If you wish to write these to file, we can wrap them in the `png(...)` and `dev.off()` lines, as below (this is true of all of the figures we will be generating):

```R
png("coverage/MAGs_BC_hclust_Samples.png", width=17, height=10, units="cm", res=300)
plot(cov_clus.avg.col, hang = -1, cex = 0.5)
dev.off()

png("coverage/MAGs_BC_hclust_MAGs.png", width=17, height=10, units="cm", res=300)
plot(cov_clus.avg.row, hang = -1, cex = 0.5)
dev.off()
```

#### 1.4 Set the colour palettes

Before generating the heat map, let's set some colour palettes to use within the plot. We will create a greyscale for the coverage values in the heat map, one colour palette to colour the rows by the taxonomy of the MAG, and an additional palette for the columns (sample groups). 

*NOTE: The list of colours for the taxonomy palette are taken from a 15-colour colour blind-friendly palette available [here](http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container). The `Group.col` palette was selected using [Chroma.js](https://gka.github.io/palettes/#/7|d|6e5300,7c8c00,00a63e|ffffe0,ff005e,93003a|1|1).*

```R
# Load greyscale palette
scalegreys <- colorRampPalette(c("white","black"), space="Lab")(100)

# Taxonomy colour palette
MAG.cols <- c(
  "#006DDB","#FF6DB6","#DBD100","#FFB677","#004949","#009292","#FFFF6D","#924900","#490092","#24FF24","#920000","#B6DBFF","#B66DFF","#6DB6FF","#000000"
)

# Data.frame containing 'sample groups' colour palette
Group.col <- data.frame(
  Group = factor(levels(map.df$Group)),
  col = c('#71dfdf', '#3690b8', '#00429d'),
  stringsAsFactors = FALSE)
```

Finally, join the `Group.col` colours to the mapping file (`map.df`) so that we can call them in the `heatmap.2` command.

```R
# Update map.df with `Group` colour code for each sample
map.col <- data.frame(full_join(map.df,Group.col))
```

#### 1.5 Build the heat map

Output the full heat map with sample (column) and variable (row) clustering.

*NOTE: Un-comment (remove the `#` symbol from) the lines `png(...)` and `dev.off()` before and after the heat map code will write the plot to file*

A number of options for visual parameters within `heatmap2` can be set. Of the ones provided below, the first six (from `Colv` to `ColSideColors`) are likely to be the most of interest at this stage. `Colv` and `Rowv` set the clustering (using the Bray-Curtis hierarcical clustering we calculated above. The `cex...` commands set the font sizes. The `...SideColors` commands set the colours placed adjacent to the rows or columns to add further information to the plot (in this case, we are using this to add MAG taxonomy information and sample grouping information). 

In some cases, you may wish to omit the column clustering (for example, if your samples are from a set transect and you want to preserve the order of samples in the plot). In this case, you can set `Colv = FALSE`.

*NOTE: Here we are using gplot::legend to add two legends to the bottom of the plot for row-side and column-side colours. This can sometimes look a little clunky (and often requires experimenting with the `margins` settings). For final publication-quality figures, it may be preferable to output each element (heatmap, legend1, and legend2) as separate image files, and then compile in a program such as Adobe Illustrator.*

```R
#png("coverage/MAGs_heatmap.png",width=17,height=21,units="cm",res=300)
heatmap.2(
  as.matrix(select(coverage.heatmap.data.MAG.filt, contains("sample"))),
  Colv = as.dendrogram(cov_clus.avg.col),
  Rowv = as.dendrogram(cov_clus.avg.row),
  cexCol = 1.2,
  cexRow = 1.2,
  RowSideColors = MAG.cols[coverage.heatmap.data.MAG.filt[,"taxonomy"]],
  ColSideColors = map.df$col,
  margins = c(28,16),
  lwid = c(1, 8),
  lhei = c(1, 8),
  col = scalegreys,
  na.color = "white",
  scale = 'none',
  srtCol = 45,
  density.info="none",
  xlab = NA,
  ylab = NA,
  trace = c('none'),
  offsetRow = 0.2, 
  offsetCol = 0.2, 
  key.title = NA,
  key.xlab = "log2(coverage)",
  key.par=list(mgp=c(1.5, 0.5, 0),
               mar=c(3, 1, 3, 0),
               cex = 0.5),
)
legend("bottomleft",
       legend = levels(coverage.heatmap.data.MAG.filt[,"taxonomy"]),
       border = NULL,
       col = MAG.cols[1:length(levels(coverage.heatmap.data.MAG.filt[,"taxonomy"]))],
       lty = 1,
       lwd = 6,
       bty = "n",
       ncol = 1,
       cex = 0.8,
       title = "Y-axis key:\nMAG taxonomy"
)
legend("bottomright",
       legend = levels(map.df$Group),
       border = NULL,
       col = unique(map.df$col),
       lty = 1,
       lwd = 6,
       bty = "n",
       ncol = 1,
       cex = 0.8,
       title = "X-axis key:\nSample groups"
)
#dev.off()
```

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_MAGs_heatmap.png)

Some quick take-aways from looking at this plot:

* We can see that three MAGs are from the same class (*Gammaproteobacteria*)
* Sample 4 (from sample group *Group_C*) is somewhat of an outlier, with only two bins (MAGs) recovered. Furthermore *bin_7_Nitrobacter* was *only* recovered from this sample.
* Samples from sample groups A and B each recovered the same sets of bins. 
* With more samples, we would also be able to better discern co-occurence patterns of the MAGs across the samples by examining the clustering patterns in the dendrogram on the left.

Feel free to experiment with some of the settings above. For example:

* What happens if you set `Colv = TRUE, Rowv = TRUE` instead of providing our pre-calculated Bray-Curtis dissimilarity-based clustering? Why might this cluster differently? (Hint: have a look at the documentation for `heatmap.2` to see it's default settings for clustering: `?heatmap.2()`). 
  * *NOTE: in the `Jupyter` environment, leaving a manual for a function open can be cumbersome. To clear the manual, right-click on the code block containing `?heatmap.2()` and click `Clear Outputs`.*
* Try removing the column clustering (but retaining the row clustering) by setting `Colv = FALSE, Rowv = as.dendrogram(cov_clus.avg.row)`.
* Play around with the `margins = ...` settings to see what difference each makes.


### Part 2 - Building a heatmap of viral contigs per sample.

We can run through the same process for the viral contigs. Many of the steps are as outlined above, so we will work through these a bit quicker and with less commentary along the way. However, we will highlight a handful of differences compared to the commands for the MAGs above, for example:  steps for selecting and/or formatting the taxonomy; importing the quality output from `CheckV`; and the (optional) addition of filtering out low quality contigs.

#### 2.1 Set working directory, load *R* libraries, and import data

In this case, import the coverage, taxonomy, mapping, *and checkv* files. For taxonomy we will select the `Genome` column (converted to `Contig`), `Order_VC_predicted` (converted to `taxonomy`), and `Genus_VC_predicted` (converted to `taxonomy_genus`) (we will use the `Order_VC_predicted` to colour code viral contigs in the heat map, and `Genus_VC_predicted` to add to viral contig names in the plot). 

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/9.data_presentation/')

# tidyverse libraries 
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Other libraries
library(gplots)
library(vegan)

# coverage table
cov_vir <- read_tsv("coverage/viruses_cov_table.txt") %>% 
  mutate(Contig = contigName) %>%
  select(c('Contig', ends_with('.bam'))) 

# taxonomy table
tax_vir <- read_tsv("coverage/vcontact2_tax_predict_table.txt") %>% 
  mutate(Contig = Genome) %>% 
  mutate(taxonomy = factor(Order_VC_predicted)) %>% 
  mutate(taxonomy_genus = factor(Genus_VC_predicted)) %>% 
  select(c('Contig', 'taxonomy', 'taxonomy_genus'))

# CheckV quality
checkv <- read_tsv("coverage/checkv_quality_summary.tsv") %>%
  mutate(Contig = contig_id) %>%
  select(c('Contig', 'checkv_quality')) %>%
  mutate(checkv_quality = factor(checkv_quality, levels = c("Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"))) 

# mapping file (import both columns as factors: col_types set to factor, factor)
map.df <- read_tsv("coverage/mapping_file.txt", col_types = "ff")
```

#### 2.2 wrangle data

This time we will collate the coverage, checkv, and taxonomy tables into a single table for downstream use, and create `vir_ID_*` unique identifers.

```R
# Strip .bam from sample names
cov_vir <- cov_vir %>%
  rename_at(
    vars(contains("sample")),
    list(~ str_replace(., ".bam", ""))
  )

# Collate data, add unique 'vir_ID' identifer, add 'vir_ID_tax' column
collate_data_vir <- cov_vir %>% 
  left_join(tax_vir, by = 'Contig') %>% 
  left_join(checkv, by = 'Contig') %>% 
  mutate(vir_ID = paste("vir", seq(1, length(Contig)), sep="_")) %>% 
  relocate('vir_ID') %>% 
  mutate(vir_ID_tax = paste(vir_ID, taxonomy_genus, sep="_")) %>% 
  filter(!is.na(checkv_quality)) %>%
  select(-taxonomy_genus)

# Log(2)-transform coverage data
collate_data_vir_log2 <- collate_data_vir
collate_data_vir_log2[,names(select(collate_data_vir_log2, contains("sample")))] <- log2(collate_data_vir_log2[,names(select(collate_data_vir_log2, contains("sample")))] + 1)
```

Now let's add an (optional) filtering step here to remove any contigs that `CheckV` returned a quality category of `Not-determined` or `Low-quality`. We will use `!=` in the `filter()` function here to return only those rows that do *not* include 'Not-determined' or 'Low-quality').

```R
collate_data_vir_log2 <- collate_data_vir_log2 %>%
  filter(checkv_quality != 'Not-determined' & checkv_quality != 'Low-quality')
```

A few final data wrangling steps:

```R
coverage.heatmap.data.vir <- collate_data_vir_log2 %>%
  mutate(sum = rowSums(select(collate_data_vir_log2, contains("sample")))) %>% 
  arrange(desc(sum)) %>%
  select("vir_ID", contains("sample"), "taxonomy") %>%
  mutate(taxonomy = factor(taxonomy)) %>% 
  droplevels() %>% 
  column_to_rownames(var = "vir_ID") %>%
  as.data.frame() 
```

*NOTE: In this case, the script only keeps the `vir_ID*` identifiers without the genus taxonomy appended for use in the plots, since the the long labels (due to numerous genus predictions for some contigs) make the plots illegible otherwise. To retain the taxonomy, change the two `vir_ID` entries above to `vir_ID_tax`.*

#### 2.3 Calculate hierarchical clustering of columns (samples) and rows (viral contigs).

Filter to remove any potential problem rows from the data, and then run clustering and plot the dendrograms.

```R
# Subset the relevant columns
coverage.heatmap.data.vir.filt <- select(coverage.heatmap.data.vir, contains("sample"), "taxonomy")

# Filter out zero-variance rows
coverage.heatmap.data.vir.filt <- coverage.heatmap.data.vir.filt[!apply(select(coverage.heatmap.data.vir.filt, -"taxonomy")==select(coverage.heatmap.data.vir.filt, -"taxonomy")[[1]],1,all),]

# hclust
vir_cov_clus.avg.col <- hclust(vegdist(t(select(coverage.heatmap.data.vir.filt, contains("sample"))), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE), "aver")
vir_cov_clus.avg.row <- hclust(vegdist(select(coverage.heatmap.data.vir.filt, contains("sample")), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE), "aver")

# plot
#png("coverage/Viruses_BC_hclust_Samples.png", width=17, height=10, units="cm", res=300)
plot(vir_cov_clus.avg.col, hang = -1, cex = 1.5)
#dev.off()

#png("coverage/Viruses_BC_hclust_MAGs.png", width=17, height=10, units="cm", res=300)
plot(vir_cov_clus.avg.row, hang = -1, cex = 1.5)
#dev.off()
```

#### 2.4 Set the colour palettes

```R
# Load greyscale palette
scalegreys <- colorRampPalette(c("white","black"), space="Lab")(100)

# Taxonomy colour palette
vir.cols <- c(
  "#006DDB","#FF6DB6","#DBD100","#FFB677","#004949","#009292","#FFFF6D","#924900","#490092","#24FF24","#920000","#B6DBFF","#B66DFF","#6DB6FF","#000000"
)

# Data.frame containing 'sample groups' colour palette
Group.col <- data.frame(
  Group = factor(levels(map.df$Group)),
  col = c('#71dfdf', '#3690b8', '#00429d'),
  stringsAsFactors = FALSE)

# Update map.df with `Group` colour code for each sample
map.df <- data.frame(full_join(map.df,Group.col))
```

#### 2.5 Build the heat map

Output the full heat map with sample (column) and variable (row) clustering.

*NOTE: Un-comment (remove the `#` symbol from) the lines `png(...)` and `dev.off()` before and after the heat map code will write the plot to file*

```R
#png("coverage/Viruses_heatmap.png",width=17,height=21,units="cm",res=300)
heatmap.2(
  as.matrix(select(coverage.heatmap.data.vir.filt, contains("sample"))),
  Colv = as.dendrogram(vir_cov_clus.avg.col),
  Rowv = as.dendrogram(vir_cov_clus.avg.row),
  cexCol = 1.2,
  cexRow = 1.2,
  RowSideColors = MAG.cols[coverage.heatmap.data.vir.filt[,"taxonomy"]],
  ColSideColors = map.df$col,
  margins = c(28,16),
  lwid = c(1, 8),
  lhei = c(1, 8),
  col = scalegreys,
  na.color = "white",
  scale = 'none',
  srtCol = 45,
  density.info="none",
  xlab = NA,
  ylab = NA,
  trace = c('none'),
  offsetRow = 0.2, 
  offsetCol = 0.2, 
  key.title = NA,
  key.xlab = "log2(coverage)",
  key.par=list(mgp=c(1.5, 0.5, 0),
               mar=c(3, 1, 3, 0),
               cex = 0.5),
)
legend("bottomleft",
       legend = levels(coverage.heatmap.data.vir.filt[,"taxonomy"]),
       border = NULL,
       col = MAG.cols[1:length(levels(coverage.heatmap.data.vir.filt[,"taxonomy"]))],
       lty = 1,
       lwd = 6,
       bty = "n",
       ncol = 1,
       cex = 0.8,
       title = "Y-axis key:\nVirus taxonomy (Order)"
)
legend("bottomright",
       legend = levels(map.df$Group),
       border = NULL,
       col = unique(map.df$col),
       lty = 1,
       lwd = 6,
       bty = "n",
       ncol = 1,
       cex = 0.8,
       title = "X-axis key:\nSample groups"
)
#dev.off()
```

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_Viruses_heatmap.png)

---

### Build a sulfur assimilation gene alignment figure to investigate gene synteny using `R`

When investigating the evolution of genomes, we sometimes want to consider not only the presence/absence of genes in a genome, but also how they are arranged in an operon. For this exercise, we are going to visualise several sulfur assimilation genes from bin_2, bin_3 and bin_4, comparing their arrangement among the organisms.

For this exercise, navigate to the folder `9.data_presentation/gene_synteny/`. Here, you have been provided with a copy of the `prodigal` gene predictions for each of the bins (`.faa` files), an annotation output table using multiple databases (`.aa` files), a small table of the annotation of some key genes of interest (`cys.txt` files), and blastn output (`blast*.txt`) comparing the genes of interest from these organisms. The annotation files were created by manually searching the annotations obtained in the previous exercises.

*NOTE: Refer to [gene_synteny_grab_GOI.md](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/gene_synteny_grab_GOI.md) for more information on how the `cys.txt` files were produced.*

#### Part 1 - Parsing files in bash

We will be performing this exercise in two stages. Firstly, in `bash`, we will use `cut` and `tail` to pull out the genes of interest listed in the `*cys.txt` files from the `prodigal` files. The gene names will then be used to create a table of gene coordinates from the `prodigal` output using `grep`, `cut`, and `sed`.

For these `bash` steps, we will need to return to our logged in NeSI terminal. You can use the same terminal you have been using for the rest of the workshop if you wish. But fortunately, the NeSI `Jupyter hub` includes a terminal window as well. 

Launch a standard **terminal** within the NeSI [Jupyter hub](https://jupyter.nesi.org.nz/hub/login):

* click on the `+` button in the top left of the screen (directly below the `File` drop-down menu) to bring up the `Jupyter` launcher window
* launch the terminal using the bottom at the bottom left of the pane (under "Other"). 

From here, we can run the `bash` commands below.

First, navigate to the `9.data_presentation/gene_synteny/` folder, and then perform the `cut` and `tail` steps outlined above.

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/9.data_presentation/gene_synteny/

cut -f3 bin_2_cys.txt | tail -n+2 > bin_2_cys.genes
cut -f3 bin_3_cys.txt | tail -n+2  > bin_3_cys.genes
cut -f3 bin_4_cys.txt | tail -n+2  > bin_4_cys.genes
```

We now have three new files, ending with the `.genes` suffix which are simply a list of the genes that we wish to extract from the `prodigal` files. We can then use each of these files as input in for `grep` to pull out the *fastA* entries that correspond to these genes.

```bash
grep -f bin_2_cys.genes bin_2.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_2_cys.coords
grep -f bin_3_cys.genes bin_3.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_3_cys.coords
grep -f bin_4_cys.genes bin_4.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_4_cys.coords
```

As previously, we will quickly go through the steps of this command:

```bash
grep -f bin_2_cys.genes bin_2_cys.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_2_cys.coords
|                                     |                      |              |                |
|                                     |                      |              |                Write the output
|                                     |                      |              |
|                                     |                      |              Replace each '#' character with tab
|                                     |                      |
|                                     |                      For each line, replace the '>' with empty text
|                                     |
|                                     Split the results into columns delimited by the # character,
|                                     then take columns 1 - 4.
| Select the lines of the bin_2_cys.faa file that contain entries found in the bin_2_cys.genes file.
```

Check the content of the `.coords` files now. You should see something like the following:

```bash
cat bin_2_cys.coords
# caef037ad318582c885d129aee6008b1_3761    3966783         3967781         1
# caef037ad318582c885d129aee6008b1_3762    3967943         3968761         1
# caef037ad318582c885d129aee6008b1_3763    3968772         3969641         1
# caef037ad318582c885d129aee6008b1_3764    3969645         3970634         1
```

If you recall from the previous exercise on gene prediction, we have taken the first four entries from each line of the `prodigal` output, which consists of:

1. The gene name, written as [CONTIG]\_[GENE]
1. The start position of the gene
1. The stop position of the gene
1. The orienation of the gene

We will now use these tables, together with the annotation tables to create the gene synteny view in `R`.

#### Part 2 - Producing the figure in *R*

First, move back to the `Jupyter Notebook` pane where we have `R` running as the kernel (or relaunch a new `R 4.0.1` Notebook).

There are two `R` libaries we need to load for this exercise,

##### Set working directory and load *R* libraries

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/9.data_presentation/gene_synteny/')

library(dplyr)
library(genoPlotR)
```
    
    Attaching package: ‘dplyr’
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    Loading required package: ade4
    Loading required package: grid


##### Part 2.1 - Load coordinate files

We can now begin importing the data. First, we will import the `.coords` files, and set column names to the files.


```R
bin_2_coords = read.table('bin_2_cys.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_2_coords) = c('name', 'start', 'end', 'strand')
bin_3_coords = read.table('bin_3_cys.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_3_coords) = c('name', 'start', 'end', 'strand')
bin_4_coords = read.table('bin_4_cys.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_4_coords) = c('name', 'start', 'end', 'strand')
```

Take a look at the content of each of these data.frames, by entering their names into the terminal. You should notice that the coordinates occur at quite different positions between the genomes. If we were looking at complete genomes, this would indicate their position relative to the *origin of replication*, but as these are unclosed genomes obtained from MAGs, they reflect the coordinates upon their particular *contig*.

We now parse these data.frames into the *dna_seg* data class, which is defined by the `genoPlotR` library.


```R
bin_2_ds = dna_seg(bin_2_coords)
bin_3_ds = dna_seg(bin_3_coords)
bin_4_ds = dna_seg(bin_4_coords)
bin_4_ds
```


<table>
<caption>A dna_seg: 7 × 11</caption>
<thead>
	<tr><th scope=col>name</th><th scope=col>start</th><th scope=col>end</th><th scope=col>strand</th><th scope=col>col</th><th scope=col>fill</th><th scope=col>lty</th><th scope=col>lwd</th><th scope=col>pch</th><th scope=col>cex</th><th scope=col>gene_type</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_588 </td><td>625636</td><td>626724</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_589 </td><td>626721</td><td>627599</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_590 </td><td>627609</td><td>628442</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_591 </td><td>628508</td><td>629296</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_592 </td><td>629387</td><td>630298</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_593 </td><td>630410</td><td>630715</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>4f94130388f7df8dfba7383a12b9de4d_594 </td><td>630774</td><td>631781</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
</tbody>
</table>


By looking at the table, we can see that the genes in bin_4 are in reversed order (strand = -1), so we might want to change the color of the gene to make sure they are different than bin_2 and bin_3.


```R
bin_4_ds$col = "#1a535c"
bin_4_ds$fill = "#1a535c"
```

##### Part 2.3 - Load annotation tables

Then, we can load the annotation tables we have into `R` and take a look at them.


```R
bin_2_ann = read.table('bin_2_cys.txt', header=T, sep='\t', stringsAsFactors=F) %>%
cbind(., bin_2_coords) %>%
select(Annotation, start1=start, end1=end, strand1=strand)
bin_3_ann = read.table('bin_3_cys.txt', header=T, sep='\t', stringsAsFactors=F) %>%
cbind(., bin_3_coords) %>%
select(Annotation, start1=start, end1=end, strand1=strand)
bin_4_ann = read.table('bin_4_cys.txt', header=T, sep='\t', stringsAsFactors=F) %>%
cbind(., bin_4_coords) %>%
select(Annotation, start1=start, end1=end, strand1=strand)
bin_2_ann
bin_4_ann
```


<table>
<caption>A data.frame: 4 × 4</caption>
<thead>
	<tr><th scope=col>Annotation</th><th scope=col>start1</th><th scope=col>end1</th><th scope=col>strand1</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>sbp; sulfate-binding protei                                       </td><td>3966783</td><td>3967781</td><td>1</td></tr>
	<tr><td>cysT; sulfate transporter CysT                                    </td><td>3967943</td><td>3968761</td><td>1</td></tr>
	<tr><td>cysW; sulfate transporter CysW                                    </td><td>3968772</td><td>3969641</td><td>1</td></tr>
	<tr><td>cysA; sulfate.thiosulfate ABC transporter ATP-binding protein CysA</td><td>3969645</td><td>3970634</td><td>1</td></tr>
</tbody>
</table>




<table>
<caption>A data.frame: 7 × 4</caption>
<thead>
	<tr><th scope=col>Annotation</th><th scope=col>start1</th><th scope=col>end1</th><th scope=col>strand1</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>cysA; sulfate transport ATP-binding ABC transporter protein</td><td>625636</td><td>626724</td><td>-1</td></tr>
	<tr><td>cysW; sulfate transport ABC transporter protein            </td><td>626721</td><td>627599</td><td>-1</td></tr>
	<tr><td>cysU; sulfate transport ABC transporter protein            </td><td>627609</td><td>628442</td><td>-1</td></tr>
	<tr><td>Domain of unknown function 2                               </td><td>628508</td><td>629296</td><td>-1</td></tr>
	<tr><td>Domain of unknown function 2                               </td><td>629387</td><td>630298</td><td>-1</td></tr>
	<tr><td>possible predicted diverged CheY-domain                    </td><td>630410</td><td>630715</td><td>-1</td></tr>
	<tr><td> sbp1; Prokaryotic sulfate-/thiosulfate-binding protein    </td><td>630774</td><td>631781</td><td>-1</td></tr>
</tbody>
</table>



We need to create one more table with descriptive information. This is an *annotation* object, which contains the name of each gene in the figure along with the coordinates to write the name. x1 = starts of the gene, x2 = end of the gene


```R
annot1 <- annotation(x1 = c(bin_2_ds$start[1], bin_2_ds$start[2],bin_2_ds$start[3],bin_2_ds$start[4]),
                    x2 = c(bin_2_ds$end[1], bin_2_ds$end[2], bin_2_ds$end[3], bin_2_ds$end[4]),
                    text = c("sbp", "cysT","cysW","cysA"),
                    rot = 0, col = "black")
annot2 <- annotation(x1 = c(bin_4_ds$start[1], bin_4_ds$start[2],bin_4_ds$start[3],bin_4_ds$start[4],bin_4_ds$start[7]),
                    x2 = c(bin_4_ds$end[1], bin_4_ds$end[2], bin_4_ds$end[3], bin_4_ds$end[6], bin_4_ds$end[7]),
                    text = c("cysA","cysW","cysU","unknown domain" ,"sbp"),
                    rot = 0, col = "black")
```

##### Part 2.3 - Creating a comparison table

Then, we can parse the blast output as comparison file among the genomes. genoPlotR can read tabular files, either user-generated tab files (read_comparison_from_tab), or from BLAST output (read_comparison_from_blast). To produce files that are readable by genoPlotR, the -m 8 or 9 option should be used in blastall, or -outfmt 6 or 7 with the BLAST+ suite.


```R
blast1 = read_comparison_from_blast("blast_bin2_bin3.txt")
blast2 = read_comparison_from_blast("blast_bin3_bin4.txt")
```

What it does here is to set the line color according to the direction of the gene match and the color gradient is based on the percent identity of the matches. Lighter color indicates weaker match and darker color indicates strong match.

Now we can generate the plot. Running this command in `RStudio` or our `Jupyter Notebook` running the `R` kernel interactively loads the figure. 

*NOTE: the commented-out lines below (the two lines starting with `#`) will not run. Un-comment these if you wish to save the figure to file rather than opening it in the `Jupyter` or `RStudio` viewer. (N.b. alternatives using a similar command are also available other formats, including `tiff` and `png`).

```R
#pdf("200926_genoplot_ori.pdf",colormodel = "cmyk",width = 8,height = 4,paper = 'special')
plot_gene_map(dna_segs = list(bin_2_ds,bin_3_ds,bin_4_ds), 
              gene_type = "arrows", dna_seg_labels = c("bin_2", "bin_3","bin_4"), 
              comparisons = list(blast1,blast2), dna_seg_scale = TRUE, 
              override_color_schemes=TRUE,annotations=list(annot1,NULL,annot2),
              annotation_height=1.7, dna_seg_label_cex = 1,main = "Sulfur assimilation")
#dev.off()
```


![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_gene_synteny_fig1.png)


Careful analysis would be needed to determine whether this is a genuine rearrangement relative to the rest of the genome, as these are draft genomes and contig orientation can either be forward or reverse. In this case, you can see that genes in bin_4 are in reversed order relative to the other bin contigs, hence, we can manually rotate the contig.

##### Rotate the contig and update the annotation

Rotate the orientation of the contig from bin_4:

```R
blast2 = mutate(blast2, direction = 1)

annot3 <- annotation(x1 = c(bin_4_ds$start[1], bin_4_ds$start[2],bin_4_ds$start[3],bin_4_ds$start[4],bin_4_ds$start[7]),
                    x2 = c(bin_4_ds$end[1], bin_4_ds$end[2], bin_4_ds$end[3], bin_4_ds$end[6], bin_4_ds$end[7]),
                    text = c("cysA","cysW","cysU","unknown domain" ,"sbp"),
                    rot = 0, col = "black")
#edit the color
bin_4_ds$col = "blue"
bin_4_ds$fill = "blue"
```

Now we can regenerate the figure:

```R
#AFTER ROTATING
#pdf("200926_genoplot_rotated.pdf",colormodel = "cmyk",width = 8, height = 4, paper = 'special')
plot_gene_map(dna_segs = list(bin_2_ds,bin_3_ds,bin_4_ds), 
              gene_type = "arrows", dna_seg_labels = c("bin_2", "bin_3","bin_4"), 
              comparisons = list(blast1,blast2), dna_seg_scale = TRUE, 
              override_color_schemes=TRUE,annotations=list(annot1,NULL,annot3),
              annotation_height=1.7, dna_seg_label_cex = 1,xlims=list(NULL, NULL, 
              c(Inf, -Inf)),main = "Sulfur assimilation")
#dev.off()
```


![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_gene_synteny_fig2.png)


All done! We can see here that compared to bin_2 and bin_3 the following differences are apparent in the gene operon for bin_4:

1. Three genes with unknown domain/function present in bin_4 in between *cysU* and *SBP*
2. *cysW* and *cysA* are highly conserved among the genomes and have higher similarity compared to *cysU* and *SBP*.

---

### Build a KEGG pathway map using *R*

In this exercise, we will generate KEGG a pathways map from genome annotations to visualize potential pathways present in our assembled, binned genome sequences. 

The key package used here is [pathview](https://academic.oup.com/bioinformatics/article/29/14/1830/232698), available from `Bioconductor` (for installation on your local version of `RStudio`, see the [intro section](#getting-started-logging-in-to-the-nesi-jupyter-hub) above). `pathview` is mainly used in transcriptomic analyses but can be adapted for any downstream analyses that utilise the KEGG database. For this exercise, we are interested in visualising the prevalence of genes that we have annotated in a pathway of interest.

#### Set the working directory and load packages and files into *R*

Set the working directory and load the required packages

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/9.data_presentation/gene_synteny/')

library(pathview)
library(KEGGREST)
library(tidyverse)
```

Load your files into R. Here, we are loading them into a list. Given that there are ten files, it is easier to load, clean, and analyze the data using list methods available in tidyverse.

```R
tmp_files <- dir(pattern = "*.aa")
tmp_names <- str_replace(tmp_files, "(.*).annotation.aa", "\\1")
bin_list <- as.list(tmp_files) %>%
  map(
    function(x) {
      read_delim(file = x, delim = "\t")
    }
  )
names(bin_list) <- tmp_names

rm(list = ls(pattern = "tmp*"))
```

*NOTE: R will print some warning messages about column types. However, we do not need all the columns for this analysis and it is not necessary to reformat them.*

#### Subset the data

For this section, we only require the KEGG orthology IDs (KO numbers). Use the code below to subset your data based on these IDs.

Check data headers to determine which columns to subset:

```R
names(bin_list$bin_0)
```

We can see that the headers related to KEGG annotations are in columns 22 to 27. We will subset the relevant columns and concatenate the rows.

```R
kegg_annotations <- map(
  bin_list,
  function(x) {
    x[, 22:27]
  }
) %>%
  bind_rows(.id = "bin")
```

We also need the KO numbers. Here, we will create a column with that information where available (it is in the column Description_2). Note that some annotations do not have KO numbers attached to them. In these cases, we will filter these data out.

```R
kegg_annotations <- kegg_annotations %>%
  mutate(
    "KO" = str_replace(Description_2, "([^:]+):.*", "\\1")
  )

kegg_filtered_annotations <- kegg_annotations %>%
  filter(
    str_detect(KO, "^K\\d+")
  )
```

Check if each query was only annotated with a single KO number. The length of KO numbers are 6 (K followed by 5 numbers). Multiple annotations are possible and will need to be separated into individual rows. For this exercise, we are only interested in identifying the genomic potential for pathways. You can choose to perform further refinement of the annotations later. Let's check for any cases where the `$KO` variable is longer than six characters long:

```R
length(which(str_length(kegg_filtered_annotations$KO) > 6))
```

We can see that there are a total of 62 genes with multiple annotations. We will split them into their own rows. 

```R
kegg_filtered_annotations <- kegg_filtered_annotations %>%
  separate_rows(KO, sep = ",")
```

#### Reshaping the data for input into pathview

`pathview` needs the data as a numeric matrix with IDs as row names and samples/experiments/bins as column names. Here, we will reshape the data into a matrix of counts per KO number in each bin.

We first create a data frame of counts of KO in each bin.

```R
kegg_tally <- kegg_filtered_annotations %>%
  select(bin, KO) %>%
  group_by(bin, KO) %>%
  tally()
```

Turn the tally into a numeric matrix with KO numbers as row names. Do not worry about NAs, `pathview` can deal with that.

```R
kegg_matrix <- kegg_tally %>%
  pivot_wider(names_from = "bin", values_from = "n") %>%
  column_to_rownames("KO") %>%
  as.matrix()
```

#### Creating pathway map of genes related to nitrogen metabolism

Now we can generate images of the KEGG pathway maps using the matrix we just made. For this section, we will try to find genes related nitrogen metabolism. 

We need to first identify the KEGG pathway ID. This is where `KEGGREST` comes in handy. `KEGGREST` can also help you identify other information stored in the KEGG database. For more information, the `KEGGREST` vignette can be viewed using the `vignette` function in `R`:

```R
vignette("KEGGREST-vignette")
```

Use the `keggFind` function to identify the pathway ID for nitrogen metabolism:

```R
keggFind(database = "pathway", query = "nitrogen")
```

The result tells us that the nitrogen metabolism pathway ID is **00910**.

First, let's generate a map using just the *bin_0* data. The output will be a png image file in your working directory.

```R
pv_bin_0 <- pathview(
  gene.data = kegg_matrix[, colnames(kegg_matrix) == "bin_0"], 
  pathway.id = "00910", 
  species = "ko", 
  kegg.native = T, 
  limit = list(
    gene = c(0, max(kegg_matrix[, 1], na.rm = T)),
    cpd = c(0, max(kegg_matrix[, 1], na.rm = T))
  ), 
  bins = list(
    gene = 4,
    cpd = 4
  ),
  both.dirs = list(
    gene = F, 
    cpd = F
  ), 
  out.suffix = "nitrogen_metabolism.bin_0",
  na.col = "white",
  low = "grey",
  mid = "green",
  high = "red"
)
```

You can modify the bins and limit arguments to make it look more sensible for your data. In the meantime, here we've just set it from zero to the maximum of the data.

We can create maps for all bins at the same time using `R`'s implementation of a `for` loop, looping over all of the entries in our `kegg_matrix` `R` object (where each column relates to one of the bins).

```R
pv_all_bins <- list()
for (i in 1:ncol(kegg_matrix)) {
  out.suffix <- paste0("nitrogen_metabolism.", colnames(kegg_matrix)[i])
  pv_all_bins[[colnames(kegg_matrix)[i]]] <- pathview(
    gene.data = kegg_matrix[, i], 
    pathway.id = "00910", 
    species = "ko", 
    kegg.native = T, 
    limit = list(
      gene = c(0, max(kegg_matrix, na.rm = T)),
      cpd = c(0, max(kegg_matrix, na.rm = T))
    ), 
    bins = list(
      gene = 4,
      cpd = 4
    ),
    both.dirs = list(
      gene = F, 
      cpd = F
    ), 
    out.suffix = out.suffix,
    na.col = "white",
    low = "grey",
    mid = "green",
    high = "red"
  )
  rm(i, out.suffix)
}
```

The above code will generate one png image file per bin, and will save it to the current working directory. The list `pv_all_bins` in the first line contains some plotting information. Importantly, it has the KO numbers it used to match the data provided to the pathway. If you wish to, you can subset your original data to only obtain those KO numbers related to the pathway that you have plotted for further analyses with other annotations.

---

### *Optional*: Build a basic heatmap from annotation data using *R*

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
