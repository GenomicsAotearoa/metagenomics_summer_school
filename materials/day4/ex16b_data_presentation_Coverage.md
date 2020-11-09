# Presentation of data: Per-sample coverage heatmaps

### Objectives

* Build heatmaps of average coverage per sample using `R`

---

### Build a heatmap of average coverage per sample using *R*

One of the first questions we often ask when studying the ecology of a system is: What are the pattens of abundance and distribution of taxa across the different samples? In the previous [coverage and taxonomy](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md) exercises we generated per-sample coverage tables by mapping the quality-filtered unassembled reads back to the refined bins and the viral contigs to then generate coverage profiles for each. 

As a reminder:

> Genomes in higher abundance in a sample will contribute more genomic sequence to the metagenome, and so the average depth of sequencing coverage for each of the different genomes provides a proxy for abundance in each sample.

A simple way to present this information is via a heatmap. In this exercise we will build a clustered heatmap of these coverage profiles in `R`. Since we also have tables of taxonomy assignments (via `gtdb-tk` for MAGs) and/or predictions (via `vContact2` for viral contigs), we will also use these to add taxonomy information to the plot.

The coverage and taxonomy tables generated in earlier exercises have been copied to `11.data_presentation/coverage/` a for use in these exercises.

In addition to this, a simple mapping file has also been created (`11.data_presentation/coverage/mapping_file.txt`). This is a tab-delimited file listing each sample ID in the first column, and the sample "Group" in the second column (*Group_A*, *Group_B*, and *Group_C*). This grouping might represent, for example, three different sampling sites that we want to compare between. If you had other data (such as environmental measures, community diversity measures, etc.) that you wish to incorporate in other downstream analyses (such an fitting environmental variables to an ordination, or correlation analyses) you could also add new columns to this file and load them at the same time.

*NOTE: It is usually necessary to normalise coverage values based on sample depth. For example, by normalising all coverages for each sample based on either minimum or average library size (the number of sequencing reads per sample). In this case, the mock metagenome data we have been working with were already of equal depth, and so we will omit any normalisation step here. But this would not be the case with your own sequencing data sets.*

---

### Part 1 - Building a heatmap of MAG coverage per sample.

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login) and open a `Notebook` running the `R 4.0.1` module as the kernel (or, outside the context of this workshop, open `RStudio` with the required packages installed (see the [data presentation intro](https://github.com/GenomicsAotearoa/metagenomics_summer_school/edit/master/materials/day4/ex16a_data_presentation_Intro.md) docs for more information).

#### 1.1 Set working directory, load *R* libraries, and import data

First, set the working directory and load the required libraries.

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/')

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

As noted during the [coverage and taxonomy](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md) exercises, it is important to remember that we currently have a table of coverage values for all *contigs* contained within each MAG. Since we're aiming to present coverage for each *MAG*, we need to reduce these contig coverages into a single mean coverage value per MAG per sample. 

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

Perform a log(2)-transformation on the coverage data (those values in the columns that `contains("sample")`). Note, here we are calculating log2(x + 1) to allow for any cases where coverage values are 0 in any of the samples (log2(0 + 1) = 0). 

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

--- 

### Part 2 - Building a heatmap of viral contigs per sample.

We can run through the same process for the viral contigs. Many of the steps are as outlined above, so we will work through these a bit quicker and with less commentary along the way. However, we will highlight a handful of differences compared to the commands for the MAGs above, for example:  steps for selecting and/or formatting the taxonomy; importing the quality output from `CheckV`; and the (optional) addition of filtering out low quality contigs.

#### 2.1 Set working directory, load *R* libraries, and import data

In this case, import the coverage, taxonomy, mapping, *and checkv* files. For taxonomy we will select the `Genome` column (converted to `Contig`), `Order_VC_predicted` (converted to `taxonomy`), and `Genus_VC_predicted` (converted to `taxonomy_genus`) (we will use the `Order_VC_predicted` to colour code viral contigs in the heat map, and `Genus_VC_predicted` to add to viral contig names in the plot). 

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/')

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

Now let's add an (optional) filtering step here to remove any contigs that `CheckV` flagged as poor quality.

First, take a look at all the categories of `checkv_quality` by using the `levels()` function.

```R
levels(collate_data_vir_log2$checkv_quality)

# 'Not-determined''Low-quality''Medium-quality''High-quality''Complete'
```

Let's filter out anything with a quality category of `Not-determined` or `Low-quality`. We will use `!=` in the `filter()` function here to return only those rows that do *not* include 'Not-determined' or 'Low-quality').

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
  RowSideColors = vir.cols[coverage.heatmap.data.vir.filt[,"taxonomy"]],
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
       col = vir.cols[1:length(levels(coverage.heatmap.data.vir.filt[,"taxonomy"]))],
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
