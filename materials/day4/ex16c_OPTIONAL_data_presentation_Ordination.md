# (*Optional*) Presentation of data: Ordinations

### Objectives

* Calculate weighted and unweighted Bray-Curtis dissimilarity metrics and nMDS analyses in `R` 
* Build nMDS plots in `R` using the `ggplot2` package.

---

### Introduction

A common method to investigate the relatedness of samples to one another is to calculate [ordinations](https://en.wikipedia.org/wiki/Ordination_(statistics)) and visualise in the form of a principal components analysis (PCA) or non-metric multidimensional scaling (nMDS) plot. In this exercise, we will calculate ordinations based on weighted and unweighted (binary) [Bray-Curtis dissimilarity](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity) and present them in nMDS plots.

The coverage tables generateed in earlier exercises have been copied to `11.data_presentation/coverage/` a for use in these exercises.

In addition to this, a simple mapping file has also been created (`11.data_presentation/coverage/mapping_file.txt`). This is a tab-delimited file listing each sample ID in the first column, and the sample "Group" in the second column (*Group_A*, *Group_B*, and *Group_C*). This grouping might represent, for example, three different sampling sites that we want to compare between. If you had other data (such as environmental measures, community diversity measures, etc.) that you wish to incorporate in other downstream analyses (such an fitting environmental variables to an ordination) you could also add new columns to this file and load them at the same time.

*NOTE: It is usually necessary to normalise coverage values based on sample depth. For example, by normalising all coverages for each sample based on either minimum or average library size (the number of sequencing reads per sample). In this case, the mock metagenome data we have been working with were already of equal depth, and so we will omit any normalisation step here. But this would not be the case with your own sequencing data sets.*

---

### 1. Import and wrangle data in *R*

*NOTE: You will recognise that the first few steps will follow the same process as the previous exercise on [generating coverage heatmaps](https://github.com/GenomicsAotearoa/metagenomics_summer_school/edit/master/materials/day4/ex16b_data_presentation_Coverage.md). In practice, these two workflows can be combined to reduce the repitition.*

#### 1.1 Set working directory, load *R* libraries, and import data

First, set the working directory and load the required libraries.

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/')

# tidyverse libraries 
#library(tidyr)
library(dplyr)
library(readr)
#library(stringr)
#library(tibble)
library(ggplot2)

# Other libraries
#library(gplots)
library(vegan)
```

Import the coverage tables and mapping file. When importing the files, we will select only the information of use here and will do a few basic first data wrangling steps. From the coverage table, we will select the `contigName` column and each of the columns of coverage values (columns `sample[1-4].bam`). 

~~For taxonomy, we will select the `user_genome` column (converted to `Bin` ID) and `classification` (coverated to `taxonomy`), and will also use `gsub` to extract just the taxonomic ranks of interest (in this case, we will extract the *class* to colour code MAGs in the heat map, and *genus* to add to MAG names in the plot), and to add `Unassigned` to any MAGs not assigned at these ranks. (*NOTE: to view a different taxonomic rank, you will need to change the two `mutate(taxonomy_class = gsub(...))` rows below accordingly*).~~

```R
# bins coverage table
cov_MAG <- read_tsv("coverage/bins_cov_table.txt") %>% 
  select(c('contigName', ends_with('.bam'))) 

# viruses coverage table
cov_vir <- read_tsv("coverage/viruses_cov_table.txt") %>% 
  mutate(Contig = contigName) %>%
  select(c('Contig', ends_with('.bam'))) 

## bins taxonomy table
#tax_MAG <- read_tsv("coverage/gtdbtk.bac120.summary.tsv") %>% 
#  mutate(Bin = gsub("(.*).filtered", "\\1", .$user_genome)) %>% 
#  mutate(taxonomy = gsub(".*;c(.*);o.*", "class\\1", classification)) %>% 
#  mutate(taxonomy = gsub("^class__$", "Unassigned", taxonomy)) %>% 
#  mutate(taxonomy_genus = gsub(".*;g__(.*);.*", "\\1", classification)) %>% 
#  mutate(taxonomy_genus = gsub("^$", "Unassigned", taxonomy_genus)) %>% 
#  select(c('Bin', 'taxonomy', 'taxonomy_genus'))

# mapping file (import both columns as factors: col_types set to factor, factor)
map.df <- read_tsv("coverage/mapping_file.txt", col_types = "ff")
```

#### 1.2 wrangle data

As noted during the [coverage and taxonomy](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md) exercises, for the bin data it is important to remember that we currently have a table of coverage values for all *contigs* contained within each bin (MAG). Since we're aiming to present coverage for each *MAG*, we need to reduce these contig coverages into a single mean coverage value per MAG per sample. 

In the following code, we first strip the `.bam` extensions off of our sample names. For the bin data, we will then leverage the fact that we added bin IDs to each of the contig headers earlier to re-extract the bin ID for each using `gsub`, use the `group_by()` function to group by `Bin`, and the `summarise()` function to return the per-sample mean coverage of each set of contigs contained within each bin (these steps are not necessary for the viral contig data).

```R
# Bins data: Extract BinID, group by Bin, calculate mean coverages for each set of contigs per Bin
cov_MAG <- cov_MAG %>%
  rename_at(
    vars(contains("sample")),
    list(~ str_replace(., ".bam", ""))
  ) %>%
  mutate(Bin = gsub("(.*)_NODE.*", "\\1", .$contigName)) %>% 
  group_by(Bin) %>% 
  summarise(across(where(is.numeric), mean))

# Viral contigs data: Strip .bam from sample names
cov_vir <- cov_vir %>%
  rename_at(
    vars(contains("sample")),
    list(~ str_replace(., ".bam", ""))
  )
```

Finally, generate a data.frame in the format required to calculate the Bray-Curtis dissimilarities via the `vegan` function `vegdist`. Here we are selecting only the coverage columns (`contains("sample")`), and then transposing so that each sample is a row and variables (bins or viral contigs) are in columns.  

```R
# Make data.frame for nMDS: MAGS
coverage.nmds.data.MAG <- t(as.data.frame(select(cov_MAG, contains("sample"))))

# Make data.frame for nMDS: viruses
coverage.nmds.data.vir <- t(as.data.frame(select(cov_vir, contains("sample"))))
```

### 2. Calculate weighted and unweighted Bray-Curtis dissimilarity and nMDS using *R*

It is often useful to examine ordinations based on both weighted and unweighted (binary) dissimilarity (or distance) metrics. Weighted metrics take into account the proportions or abundances of each variable (in our case, the coverage value of each bin or viral contig). This can be particularly useful for visualising broad shifts in overal community structure (while the membership of the community may remain relatively unchanged). Unweighted metrics are based on presence/absence alone, and can be useful for highlighting cases where the actual membership of communities differs (ignoring their relative proportions within the communities). 

Here we will use the functions `vegdist()` and `metaMDS()` from the `R` package `vegan` to generate weighted and unweighted Bray-Curtis dissimilarity matrices and nMDS solutions for the microbial bin data and viral contigs data. 

*NOTE: you may also wish to make use of the `set.seed()` function before each calculation to ensure that you obtain consistent results if the same commands are re-run at a later date.*

```bash 
# Bins: weighted Bray-Curtis
cov.bray.MAG <- vegdist(coverage.nmds.data.MAG, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.sol.MAG <- metaMDS(cov.bray.MAG, k=2, trymax=999)

# Bins: unweighted (binary) Bray-Curtis
cov.bray.binary.MAG <- vegdist(coverage.nmds.data.MAG, method="bray", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.binary.sol.MAG <- metaMDS(cov.bray.binary.MAG, k=2, trymax=999)

# Bins: weighted Bray-Curtis
cov.bray.vir <- vegdist(coverage.nmds.data.vir, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.sol.vir <- metaMDS(cov.bray.vir, k=2, trymax=999)

# Bins: unweighted (binary) Bray-Curtis
cov.bray.binary.vir <- vegdist(coverage.nmds.data.vir, method="bray", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.binary.sol.vir <- metaMDS(cov.bray.binary.vir, k=2, trymax=999)
```

---

### 3. Build nMDS plots in *R* using the *ggplot2* package.

#### 3.1 Select the data set

From here, the process is identical for each of the different analyses calculated above (based on weighted and unweighted Bray-Curtis for both microbial bin (MAG) data and viral contig data). In the first step, we are simply setting up which data set we are wishing to plot (bin or viral data, weighted or unweighted Bray-Curtis dissimilarity). The same process can then be re-run for each data set by changing *both* of the data.frames being read into these `bray.dist` and `bray.sol` variables.

```bash
bray.dist <- cov.bray.MAG
bray.sol <- cov.bray.sol.MAG
```

#### 3.2 Create the nMDS data.frame

Now, create a data.frame that includes each of the nMDS X and Y points (scores) for each sample, and merge with the mapping file to add the sample group information (i.e. whether samples are from GroupA, GroupB, and GroupC).

```bash
sol.scrs <- data.frame(
  SampleID = dimnames(bray.sol$points)[[1]],
  NMDS1 = bray.sol$point[,1],
  NMDS2 = bray.sol$point[,2]
) %>%
merge(., map, by = "SampleID", type = "full") 
```

#### 3.3 build the plot

Finally, build the nMDS using the [ggplot2](https://ggplot2.tidyverse.org/) package. Note that `ggplot`s are based on the principle of layering various aspects of a plot on top of each other in sequential calls. Getting familair with the functionality of `ggplot`s is incredibly useful for visualising many different types of data sets in a variety of different formats. 

First, set the plot theme elements (via the `theme()` function) and the colour palette for the sample grouping. 

```
theme_Ordination <- theme(
  axis.text = element_text(size=7),
  axis.title = element_text(size=7),
  axis.line = element_line(color="grey40", size=0.3, linetype=1),
  axis.ticks = element_line(color="grey40", size=0.3, linetype=1),
  legend.text = element_text(size=7),
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.position="bottom",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill="white", colour="grey40", size=0.3, linetype=1),
  plot.title = element_text(size=7, face="bold"),
  strip.background = element_rect(fill="grey90", colour="grey40", size=0.3, linetype=0),
  strip.text.x = element_text(size=7),
  panel.spacing = unit(0.2, "lines"),
  plot.margin=unit(c(5, 5, 5, 5), "points")
)

group.cols = c('#71dfdf', '#3690b8', '#00429d')
```

Now, build the plot. In this function, we:

* load the data into `ggplot()`, setting `x` and `y` axes data as the `NMDS1` and `NMDS2` columns from the `sol.scrs` data.frame
* plot the sample points using the `geom_point()` function, and also set the colour aesthetic to be based on the sample group (`Group`) from the mapping file information
* set `coord_fixed()` to ensure the the scales of the x- and y-axes are kept consistent
* use the `geom_text` function to add a text label of the [stress value](https://www.researchgate.net/post/What_is_the_importanceexplanation_of_stress_values_in_NMDS_Plots) for this ordination (extracted from `bray.sol$stress`) to the top right corner of the plot
* set the colours for the sample groups using the `scale_colour_manual()` function. Here we will pass the `group.cols` variable generated above (which is set to have a consistent colour scheme with the figures generated in the previous exercise on coverage heatmaps).
* Add x- and y-axis labels
* Modify various visual aspects of the plot using the `theme()` function (as set in the `theme_Ordination` variable above)
* Save this into the variable `NMDS.plot`
* *NOTE: here we are also surrounding the entire call in parentheses (`(...)`). This tells `R` to both save the plot in the variable `NMDS.plot` **and** also output the plot to the visualisation pane in `RStudio` or the `Jupyter Notebook`. Omitting the parentheses would result in saving to `NMDS.plot` but not also opening the plot.*

```bash
(NMDS.plot <- ggplot(sol.scrs, aes(NMDS1, NMDS2)) +
    geom_point(size = 2.5, aes(col = Group)) +
    coord_fixed() +
    geom_text(data=NULL, x=Inf, y =Inf, vjust=2.4, hjust=1.2, size=2.5, colour="black", fontface="italic", label=paste("Stress =", round(bray.sol$stress, digits=4))) +
    scale_color_manual(values = salinity.col) +
    ylab("nMDS2") + xlab("nMDS1") +
    theme_Ordination)
```

#### 3.4 Write the plot to file

*NOTE: ensure to change the file output name in the first line below for each different data set run through the above commands*

```bash
tiff("nMDS_MAGs_weighted.tiff", width=17, height=17, units="cm", res=300)
NMDS.plot
dev.off()
```

Repeat the steps in "**3. Build nMDS plots in *R* using the *ggplot2* package**" above, each time inputting a different data set (bin or viral data, weighted or unweighted Bray-Curtis dissimilarity) into the `bray.dist` and `bray.sol` variables in the section, "**3.1 Select the data set**".

*NOTE: It is often valuable to follow these visualisations up with tests for beta-dispersion (whether or not sample groups have a comparable **spread** to one another) and, provided that beta-dispersion is not significantly different between groups, PERMANOVA tests (the extent to which the variation in the data can be explained by a given variable (such as sample groups or other environmental factors). Beta-dispersion can be calculated using the  `betadisper()` function from the `vegan` package (passing the `bray.dist` data and the `map$Group` variable to group by), followed by `anova()`, `permutest()`, or `TukeyHSD()` tests of differences between the groups on the generated `betadisper` output. PERMANOVA tests can be conducted via the `adonis()` function from the `vegan` package (e.g. via: `adonis(bray.dist ~ Group, data=map, permutations=999, method="bray")`.*

---
