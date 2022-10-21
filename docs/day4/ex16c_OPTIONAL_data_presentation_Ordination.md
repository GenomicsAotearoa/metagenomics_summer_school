# (*Optional*) Presentation of data: Ordinations

!!! info "Objectives"

    * [Import and wrangle data in `R`](#1-import-and-wrangle-data-in-r)
    * [Calculate weighted and unweighted Bray-Curtis dissimilarity metrics and nMDS analyses](#2-calculate-weighted-and-unweighted-bray-curtis-dissimilarity-and-nmds-using-r)
    * [Build nMDS plots in `R` using the `ggplot2` package](#3-build-nmds-plots-in-r-using-the-ggplot2-package)
    * [Discussion: Follow-up analyses](#follow-up-analyses)

---

### Introduction

A common method to investigate the relatedness of samples to one another is to calculate [ordinations](https://en.wikipedia.org/wiki/Ordination_(statistics)) and to visualise this in the form of a principal components analysis (PCA) or non-metric multidimensional scaling (nMDS) plot. In this exercise, we will calculate ordinations based on weighted and unweighted (binary) [Bray-Curtis dissimilarity](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity) and present these in nMDS plots.

The coverage tables generated in earlier exercises have been copied to `11.data_presentation/coverage/` a for use in these exercises.

In addition to this, a simple mapping file has also been created (`11.data_presentation/coverage/mapping_file.txt`). This is a tab-delimited file listing each sample ID in the first column, and the sample "Group" in the second column (*Group_A*, *Group_B*, and *Group_C*). This grouping might represent, for example, three different sampling sites that we want to compare between. If you had other data (such as environmental measures, community diversity measures, etc.) that you wish to incorporate in other downstream analyses (such an fitting environmental variables to an ordination) you could also add new columns to this file and load them at the same time.

*NOTE: As discussed in the [coverage and taxonomy exercises](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md), it is usually necessary to normalise coverage values across samples based on equal sequencing depth. This isn't necessary with the mock metagenome data we're working with, but if you include this step in your own work you would read the **normalised** coverage tables into the steps outlined below.*

---

### 1. Import and wrangle data in *R*

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login) and open a `Notebook` running the `R 4.0.1` module as the kernel (or, outside the context of this workshop, open `RStudio` with the required packages installed (see the [data presentation intro](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day4/ex16a_data_presentation_Intro.md) docs for more information)).

*NOTE: You will recognise that the first few steps will follow the same process as the previous exercise on [generating coverage heatmaps](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day4/ex16b_data_presentation_Coverage.md). In practice, these two workflows can be combined to reduce the repititive aspects.*

#### 1.1 Set working directory, load *R* libraries, and import data

First, set the working directory and load the required libraries.

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/')

# tidyverse libraries 
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# Other libraries
library(vegan)
```

Import the coverage tables and mapping file. When importing the files, we will select only the information of use here and will do a few basic first data wrangling steps. From the coverage table, we will select the `contigName` column and each of the columns of coverage values (columns `sample[1-4].bam`). 

```R
# bins coverage table
cov_MAG <- read_tsv("coverage/bins_cov_table.txt") %>% 
  select(c('contigName', ends_with('.bam'))) 

# viruses coverage table
cov_vir <- read_tsv("coverage/viruses_cov_table.txt") %>% 
  mutate(Contig = contigName) %>%
  select(c('Contig', ends_with('.bam'))) 

# mapping file (import both columns as factors: col_types set to factor, factor)
map.df <- read_tsv("coverage/mapping_file.txt", col_types = "ff")
```

#### 1.2 wrangle data

As noted during the [coverage and taxonomy](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day3/ex11_coverage_and_taxonomy.md) exercises, for the bin data it is important to remember that we currently have a table of coverage values for all *contigs* contained within each bin (MAG). Since we're aiming to present coverage for each *MAG*, we need to reduce these contig coverages into a single mean coverage value per MAG per sample. 

In the following code, we first strip the `.bam` extensions off of our sample names. For the bin data, we will then leverage the fact that we added bin IDs to each of the contig headers earlier to re-extract the bin ID for each using `gsub`, use the `group_by()` function to group by `Bin`, and the `summarise()` function to return the per-sample mean coverage of each set of contigs contained within each bin (these grouping steps are not necessary for the viral contig data).

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

It is often useful to examine ordinations based on both weighted and unweighted (binary) dissimilarity (or distance) metrics. Weighted metrics take into account the proportions or abundances of each variable (in our case, the coverage value of each bin or viral contig). This can be particularly useful for visualising broad shifts in overall community structure (while the membership of the community may remain relatively unchanged). Unweighted metrics are based on presence/absence alone, and can be useful for highlighting cases where the actual membership of communities differs (ignoring their relative proportions within the communities). 

Here we will use the functions `vegdist()` and `metaMDS()` from the `R` package `vegan` to generate weighted and unweighted Bray-Curtis dissimilarity matrices and nMDS solutions for the microbial bin data and viral contigs data. 

*NOTE: you may also wish to make use of the `set.seed()` function before each calculation to ensure that you obtain consistent results if the same commands are re-run at a later date.*

*NOTE: in the `Jupyter` environment, these commands will create long outputs below each code block. If need be, these outputs can be cleared from the screeen via `right-click on the code block > Clear Outputs`.*

```bash 
# Bins: weighted Bray-Curtis
cov.bray.MAG <- vegdist(coverage.nmds.data.MAG, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.sol.MAG <- metaMDS(cov.bray.MAG, k=2, trymax=999)
```

```bash
# Bins: unweighted (binary) Bray-Curtis
cov.bray.binary.MAG <- vegdist(coverage.nmds.data.MAG, method="bray", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.binary.sol.MAG <- metaMDS(cov.bray.binary.MAG, k=2, trymax=999)
```

```bash
# Viruses: weighted Bray-Curtis
cov.bray.vir <- vegdist(coverage.nmds.data.vir, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.sol.vir <- metaMDS(cov.bray.vir, k=2, trymax=999)
```

```bash
# Viruses: unweighted (binary) Bray-Curtis
cov.bray.binary.vir <- vegdist(coverage.nmds.data.vir, method="bray", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE)
cov.bray.binary.sol.vir <- metaMDS(cov.bray.binary.vir, k=2, trymax=999)
```

---

### 3. Build nMDS plots in *R* using the *ggplot2* package

#### 3.1 Set the *ggplot* plot theme and sample groups colour palette

We will build the nMDS plot using the [ggplot2](https://ggplot2.tidyverse.org/) package. Note that `ggplot`s are based on the principle of layering various aspects of a plot on top of each other in sequential calls. Getting familair with the functionality of `ggplot`s is incredibly useful for visualising many different types of data sets in a variety of different formats. 

First, set the plot theme elements (via the `theme()` function), and the colour palette for the sample grouping. 

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

#### 3.2 Select the data set

From here, the process is identical for each of the different analyses calculated above (based on weighted and unweighted Bray-Curtis for both microbial bin (MAG) data and viral contig data). 

In this first step, we are simply setting up which data set we are wishing to plot (bin or viral data, weighted or unweighted Bray-Curtis dissimilarity). The same process can then be re-run for each data set by changing *both* of the data.frames being read into these `bray.dist` and `bray.sol` variables (copy the relevant data.frame names from the section above, or un-comment the relevant lines below).

```bash
bray.dist <- cov.bray.MAG
bray.sol <- cov.bray.sol.MAG

#bray.dist <- cov.bray.binary.MAG
#bray.sol <- cov.bray.binary.sol.MAG

#bray.dist <- cov.bray.vir
#bray.sol <- cov.bray.sol.vir

#bray.dist <- cov.bray.binary.vir
#bray.sol <- cov.bray.binary.sol.vir
```

#### 3.3 Create the nMDS data.frame

Now, create a data.frame that includes each of the nMDS X and Y points (NMDS1 and NMDS2 scores) for each sample, and merge with the mapping file to add the sample group information (i.e. whether samples are from *GroupA*, *GroupB*, or *GroupC*).

```bash
sol.scrs <- data.frame(
  SampleID = dimnames(bray.sol$points)[[1]],
  NMDS1 = bray.sol$point[,1],
  NMDS2 = bray.sol$point[,2]
) %>%
merge(., map.df, by = "SampleID", type = "full") 
```

#### 3.4 build the plot

Finally, build the plot using `ggplot`. In this function, we:

* load the data into `ggplot()`, setting `x` and `y` axes data as the `NMDS1` and `NMDS2` columns from the `sol.scrs` data.frame
* plot the sample points using the `geom_point()` function, and also set the colour aesthetic to be based on the sample group (`Group`) from the mapping file information
* set `coord_fixed()` to ensure the the scales of the x- and y-axes are kept consistent
* use the `geom_text` function to add a text label of the [stress value](https://www.researchgate.net/post/What_is_the_importanceexplanation_of_stress_values_in_NMDS_Plots) for this ordination (extracted from `bray.sol$stress`) to the top right corner of the plot
* set the colours for the sample groups using the `scale_colour_manual()` function. Here we will pass the `group.cols` variable generated above (which is set to have a consistent colour scheme with the figures generated in the previous exercise on coverage heatmaps).
* Add x- and y-axis labels
* Modify various visual aspects of the plot using the `theme()` function (as set in the `theme_Ordination` variable above)
* Save this into the variable `NMDS.plot`
* *NOTE: here we are also surrounding the entire call in parentheses (`(...)`). This tells `R` to both save the plot in the variable `NMDS.plot` and **also** output the plot to the visualisation pane in `RStudio` or the `Jupyter Notebook`. Omitting the parentheses would result in saving to `NMDS.plot` but not viewing the plot.*

```bash
(NMDS.plot <- ggplot(sol.scrs, aes(NMDS1, NMDS2)) +
    geom_point(size = 2.5, aes(col = Group)) +
    coord_fixed() +
    geom_text(data=NULL, x=Inf, y =Inf, vjust=2.4, hjust=1.2, size=2.5, colour="black", fontface="italic", label=paste("Stress =", round(bray.sol$stress, digits=4))) +
    scale_color_manual(values = group.cols) +
    ylab("nMDS2") + xlab("nMDS1") +
    theme_Ordination)
```

#### 3.4 Save the plot to file

*NOTE: change the file output name in the first line below for each different data set run through the above commands*

```bash
png("nMDS_MAGs_weighted.png", width=17, height=17, units="cm", res=300)
NMDS.plot
dev.off()
```

Repeat the steps above (from 3.2 onwards), each time inputting a different data set (bin or viral data, weighted or unweighted Bray-Curtis dissimilarity) into the `bray.dist` and `bray.sol` variables in the section "**3.2 Select the data set**".

---

### nMDS figures: example outputs

*NOTE: How informative these types of analyses are depends in part on the number of samples you actually have and the degree of variation between the samples. As you can see in the nMDS plots based on unweighted (binary) Bray-Curtis dissimilarities (especially for the MAGs data) there are not enough differences between any of the samples (in this case, in terms of community membership, rather than relative abundances) for this to result in a particularly meaningful or useful plot in these cases.*

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex15_fig1_nMDS_MAGs_weighted.png)

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex15_fig2_nMDS_MAGs_unweighted.png)

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex15_fig3_nMDS_Vir_weighted.png)

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex15_fig4_nMDS_Vir_unweighted.png)

---

### Follow-up analyses

It is often valuable to follow these visualisations up with tests for beta-dispersion (whether or not sample groups have a comparable *spread* to one another) and, provided that beta-dispersion is not significantly different between groups, PERMANOVA tests (the extent to which the variation in the data can be explained by a given variable (such as sample groups or other environmental factors, based on differences between the *centroids* of each group). 

Beta-dispersion can be calculated using the  `betadisper()` function from the `vegan` package (passing the `bray.dist` data and the `map$Group` variable to group by), followed by `anova()`, `permutest()`, or `TukeyHSD()` tests of differences between the groups (by inputting the generated `betadisper` output). PERMANOVA tests can be conducted via the `adonis()` function from the `vegan` package (for example, via: `adonis(bray.dist ~ Group, data=map, permutations=999, method="bray")`.

---
