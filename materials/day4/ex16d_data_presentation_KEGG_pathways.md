# Presentation of data: KEGG pathway maps

### Objectives

* Build a KEGG pathway map using `R`

---

### Build a KEGG pathway map using *R*

In this exercise, we will generate KEGG a pathways map from genome annotations to visualize potential pathways present in our assembled, binned genome sequences. 

The key package used here is [pathview](https://academic.oup.com/bioinformatics/article/29/14/1830/232698), available from `Bioconductor` (for installation on your local version of `RStudio`, see the previous [ex16a_data_presentation_Intro.md](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day4/ex16a_data_presentation_Intro.md) section). `pathview` is mainly used in transcriptomic analyses but can be adapted for any downstream analyses that utilise the KEGG database. For this exercise, we are interested in visualising the prevalence of genes that we have annotated in a pathway of interest.

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login) and open a `Notebook` running the `R 4.0.1` module as the kernel (or, outside the context of this workshop, open `RStudio` with the required packages installed (see the [data presentation intro](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day4/ex16a_data_presentation_Intro.md) docs for more information)).

#### Set the working directory and load packages and files into *R*

Set the working directory and load the required packages

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/kegg_map')

library(pathview)
library(KEGGREST)

# tidyverse libraries
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
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

We can see that there are a total of 66 genes with multiple annotations. We will split them into their own rows. 

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

We need to first identify the KEGG pathway ID. This is where `KEGGREST` comes in handy. `KEGGREST` can also help you identify other information stored in the KEGG database. For more information, the `KEGGREST` vignette can be viewed using the `vignette` function in `R`: `vignette("KEGGREST-vignette")`

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
  low ='#ffdfdc',
  mid = '#d57c69',
  high = '#980000'
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
    low ='#ffdfdc',
    mid = '#d57c69',
    high = '#980000'
  )
  rm(i, out.suffix)
}
```

The above code will generate one png image file per bin, and will save it to the current working directory. The list `pv_all_bins` in the first line contains some plotting information. Importantly, it has the KO numbers it used to match the data provided to the pathway. If you wish to, you can subset your original data to only obtain those KO numbers related to the pathway that you have plotted for further analyses with other annotations.

A few example outputs for bins 3 and 5:

#### Bin 3 nitrogen metabolism KEGG map

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex15_KEGG_maps_ko00910_bin_3.png)

#### Bin 5 nitrogen metabolism KEGG map

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex15_KEGG_maps_ko00910_bin_5.png)

---
