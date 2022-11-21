# Presentation of data: KEGG pathway maps

!!! info "Objectives"

    * [Build a KEGG pathway map using `R`](#build-a-kegg-pathway-map-using-r)
    * [Import and wrangle data in `R`](#set-the-working-directory-and-load-packages-and-files-into-r)
    * [Subset the data for KO IDs of interest](#subset-the-data)
    * [Reshape the data for input into pathview](#reshaping-the-data-for-input-into-pathview)
    * [Creating pathway map of genes related to nitrogen metabolism](#creating-pathway-map-of-genes-related-to-nitrogen-metabolism)

---

### Build a KEGG pathway map using *R*

In this exercise, we will generate KEGG a pathways map from genome annotations to visualize potential pathways present in our assembled, binned genome sequences.

The key package used here is [pathview](https://academic.oup.com/bioinformatics/article/29/14/1830/232698), available from `Bioconductor` (for installation on your local version of `RStudio`, see the previous [intro to data presentation](../day4/ex16a_data_presentation_Intro.md) section). `pathview` is mainly used in transcriptomic analyses but can be adapted for any downstream analyses that utilise the KEGG database. For this exercise, we are interested in visualising the prevalence of genes that we have annotated in a pathway of interest.

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login) and open `RStudio`.

#### 1. Prepare environment

Set the working directory, load the required packages, and import data.

```R
# Set working directories ----
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/kegg_map')

# Load libraries ----
# Tidyverse packages
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

# Colour package
library(viridis)

# KEGG maps
library(pathview)
library(KEGGREST)
```

Load your files into R. Here, we are loading them into a list. Given that there are ten files, it is easier to load, clean, and analyze the data using list methods available in `tidyverse`.

```R
filenames <- list.files(pattern = ".*.aa")
bin_ids <- str_extract(filenames, "bin_\\d+")
annotations <- map(filenames, function(file) read_tsv(file)) %>% 
  setNames(., bin_ids)
```

!!! note "Note"
    R will print some messages about column types and renaming columns. However, we do not need all the columns for this analysis and it is not necessary to reformat them.

#### 2. Subset the data

For this exercise, we really only require the KEGG orthology IDs (KO numbers) for each bin. There are many columns in the table with annotations pulled from many databases. Lets begin by finding out which columns correspond to each predicted gene and relevant KEGG annotations in the unified annotations table. We will use annotations for `bin_0` as an example.

```R
# Find out what columns are in the table.
names(annotations$bin_0)

# names(annotations$bin_0)
#  [1] "Query Gene"              "GC%"                     "Contig name"             "Start position"         
#  [5] "Stop position"           "Orientation"             "Query sequence"          "Signalling"             
#  [9] "Signal confidence"       "Target gene (UniProt)"   "Identity...11"           "Coverage...12"          
# [13] "E-value...13"            "Description...14"        "Taxonomy...15"           "Target gene (UniRef100)"
# [17] "Identity...17"           "Coverage...18"           "E-value...19"            "Description...20"       
# [21] "Taxonomy...21"           "Target gene (KEGG)"      "Identity...23"           "Coverage...24"          
# [25] "E-value...25"            "Description...26"        "Taxonomy...27"           "Target gene (Pfam)"     
# [29] "E-value...29"            "Description...30"        "Target gene (TIGRfam)"   "E-value...32"           
# [33] "Description...33"
```

We can see that the headers related to KEGG annotations are in columns 22 to 27, with column headers `Target gene (KEGG)`, `Identity...23`, `Coverage...24`, `E-value...25`, `Description...26`, and `Taxonomy...27`. However, lets be surgical and only get columns that match the pattern of a KO number. All KO numbers start with K and are followed by 5 digits. We can use this pattern to find the required column.

```R
map(annotations$bin_0, function(column) {
  any(str_detect(column, "K\\d{5}"))
}) %>% 
  keep(isTRUE)

# $Description...26
# [1] TRUE
```

The code above goes through each column in the annotation table of bin 0, then detects if any string in the column contains the KO pattern. We then follow through by keeping only those columns that report `TRUE` in detecting the requisite pattern.

We can move on by only selecting for the gene ID (here, `Query Gene`) and the column where KO numbers are. The code below also creates a separate column `KO` that contains only string that matches the KO pattern.

```R
KEGG_annotations <- map(annotations, function(data) {
  data %>% 
    # Select columns for gene ID and KEGG annotations
    select(`Query Gene`, `Description...26`) %>% 
    # Create a column for KO numbers
    mutate(
      KO = str_extract_all(`Description...26`, "K\\d{5}")
    )
})
```

#### 3. Summarise KO per bin

Here, we are interested in the available KO in each bin. Thus, we can summarise the data by the bin to generate a list of KO per bin. Note that some annotations do not have KO numbers attached to them. In these cases, we will filter these data out.

!!! info "Multiple KOs per bin"

    Multiple annotations per bin is possible and not entirely rare, even if you did filter by E-value/bitscore. Some genes may just be very difficult to tell apart based on pairwise sequence alignment annotations. In this case, we are looking for overall trends. Our question here is: Does this MAG have this pathway? We can further refine annotations by comparing domains and/or gene trees to known, characterised gene sequences if gene annotations look suspicious.

```R
KO_bins <- map(annotations, function(data) {
  data %>% 
    # Selecting the relevant columns
    select(`Query Gene`, KO)
}) %>% 
  # Concatenate data frames into one big data frame
  bind_rows(.id = "bin_id") %>% 
  # Separate multiple KO per row into their own row
  unnest(KO) %>% 
  # Tally hits by KO and MAG
  group_by(KO, bin_id) %>% 
  tally(name = "hits") %>% 
  # Filter tallies that are not based on KO numbers (some KEGG annotations do
  # not have an assigned KO number)
  filter(str_detect(KO, "K\\d{5}")) %>% 
  # Arrange columns by MAG/bin ID and KO numbers (aesthetic purposes)
  arrange(bin_id, KO)
```

#### Identifying pathway maps of interest

Before moving on, we must first identify the pathway map ID of our pathway of interest. Lets say, for this exercise, we are interested in the TCA cycle. Here, we will use `KEGGREST` to access the KEGG database and query it with a search term. 

`KEGGREST` can also help you identify other information stored in the KEGG database. For more information, the `KEGGREST` vignette can be viewed using the `vignette` function in `R`: `vignette("KEGGREST-vignette")`

```R
keggFind(database = "pathway", query = "TCA cycle")
#               path:map00020 
# "Citrate cycle (TCA cycle)"

# We find the map ID is 00020 and assign it to an object.
tca_map_id <- "00020"
```

#### Reshaping the data for input into pathview

`pathview` needs the data as a numeric matrix with IDs as row names and samples/experiments/bins as column names. Here, we will reshape the data into a matrix of counts per KO number in each bin.

```R
KO_matrix <- pivot_wider(
  KO_bins,
  names_from = "bin_id", 
  values_from = "hits", 
  values_fill = NA
) %>% 
  column_to_rownames("KO") %>% 
  as.matrix()
```

If you click on `KO_matrix` in the Environment pane, you can see that it is now a matrix of counts per KO per bin. Bins that do not possess a particular KO number is given NA. Do not worry about that as `pathview` can deal with that.

#### Creating pathway map of genes related to TCA cycle

Now we can generate images of the KEGG pathway maps using the matrix we just made. For this section, we will try to find genes invovled in the TCA cycle.

```R
pv_bin_0 <- pathview(
  gene.data = KO_matrix[, "bin_0"],
  pathway.id = tca_map_id,
  species = "ko",
  out.suffix = "pv_bin_0"
)
```

There is no plot output for this command as it automatically writes the results into the current working directory. By default, it names the file as `<species><mapID>.<out.suffix>.png`. If this is the first time this is run, it will also write the pathway map's original image file `<species><mapID>.png` and the `<species><mapID>.xml` with information about how the pathway is connected.

Lets take a look at our first output.

![image](../figures/day4_keggmap.ko00020.pv_bin_0.png)

Boxes highlighted in red means that our MAG has this gene. However, the colour scale is a little strange seeing as there cannot be negative gene annotation hits (its either NA or larger than 0). Also, we know that there are definitely bins with more than 1 of some KO, but the colour highlights do not show that. Lets tweak the code further and perhaps pick better colours. For the latter, we will use the `viridis` colour package that is good for showing a gradient.

```R
# Set colours
path_colours <- viridis(n = 3, begin = 0.65, end = 1, direction = 1)

# For more information on the viridis package: 
# vignette("intro-to-viridis")

# Plot pathway
pv_bin_0.2 <- pathview(
  gene.data = KO_matrix[, "bin_0"],
  pathway.id = tca_map_id,
  species = "ko",
  # Lets make an arbitrary assumption that 5 copies is a lot
  limit = list(
    gene = c(1, 5),
    cpd = c(1, 5)
  ),
  bins = list(
    gene = 4,
    cpd = 4
  ),
  # We are plotting number of hits, so specify TRUE for this
  # If plotting, say, gene/transcript abundance, set this to FALSE
  discrete = list(
    gene = TRUE,
    cpd = TRUE
  ),
  # Tally colours
  low = path_colours[1],
  mid = path_colours[2],
  high = path_colours[3],
  out.suffix = "pv_bin_0.2"
)
```

![image](../figures/day4_keggmap.ko00020.pv_bin_0.2.png)

This plot looks much better. We can see that some genes do have more hits than others. Now, lets propagate this using `map(...)` based on our bin IDs.

```R
pv_bin_all <- map(bin_ids, function(bin) {
  # Get column with correct bin ID
  bin_data <- KO_matrix[, bin]
  # Prepare output suffix
  out.suffix = paste0("TCA.", bin)
  # Plot
  pathview(
    gene.data = bin_data,
    pathway.id = tca_map_id,
    species = "ko",
    # Lets make an arbitrary assumption that 5 copies is a lot
    limit = list(
      gene = c(1, 5),
      cpd = c(1, 5)
    ),
    bins = list(
      gene = 4,
      cpd = 4
    ),
    # We are plotting number of hits, so specify TRUE for this
    # If plotting, say, gene/transcript abundance, set this to FALSE
    discrete = list(
      gene = TRUE,
      cpd = TRUE
    ),
    # Tally colours
    low = path_colours[1],
    mid = path_colours[2],
    high = path_colours[3],
    out.suffix = out.suffix
  )
})
```

!!! success "Results"
    
    === "Bin 0"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_0.png)

    === "Bin 1"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_1.png)

    === "Bin 2"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_2.png)

    === "Bin 3"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_3.png)
        
    === "Bin 4"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_4.png)

    === "Bin 5"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_5.png)

    === "Bin 6"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_6.png)

    === "Bin 7"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_7.png)

    === "Bin 8"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_8.png)

    === "Bin 9"

        ![image](../figures/day4_keggmap.ko00020.TCA.bin_9.png)

Based on the plots, it seems that not all bins have a complete TCA cycle. Bins 4, 6, and 7 have a substantially truncated cycle.

Now that you know how to make pathway maps, try it using different pathways of interest!

---
