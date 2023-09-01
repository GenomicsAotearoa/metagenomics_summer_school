# Presentation of data: CAZy annotations heatmap

!!! info "Objectives"

    * [Build a basic heatmap from `BLAST` data using `R`](#build-a-basic-heatmap-from-annotation-data-using-r)
    * [Import and wrangle data in `R`](#import-the-data-into-an-r-dataframe)
    * [Build the plot in `R`](#build-the-plot-in-r)


---

## Build a basic heatmap from annotation data using *R*

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login) and open `RStudio`.

For this exercise, set `11.data_presentation/cazy_heatmap/` as the working directory. These annotation files we will be using have been pre-computed by annotating the `prodigal` gene predictions against the [CAZy](http://www.cazy.org/) database using the [dbCAN](http://bcb.unl.edu/dbCAN2/) resource. Briefly, each annotation was made by:

1. Annotating each `prodigal` output against the **dbCAN** database using `hmmer`
1. Converting the raw `hmmer` output into a table using the `hmmscan-parser.py` script that bundles with **dbCAN**

### Import the data into an R data.frame

The first thing we need to do is import these annotations into `R`. We will do this using the following workflow

!!! quote ""
    1. Import `tidyverse` libraries for manipulating data tables
    1. Create an empty data.frame (master data.frame)
    1. Importing each of the text file into an `R`, then 
        - Append a new column to the text file, consisting of the name of the bin
        - Select only the bin name and gene annotation from the imported table
        - Append the text table to the master data.frame

First, we import our `R` libraries with the `library()` command. For this workflow, we need three libraries from the `tidyverse` package:

!!! r-project "code"

    ```R
    setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/cazy_heatmap/')

    library(dplyr)
    library(tidyr)
    library(tibble)
    ```

We can then import our data using the `list.files()` function to loop over each text file, and the `mutate`, `select`, and pipe (`%>%`) functions from the `dplyr` library.

!!! r-project "code"

    ```R
    cazy_files <- list.files('.')

    # For each file, import it, drop unneeded columns, and add a column recording the bin name
    cazy_df <- data.frame()

    for( cazy_file in cazy_files ) {
        df <- read.table(cazy_file, sep='\t', stringsAsFactors=F, header=F) %>% 
             mutate('Bin' = cazy_file) %>%
             select(CAZy=V1, Bin)
        cazy_df <- rbind(cazy_df, df)
    }
    ```

We can inspect the final data.frame using the `head` command:

!!! r-project "code"

    ```R
    head(cazy_df)
    ```

!!! circle-check "Console output"

    ```
                         CAZy                 Bin
    1                 AA3.hmm bin_0_parsed.domtbl
    2                 AA4.hmm bin_0_parsed.domtbl
    3 GT2_Glycos_transf_2.hmm bin_0_parsed.domtbl
    4                 CE1.hmm bin_0_parsed.domtbl
    5                GH23.hmm bin_0_parsed.domtbl
    6                 GT8.hmm bin_0_parsed.domtbl
    ```

We can also confirm that we have imported all of the text files by looking at the unique entries in the `Bin` column:

!!! r-project "code"

    ```R
    unique(cazy_df$Bin)
    ```

!!! circle-check "Console output"

    ```
    [1] "bin_0_parsed.domtbl" "bin_1_parsed.domtbl" "bin_2_parsed.domtbl" "bin_3_parsed.domtbl"
    [5] "bin_4_parsed.domtbl" "bin_5_parsed.domtbl" "bin_6_parsed.domtbl" "bin_7_parsed.domtbl"
    [9] "bin_8_parsed.domtbl" "bin_9_parsed.domtbl"
    ```

We will now perform a summarising step, aggregating instances of multiple genes with the same annotation into a single count for each genome. We do this by

- For each bin in the data frame
    - For each annotation in the bin
        - Count the number of times the annotation is observed

For the majority of cases this will probably be one, but there will be a few cases where multiple annotations have been seen.

This process is done using the `group_by` and `tally` functions from `dplyr`, again using pipes to pass the data between functions.

!!! r-project "code"

    ```R
    cazy_counts <- cazy_df %>% 
      group_by(Bin, CAZy) %>% 
      tally()

    cazy_counts
    ```

!!! circle-check "Console output"

    ```
    A tibble: 402 × 3
    Groups:   Bin [10]
       Bin                 CAZy          n
       <chr>               <chr>     <int>
     1 bin_0_parsed.domtbl AA3.hmm       1
     2 bin_0_parsed.domtbl AA4.hmm       1
     3 bin_0_parsed.domtbl AA6.hmm       2
     4 bin_0_parsed.domtbl CBM12.hmm     1
     5 bin_0_parsed.domtbl CBM50.hmm     2
     6 bin_0_parsed.domtbl CBM78.hmm     2
     7 bin_0_parsed.domtbl CE1.hmm       3
     8 bin_0_parsed.domtbl CE11.hmm      1
     9 bin_0_parsed.domtbl CE3.hmm       1
    10 bin_0_parsed.domtbl CE4.hmm       1
    … with 392 more rows
    ```

We now have a data.frame-like object (a [tibble](https://tibble.tidyverse.org/)) with three columns. We can convert this into a gene matrix using the `pivot_wider` function from the `tidyr` library to create a genome x gene matrix in the following form:

|Bin|CAZy_1|CAZy_2|...|CAZy_n|
|:---:|:---:|:---:|:---:|:---:|
|bin_0|N. of genes|N. of genes|...|N. of genes|
|bin_1|N. of genes|...|...|...|
|...|...|...|...|...|
|bin_9|N. of genes|...|...|...|

!!! r-project "code"

    ```R
    cazy_matrix <- cazy_counts %>% 
      pivot_wider(id_cols=Bin, names_from=CAZy, values_from=n, values_fill=list(n = 0))
    ```

### Build the plot in *R*

Finally, we create the actual plot by passing this matrix into the `pheatmap` library. Before doing this, we need to take the text column `Bin` from the matrix and move it into the rownames, as this is how `pheatmap` infers the names of our samples. Also, if we left text data in the numeric input for a heatmap, the function would crash. We can quickly transfer the `Bin` column to the rownames using the `column_to_rownames` function from the `tibble` library.

!!! r-project "code"

    ```R
    library(pheatmap)

    colours <- colorRampPalette(c("#fff9e7","#920000"), space="Lab")(100)

    cazy_matrix %>% column_to_rownames('Bin') %>% as.matrix(.) %>% pheatmap(., col = colours)
    ```
    
<center>
![image](../figures/ex15_CAZy_heatmap.png){width="700"}
</center>

And there we go. This is a pretty basic heatmap, so there are a number of presentation issues with it. If you have time, try to do the following fixes to the heatmap by exploring the manual for `pheatmap` or other `tidyverse` and `R` functions.

!!! quote ""
    1. Replace the column-wise clustering with alphabetic arrangement of the gene columns
    1. Change the colour palette of the heatmap
    1. Reduce the font size in the column (gene) labels
    1. Remove the `.hmm` extensions from the gene names, and the `.txt` extensions from the bin names
    1. Add colouring to the row (bins), marking bins belonging to the same phyla

---
