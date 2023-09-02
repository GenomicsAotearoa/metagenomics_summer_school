# Introduction to data presentation

!!! info "Objectives"

    * [Overview, `RStudio`, and using `R` in the `Jupyter` environment](#overview-rrstudio-and-using-r-in-the-jupyter-environment)
    * [Data visualisation and accessibility](#data-visualisation-and-accessibility)
    * [Logging in to the NeSI `Jupyter hub`](#getting-started-logging-in-to-the-nesi-jupyter-hub)
---

## Overview, *R/RStudio*, and using *R* in the *Jupyter* environment

There are a number of powerful packages within the `R` software environment which can be used to create high quality images of genomic and metagenomic data. While each of these packages comes with its own documentation, these documents and tutorials usually assume that your data is already in some already-prepared format. Our data will almost never be in this format, though, so these exercises show a few brief examples of how we can scrape data from our existing files to create useful figures. As such, these examples are more complicated than what you would get reading the tutorials and manuals of the plotting tools, but will be transferable to your own work.

In your own work, it may be preferable to download the relevant files from NeSI (e.g. via `scp ...`, download from `Jupyter` file explorer pane) and work with them on a locally installed version of `RStudio` on your own machine. That way you have control over what packages you can install. For today, to be able to run these `R` exercises in a stable environment within the NeSI platform, we will be running an `RStudio` from within a [Jupyter](https://jupyter.org/) environment.

By now you should be very familiar with running the terminal window from within the NeSI [Jupyter hub](https://jupyter.nesi.org.nz/hub/login). In addition to the terminal, NeSI has recently integrated `RStudio` into the `Jupyter` hub. You can also opt to run R code on `Jupyter Notebooks` which also provides an interactive space that allows for mixing multiple languages within a single document, including [Markdown](https://en.wikipedia.org/wiki/Markdown), `Python`, and `R` (by default, `Markdown` and one coding language such as `R` can be used within one document, but there are add-ons available to expand this, should you wish to). `Jupyter Notebooks` can be extremely useful as a workspace that is the equivalent of an electronic "lab book".

These exercises will take place with files in the `11.data_presentation/` folder.

## Data visualisation and accessibility

In this section, we will work through a number of example exercises for visualising various aspects of metagenomic data.

As the fundamental point of data visualisation is *communication*, when building figures it is important to be mindful of aspects of your figures that might affect the accessibility of what you're trying to communicate (i.e. to maximise the number of people you will be communicating effectively with). A considerable number of your intended audience will be affected by one of the forms of colour vision deficiency (colour blindness). There are a number of useful resources available online for both selecting and testing the appropriateness of your colour selection. Some include:

* [ColorBrewer2](https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3) (select 'colorblindsafe')
* [chroma.js](https://gka.github.io/palettes/#/7|d|6e5300,7c8c00,00a63e|ffffe0,ff005e,93003a|1|1)
* Selecting and checking your colour choice using [Viz Palette](https://projects.susielu.com/viz-palette?colors=[%22#ffd700%22,%22#ffb14e%22,%22#fa8775%22,%22#ea5f94%22,%22#cd34b5%22,%22#9d02d7%22,%22#0000ff%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22achromatopsia%22)
* [An article featuring tips for visual accessibility](https://www.nature.com/articles/d41586-021-02696-z)
* Several useful colour palettes designed by Martin Krzywinski are available [here](http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container)
* [Stack Overflow community suggestions](https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible)

We have been mindful to make appropriate colour selections throughout these examples, but please do let us know if you spot any we might have overlooked.

---

## Getting started: logging in to the NeSI *Jupyter hub*

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login).

Within the `Jupyter` launcher, click on the `RStudio` button to start a session.

All of the required packages for these exercises are already installed. 

??? note "Local RStudio"
    If you are running this on your local `R` or `RStudio`, you will need to run the following:

    !!! r-project "code"
    
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

If you are new to `RStudio`, spend a few minutes familiarising yourself with the environment (also known as the "workspace"). There is a pane for the `R` console on the left that prints information on the `R` version you are using on start-up. On the top right pane, there are 2 important tabs to note: `Environment` (where you can explore your data) and `History` (all code that you have run). The bottom right pane has tabs for `Files` (where you can navigate the directory structure), `Plots` (where your plots appear), `Packages` (available packages), `Help` (all help pages and manual are shown here), and `Viewer` (some packages output interactive content and it will show up here). At the very top, there are two toolbars: the first leads to other settings and options (you can explore this on your own time), the second one has icons for (starting from the far left and hover mouse pointer to see description):

* Open new file
* Create an `R` project
* Open existing file
* Save current document
* Save all open documents
* Print current file
* Search bar to navigate open documents
* Manage workspace panes
* Options for additional plugins

To start, open a new file to start writing code in.

---
