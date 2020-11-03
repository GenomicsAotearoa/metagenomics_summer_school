# Presentation of data

### Objectives

* Overview, `RStudio`, and using `R` in the `Jupyter` environment
* Data visualisation and accessibility

---

### Overview, *R/RStudio*, and using *R* in the *Jupyter* environment

There are a number of powerful packages within the `R` software environment which can be used to create high quality images of genomic and metagenomic data. While each of these packages comes with its own documentation, these documents and tutorials usually assume that your data is already in some already-prepared format. Our data will almost never be in this format, though, so these exercises show two brief examples of how we can scrape data from our existing files to create useful figures. As such, these examples are more complicated than what you would get reading the tutorials and manuals of the plotting tools, but will be transferable to your own work.

In your own work, it may be preferable to download the relevant files from NeSI (e.g. via `scp ...`) and work with them on a locally installed version of `RStudio` on your own machine. For today, to be able to run these `R` exercises in a stable environment within the NeSI platform, we will be running an `R` [kernel](https://en.wikipedia.org/wiki/Kernel_(operating_system)) from within a [Jupyter](https://jupyter.org/) environment. 

By now you should be very familiar with running the terminal window from within the NeSI [Jupyter hub](jupyter.nesi.org.nz/hub/login). In addition to the terminal, `Jupyter Notebooks` more generally also provide an interactive space that allows for mixing multiple languages within a single document, including [Markdown](https://en.wikipedia.org/wiki/Markdown), `Python`, and `R` (by default, `Markdown` and one coding language such as `R` can be used within one document, but there are add-ons available to expand this, should you wish to). `Jupyter Notebooks` can be extremely useful as a workspace that is the equivalent of an electronic "lab book". Today, we will be using it as an interactive space to run `R`. Note that, while the layout will be slightly different to `RStudio`, the commands we will be working through will work the same in both environments.

These exercises will take place with files in the `11.data_presentation/` folder.

---

### Data visualisation and accessibility

In this section, we will work through a number of example exercises for visualising various aspects of metagenomic data. 

As the fundamental point of data visualisation is *communication*, when building figures it is important to be mindful of aspects of your figures that might affect the accessibility of what you're trying to communicate (i.e. to maximise the number of people you will be communicating effectively with). A considerable number of your intended audience will be affected by one of the forms of colour vision deficiency (colour blindness). There are a number of useful resources available online for both selecting and testing the appropriateness of your colour selection. Some include:

* [ColorBrewer2](https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3) (select 'colorblindsafe')
* [chroma.js](https://gka.github.io/palettes/#/7|d|6e5300,7c8c00,00a63e|ffffe0,ff005e,93003a|1|1)
* Selecting and checking your colour choice using [Viz Palette](https://projects.susielu.com/viz-palette?colors=[%22#ffd700%22,%22#ffb14e%22,%22#fa8775%22,%22#ea5f94%22,%22#cd34b5%22,%22#9d02d7%22,%22#0000ff%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22achromatopsia%22)
* [Blog post](https://bconnelly.net/posts/creating_colorblind-friendly_figures/) by Brian Connelly.
* Several useful colour palettes designed by Martin Krzywinski are available [here](http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container)

We have been mindful to make appropriate colour selections throughout these examples, but please do let us know if you spot any we might have overlooked.

---

### Getting started: logging in to the NeSI *Jupyter hub*

To get started, if you're not already, log back in to NeSI's [Jupyter hub](https://jupyter.nesi.org.nz/hub/login). 

Within the `Jupyter` launcher, this time open a `Notebook` running the `R 4.0.1` module as the kernel. 

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

Spend a few minutes familiarising yourself with the `Jupyter Notebook` workspace, and how it differs to the standard terminal we've been working in. You'll see the coding language kernel running in the background of this `Notebook` in the top right of the pane. The main document works in *blocks*; click the `+` button to add additional blocks. The *Code* drop-down menu allows you to select between whether the current block is a `Markdown` or `code` (in this case, `R`) block. 

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
