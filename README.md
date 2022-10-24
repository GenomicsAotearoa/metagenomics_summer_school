# Metagenomics Summer School


## <img src="./docs/theme_images/link2material.png" alt="drawing" width="20"> **Material** https://genomicsaotearoa.github.io/metagenomics_summer_school/

---

## Useful locations and links

#### Working directories

For all exercises after the `bash` introduction, you will be working from the file path

```bash
/nesi/nobackup/nesi02659/MGSS_U/YOUR_USERNAME/
```

#### bash/slurm cheatsheet

A few helpful commands and shortcuts for working in `bash` or with `slurm` can be found [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/command_line_shortcuts.md).

#### Snapshots of results to download

If you are having trouble downloading files using `scp`, we are providing exemplar output files, which you can download through your browser, [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/materials/resources), or via the following links: 

* [FastQC results](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/fastqc_results.zip)
* [Quast results](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/quast_results.zip)
* [Input files for VizBin](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/vizbin_files.zip)
* [Tables of gene annotations](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/example_annotation_tables.zip)
* [DRAM output files](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/DRAM_results.zip)

#### Slides for workshop

You can find a copy of the slides presented during the workshop, with published figures removed, in the [slides/](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/slides) folder.

#### Etherpad

We will be using a collaborative document to share long code, results, or even code errors. This document can be found [here](https://etherpad.wikimedia.org/p/MGSS_17_11_20) or type https://etherpad.wikimedia.org/p/MGSS_17_11_20 in your browser.



---


## Timetable

### Day 1 - 17<sup>th</sup> November 2020

|Time|Event|Session leader|
|:---|:---|:---|
|9:00 am – 9:30 am|**Introduction**<br>- Welcome<br>- Logging into NeSI|Michael Hoggard|
|9:30 am – 10:30 am|**TASK:** Bash scripting|Ngoni Faya|
|10:30 am – 10:50 am|**Morning tea break**||
|10:50 am – 11:30 am|**TASK:** Bash scripting (continued)|Ngoni Faya|
|11:30 am – 12:00 pm|**TALK:** The metagenomics decision tree|Kim Handley|
|12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
|12:45 pm – 1:45 pm|**TALK:** Quality filtering raw reads<br>**TASK:** Visualisation with *FastQC*<br>**TASK:** Read trimming and adapter removal<br>**TASK:** Diagnosing poor libraries<br>**TASK:** Common issues and best practice|Carmen Astudillo-Garcia|
|1:45 pm – 3:00 pm|**TASK:** Run IDBA_UD assembly<br>**TALK:** Assembly<br>- Choice of assemblers<br>- Considerations for parameters, and when to stop!<br>**TASK:** Exploring assembler options<br>**TASK:** Submitting jobs to NeSI via slurm<br>**TASK:** Run SPAdes assembly<br>**TASK (*Optional*):** Submitting variant assemblies to NeSI|Kim Handley|
|3:00 pm – 3:20 pm|**Afternoon tea break** (*Tea, coffee, and snacks provided*)||
|3:20 pm – 5:00 pm|**TALK:** Future considerations - co-assembly vs. single assemblies<br>**TASK:** Assembly evaluation<br>**TASK:** Short contig removal|Michael Hoggard|


### Day 2 - 18<sup>th</sup> November 2020

|Time|Event|Session leader|
|:---|:---|:---|
|9:00 am – 9:15 am|**Introduction**<br>- Overview of yesterday, questions<br>- Overview of today|Carmen Astudillo-Garcia|
|9:15 am – 9:30 am|**TALK:** Automation, reproducibility, and FAIR principles|Dan Jones|
|9:30 am – 10:30 am|**Binning (part 1)**<br>**TALK:** Overview of binning history<br>- Key parameters and strategies for binning<br>**TASK:** Read mapping|Kim Handley|
|10:30 am – 10:50 am|**Morning tea break** (*Tea, coffee, and snacks provided*)||
|10:50 am – 11:20 am|**TALK:** Overview of binning history (*continued*)<br>- Key parameters and strategies for binning|Kim Handley|
|11:20 am – 12:00 pm|**Binning (part 2)**<br>**TASK:** Multi-binning strategy (*Metabat* and *Maxbin*)|Kim Handley|
|12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
|12:45 pm – 2:00 pm|**Binning (part 3)**<br>**TASK:** Bin dereplication via *DAS_Tool*<br>**TASK:** Evaluating bins using *CheckM*|Michael Hoggard|
|2:00 pm - 3:00 pm|**Binning (part 4)**<br>- Discuss additional dereplication strategies, such as *dRep*<br>- How to work with viral and eukaryotic bins<br>- Dealing with organisms which possess minimal genomes<br>**TALK:** Bin refinement<br>- Refinement strategies - |Carmen Astudillo-Garcia<br>Michael Hoggard|
|3:00 pm – 3:20 pm|**Afternoon tea break** (*Tea, coffee, and snacks provided*)||
|3:20 pm – 5:00 pm|**TALK:** Bin refinement<br>- Refinement strategies (cont) - *VizBin* and *ESOMana*<br>**TASK:** Working with *VizBin*<br>**TASK:** Submit *VIBRANT* job|Michael Hoggard|

### Day 3 - 19<sup>th</sup> November 2020

|Time|Event|Session leader|
|:---|:---|:---|
|9:00 am – 9:30 am|**Introduction**<br>- Overview of yesterday, questions<br>- Overview of today|Michael Hoggard|
|9:30 am – 10:30 am|**TASK:** Identifying viral contigs (*VIBRANT*)<br>**TALK:** Identifying viruses from metagenomic data<br>**TASK:** QC of viral contigs (*CheckV*)<br>**TASK:** Coverage calculation (*bowtie*)<br>**TASK:** Taxonomic classification (Bin taxonomy with *GTDB-TK*; viral taxonomy predictions with *vConTACT2*)|Michael Hoggard<br><br><br>David Waite|
|10:30 am – 10:50 am|**Morning tea break** (*Tea, coffee, and snacks provided*)||
|10:50 am – 11:30 am|**TALK:** Gene prediction, using *prodigal*, and other tools (*RNAmer*, *Aragorn*, etc)<br>**TASK:** Predict open reading frames and protein sequences|David Waite|
|11:30 am – 12:00 pm|**TALK:** Gene annotation (part 1)<br>**TASK:** Gene annotation using *diamond* and *hmmer*<br>**Discussion:** Evaluating the quality of gene assignment<br>**Discussion:** Differences in taxonomies (*GTDB*, *NCBI* etc)|Carmen Astudillo-Garcia|
|12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
|12:45 pm – 2:00 pm|**TALK:** Gene annotation (part 2)<br>- Using online resources (e.g. *KEGG, BioCyc, MetaCyc, HydDB, PSORT*)<br>**TASK:** View KEGG annotation in KEGG website|Christia Straub<br>Florian Pichlmueller|
|2:00 pm – 3:00 pm|**TALK:** Bin taxonomic classification<br>- Bin and species determination<br>**TASK:** View phylogenetic trait distribution (*ANNOTREE*)|David Waite|
|3:00 pm – 3:20 pm|**Afternoon tea break**||
|3:20 pm – 4:30 pm|**TASK:** MAG annotation with *DRAM*<br>**TASK:** Introduce group project goals<br>**TASK:** Dividing into working groups / get a group name<br>**TASK:** Select a goal from your project|Carmen Astudillo-Garcia|
|4:30 pm – 5:00 pm|**End of day wrap up**|Kim Handley|

### Day 4 - 20<sup>th</sup> November 2020

|Time|Event|Session leader|
|:---|:---|:---|
|9:00 am – 9:15 am|**Introduction**<br>- Overview of yesterday, questions<br>- Overview of today|Michael Hoggard|
|9:15 am – 10:00 am|**TALK:** *DRAM* results overview<br>**TASK:** Explore *DRAM* results|Carmen Astudillo-Garcia|
|10:00 am – 10:30 am|**Presentation of data**<br>**TALK:** Visualising findings (metabolism maps, heatmaps, cell schematics, gene trees, gene maps)<br>**TASK:** Coverage heatmap / Ordination (*Optional*)|Michael Hoggard|
|10:30 am – 10:50 am|**Morning tea break** (*Tea, coffee, and snacks provided*)<br>**TASK:** Workshop survey||
|10:50 am – 12:00 pm|**Presentation of data (continued)**<br>**TALK:** Visualising findings (metabolism maps, heatmaps, cell schematics, gene trees, gene maps)<br>**TASK:** KEGG metabolic pathways<br>**TASK:** Gene synteny<br>**TASK:** CAZy heatmaps (*Optional*)|Carmen Astudillo-Garcia<br>Boey Jian Sheng<br>Hwee Sze Tee|
|12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
|12:45 pm – 2:30 pm|**TASK:** Analyse data for group work<br>**TASK:** Prepare group presentation|Kim Handley|
|2:30 pm – 3:00 pm|**Present and discuss findings**<br>**TASK:** Each group to give an informal presentation of their data|Kim Handley|
|3:00 pm – 3:20 pm|**Afternoon tea break** (*Tea, coffee, and snacks provided*)||
|3:20 pm – 3:40 pm|**Present and discuss findings (continued)**<br>**TASK:** Each group to give an informal presentation of their data|Carmen Astudillo-Garcia|
|3:40 pm – 4:00 pm|**End of day wrap up**<br>- Final discussion|Kim Handley<br>Michael Hoggard|

----
