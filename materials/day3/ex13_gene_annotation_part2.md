# Gene annotation

### Objectives

* Overview of the *BLAST XML* output
* Looking at gene networks in **MEGAN**
* Tie findings to your initial hypothesis

---

###  Overview of the *BLAST XML* output

Your job from the previous session should have now completed, so we will start to examine the results. Unlike the first job, we have specified this output to be writte in the *xml* format. What this means it not particualrly important, only than to say that *xml* is a language used for writing information with a hierarchical/nested structure. For example, a piece of *xml* that describes some of the course demonstrators could be written as the follow:

```xml
<tutor-list>
  <tutor>
    <name>Kim Handley</name>
    <position>Senior Lecturer</position>
    <affiliation>University of Auckland</affiliation>
    <previous-affiliation-list>
      <previous-affiliation>University of Chicago</previous-affiliation>
      <previous-affiliation>University of California, Berkley</previous-affiliation>
    </previous-affiliations>
  </tutor>
  <tutor>
    <name>David Waite</name>
    <position>Postdoctoral Research Fellow</position>
    <affiliation>University of Auckland</affiliation>
    <previous-affiliation-list>
      <previous-affiliation>University of Queensland</previous-affiliation>
      <previous-affiliation>Ministry for Primary Industries</previous-affiliation>
    </previous-affiliation-list>
  </tutor>
  <tutor>
    <name>Christina Straub</name>
    <position>Postdoctoral Research Fellow</position>
    <affiliation>Institute of Environmental Science and Research</affiliation>
    <previous-affiliation-list>
      <previous-affiliation>Massey University</previous-affiliation>
    </previous-affiliation-list>
  </tutor>
  <tutor>
    <name>Hwee Sze Tee</name>
    <position>Doctoral Candidate</position>
    <affiliation>University of Auckland</affiliation>
  </tutor>
</tutor-list>
```

This is not necessarily easy to read by eye, especially as the number of entries and categories grow. However, consider the case for the 'previous affiliation' entry in the above. Since tutors in the list can have 0, 1, or multiple previous affiliations this data would be difficult to capture in a standardised table format.

In the case of our annotation data, NCBI's *xml* schema is a good way to store information surrounding the query/target alignment itself, and also additional metadata for the target sequence. This information is used by **MEGAN** in the next exercise.

---

### Looking at gene networks in *MEGAN*

**MEGAN** is a really powerful GUI-based tool for exploring and contextualising gene annotations. The current version is **MEGAN6** but we are going to work with **MEGAN5** in this exercise. This is because **MEGAN** provides a free and licensed version, and in **MEGAN6** one of the key features we will use today was migrated into the licensed version. If you are considering using **MEGAN** for your own analysis, we do recommend obtaining the latest version.

##### Executing *MEGAN5* from the command line

<TO DO>

##### Loading your gene and annotation files

Once **MEGAN** has loaded, we are going to use the 'File' -> 'Import From BLAST...' option.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_import_menu.PNG)

Use the file browser icon to select your *xml* file for the first entry, and *fastA* protein sequences for the second. If you have a consistent naming scheme, **MEGAN** might be able to auto-detect the correct *fastA* file, but make sure you check that it is correct.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_browser.PNG)

##### Parsing hits to KEGG ontology

Navigate to the 'KEGG' tab and turn on the 'Analyse KEGG content' feature. Make sure you also enable 'Use Built-In RefSeq Map'. Also click the 'Load GI Mapping File', and select the file *gi2kegg-Feb2015X.bin* located in the same location as the **MEGAN** executable.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_parse_kegg.PNG)

We will not toggle any extra features for now. Click the 'Apply' button to start the import process.

##### Interpreting the taxonomy view

**MEGAN** will extract taxonomy information from the annotations to build up a taxonomic overview of our MAG. For each protein, taxonomy is assigned using the Lowest Common Ancestor (LCA) method, for which each protein is assigned the consensus taxonomy of its annotations.

You can view how many genes are assigned at least node of the taxonomy tree by clicking on them.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_node_assignment.PNG)

Depending on the node you click, you will see either one or two numbers reported. The 'Ass' number is the number of genes assigned at this rank specifically. In the screenshot above, this means that there are 64 genes in the MAG that cannot be reliably classified deeper than to the phylum Cyanobacteria. The 'Sum' number refers to all genes assigned to this node, including all subordinate nodes. We can expand this node by right-clicking and selecting the 'Uncollapse' option.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_expand.PNG)

Carry this on as far as you like. The deeper into the tree you go the more branches will appear as it is harder to reliably assign the deeper taxonomic ranks across a diverse set of genes.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_expanded.PNG)

If you are interested in seeing which genes were assigned to a particular node, you can export them using the 'Extract Reads...' option.

##### Exploring the KEGG view

Along the top of the **MEGAN** browser is an icon to open the KEGG view of your annotation.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_open_kegg.PNG)

Clicking on this view will open your annotation using the KEGG mappings for your annotation. This is a powerful way to explore annotations, as it provides each gene in a metabolic map of the processes they are involved in.

Similar to in the taxonomy view, you can select individual nodes and export the sequences assigned to each function. This is less useful at the high levels, but as we drill into functionality it is a nice way to identify the genes of key metabolic interest.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_kegg_view.PNG)

##### Drilling into the KEGG annotation

Open up the 'Energy Metabolism' category in the left-hand menu, then right click on the 'Carbonfixation in photosynthetic organisms' and select the 'Show Pathway' option.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_carbon_fixing_1.PNG)

This provides a fine-detailed map of the genes in our MAG, and which functionality it can perform. Hovering the mouse over a gene block will provide contextual information, which can be explored in greater detail by right clicking the node.

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex13_carbon_fixing_2.PNG)

---

### Tie findings to your initial hypothesis

It is now time to explore the genomes and try to address your original hypothesis!

Depending on what you are looking for, you will either be trying to find gene(s) of relevance to a particular functional pathway, or the omission of genes that might be critical in function. In either case, try to use the taxonomy of each MAG to determine whether it is likely to be a worthwhile candidate for exploration, or whether it is unlikely to be involved in the function you are interested in.

---
