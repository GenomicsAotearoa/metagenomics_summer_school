# Presentation of data: Gene synteny

### Objectives

* Build a sulfur assimilation gene alignment figure to investigate gene synteny using `R`

---

### Build a sulfur assimilation gene alignment figure to investigate gene synteny using `R`

When investigating the evolution of genomes, we sometimes want to consider not only the presence/absence of genes in a genome, but also how they are arranged in an operon. For this exercise, we are going to visualise several sulfur assimilation genes from bin_4, bin_5 and bin_7, comparing their arrangement among the organisms.

For this exercise, navigate to the folder `11.data_presentation/gene_synteny/`. You have been provided with a copy of the `prodigal` gene predictions for each of the bins (`.faa` files), an annotation output table using multiple databases (`.aa` files), a small table of the annotation of some key genes of interest (`cys.txt` files), and blastn output (`blast*.txt`) comparing the genes of interest from these organisms. The annotation files were created by manually searching the annotations obtained in the previous exercises.

*NOTE: Refer to [gene_synteny_grab_GOI.md](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/resources/gene_synteny_grab_GOI.md) for more information on how the `cys.txt` files were produced.*

#### Part 1 - Parsing files in bash

We will be performing this exercise in two stages. Firstly, in `bash`, we will use `cut` and `tail` to pull out the genes of interest listed in the `*cys.txt` files from the `prodigal` files. The gene names will then be used to create a table of gene coordinates from the `prodigal` output using `grep`, `cut`, and `sed`.

For these `bash` steps, we will need to return to our logged in NeSI terminal. Switch over to a NeSI `Jupyter hub` terminal or log in to a fresh session in a new terminal. 

Navigate to the `11.data_presentation/gene_synteny/` folder, and then perform the `cut` and `tail` steps outlined above.

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/gene_synteny/

cut -f3 bin_4_cys.txt | tail -n+2 > bin_4_cys.genes
cut -f3 bin_5_cys.txt | tail -n+2 > bin_5_cys.genes
cut -f3 bin_7_cys.txt | tail -n+2 > bin_7_cys.genes
```

We now have three new files, ending with the `.genes` suffix which are simply a list of the genes that we wish to extract from the `prodigal` files. We can then use each of these files as input for `grep` to pull out the *fastA* entries that correspond to these genes.

```bash
grep -f bin_4_cys.genes bin_4.filtered.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_4_cys.coords
grep -f bin_5_cys.genes bin_5.filtered.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_5_cys.coords
grep -f bin_7_cys.genes bin_7.filtered.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_7_cys.coords
```

As previously, we will quickly go through the steps of this command:

```bash
grep -f bin_4_cys.genes bin_4.filtered.genes.faa | cut -f1,2,3,4 -d "#" | sed 's/>//g' | sed 's/#/\t/g' > bin_4_cys.coords
|                                     |                      |              |                |
|                                     |                      |              |                Write the output
|                                     |                      |              |
|                                     |                      |              Replace each '#' character with tab
|                                     |                      |
|                                     |                      For each line, replace the '>' with empty text
|                                     |
|                                     Split the results into columns delimited by the # character,
|                                     then take columns 1 - 4.
| Select the lines of the bin_4_cys.faa file that contain entries found in the bin_4_cys.genes file.
```

Check the content of the `.coords` files now. You should see something like the following:

```bash
cat bin_4_cys.coords
# bin_4_NODE_55_length_158395_cov_1.135272_128     135928          136935          1
# bin_4_NODE_55_length_158395_cov_1.135272_129     136994          137299          1
# bin_4_NODE_55_length_158395_cov_1.135272_130     137411          138322          1
# bin_4_NODE_55_length_158395_cov_1.135272_131     138413          139201          1
# bin_4_NODE_55_length_158395_cov_1.135272_132     139267          140100          1
# bin_4_NODE_55_length_158395_cov_1.135272_133     140110          140988          1
# bin_4_NODE_55_length_158395_cov_1.135272_134     140985          142073          1
```

If you recall from the previous exercise on gene prediction, we have taken the first four entries from each line of the `prodigal` output, which consists of:

1. The gene name, written as [CONTIG]\_[GENE]
1. The start position of the gene
1. The stop position of the gene
1. The orienation of the gene

We will now use these tables, together with the annotation tables to create the gene synteny view in `R`.

#### Part 2 - Producing the figure in *R*

First, move back to the [Jupyter hub Notebook](https://jupyter.nesi.org.nz/hub/login) pane where we have `R` running as the kernel (or relaunch a new `R 4.0.1` Notebook). (Outside the context of this workshop, open `RStudio` with the required packages installed (see the [data presentation intro](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day4/ex16a_data_presentation_Intro.md) docs for more information)).

There are two `R` libaries we need to load for this exercise.

##### Set working directory and load *R* libraries

```R
setwd('/nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/11.data_presentation/gene_synteny/')

library(dplyr)
library(genoPlotR)
```
    
    Attaching package: ‘dplyr’
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    Loading required package: ade4
    Loading required package: grid


##### Part 2.1 - Load coordinate files

We can now begin importing the data. First, we will import the `.coords` files, and set column names to the files.


```R
bin_4_coords = read.table('bin_4_cys.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_4_coords) = c('name', 'start', 'end', 'strand')
bin_5_coords = read.table('bin_5_cys.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_5_coords) = c('name', 'start', 'end', 'strand')
bin_7_coords = read.table('bin_7_cys.coords', header=F, sep='\t', stringsAsFactors=F)
colnames(bin_7_coords) = c('name', 'start', 'end', 'strand')
```

Take a look at the content of each of these data.frames, by entering their names into the terminal. You should notice that the coordinates occur at quite different positions between the genomes. If we were looking at complete genomes, this would indicate their position relative to the *origin of replication*, but as these are unclosed genomes obtained from MAGs, they reflect the coordinates upon their particular *contig*.

We now parse these data.frames into the *dna_seg* data class, which is defined by the `genoPlotR` library.


```R
bin_4_ds = dna_seg(bin_4_coords)
bin_5_ds = dna_seg(bin_5_coords)
bin_7_ds = dna_seg(bin_7_coords)
bin_5_ds
```


<table>
<caption>A dna_seg: 4 × 11</caption>
<thead>
	<tr><th scope=col>name</th><th scope=col>start</th><th scope=col>end</th><th scope=col>strand</th><th scope=col>col</th><th scope=col>fill</th><th scope=col>lty</th><th scope=col>lwd</th><th scope=col>pch</th><th scope=col>cex</th><th scope=col>gene_type</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>bin_5_NODE_95_length_91726_cov_0.379357_18 </td><td>16268</td><td>17257</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>bin_5_NODE_95_length_91726_cov_0.379357_19 </td><td>17261</td><td>18130</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>bin_5_NODE_95_length_91726_cov_0.379357_20 </td><td>18141</td><td>18959</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
	<tr><td>bin_5_NODE_95_length_91726_cov_0.379357_21 </td><td>19121</td><td>20119</td><td>-1</td><td>blue</td><td>blue</td><td>1</td><td>1</td><td>8</td><td>1</td><td>arrows</td></tr>
</tbody>
</table>


By looking at the table, we can see that the genes in bin_5 are in reversed order (strand -1), so we might want to change the color of the gene to make sure they are different than bin_4 and bin_7.


```R
bin_5_ds$col = "#1a535c"
bin_5_ds$fill = "#1a535c"
```

##### Part 2.3 - Load annotation tables

Then, we can load the annotation tables we have into `R` and take a look at them.


```R
bin_4_ann = read.table('bin_4_cys.txt', header=T, sep='\t', stringsAsFactors=F) %>%
cbind(., bin_4_coords) %>%
select(Annotation, start1=start, end1=end, strand1=strand)
bin_5_ann = read.table('bin_5_cys.txt', header=T, sep='\t', stringsAsFactors=F) %>%
cbind(., bin_5_coords) %>%
select(Annotation, start1=start, end1=end, strand1=strand)
bin_7_ann = read.table('bin_7_cys.txt', header=T, sep='\t', stringsAsFactors=F) %>%
cbind(., bin_7_coords) %>%
select(Annotation, start1=start, end1=end, strand1=strand)
bin_5_ann
```

<table>
<caption>A data.frame: 4 × 4</caption>
<thead>
	<tr><th scope=col>Annotation</th><th scope=col>start1</th><th scope=col>end1</th><th scope=col>strand1</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> sbp; sulfate-binding protein                                      </td><td>16268</td><td>17257</td><td>-1</td></tr>
	<tr><td> cysT; sulfate transporter CysT                                    </td><td>17261</td><td>18130</td><td>-1</td></tr>
	<tr><td> cysW; sulfate transporter CysW                                    </td><td>18141</td><td>18959</td><td>-1</td></tr>
	<tr><td> cysA; sulfate.thiosulfate ABC transporter ATP-binding protein CysA</td><td>19121</td><td>20119</td><td>-1</td></tr>
</tbody>
</table>


We need to create one more table with descriptive information. This is an *annotation* object, which contains the name of each gene in the figure along with the coordinates to write the name. x1 = starts of the gene, x2 = end of the gene


```R
annot1 <- annotation(x1 = c(bin_4_ds$start[1], bin_4_ds$start[2],bin_4_ds$start[5],bin_4_ds$start[6],bin_4_ds$start[7]),
                    x2 = c(bin_4_ds$end[1], bin_4_ds$end[4], bin_4_ds$end[5], bin_4_ds$end[6], bin_4_ds$end[7]),
                    text =  c("sbp", "unknown domain", "cysU", "cysW", "cysA"),
                    rot = 0, col = "black")
annot2 <- annotation(x1 = c(bin_7_ds$start[1], bin_7_ds$start[2],bin_7_ds$start[3],bin_7_ds$start[4]),
                    x2 = c(bin_7_ds$end[1], bin_7_ds$end[2], bin_7_ds$end[3], bin_7_ds$end[4]),
                    text = c("sbp", "cysT","cysW","cysA"),
                    rot = 0, col = "black")
```

##### Part 2.3 - Creating a comparison table

Then, we can parse the blast output as comparison file among the genomes. genoPlotR can read tabular files, either user-generated tab files (read_comparison_from_tab), or from BLAST output (read_comparison_from_blast). To produce files that are readable by genoPlotR, the -m 8 or 9 option should be used in blastall, or -outfmt 6 or 7 with the BLAST+ suite.


```R
blast1 = read_comparison_from_blast("blast_bin4_bin5.txt")
blast2 = read_comparison_from_blast("blast_bin5_bin7.txt")
```

What it does here is to set the line color according to the direction of the gene match and the color gradient is based on the percent identity of the matches. Lighter color indicates weaker match and darker color indicates strong match.

Now we can generate the plot. Running this command in `RStudio` or our `Jupyter Notebook` running the `R` kernel interactively loads the figure. 

*NOTE: the commented-out lines below (the two lines starting with `#`) will not run. Un-comment these if you wish to save the figure to file rather than opening it in the `Jupyter` or `RStudio` viewer. (N.b. alternatives using a similar command are also available other formats, including `tiff` and `png`).*

```R
#pdf("genoplot.pdf",colormodel = "cmyk",width = 8,height = 4,paper = 'special')
plot_gene_map(dna_segs = list(bin_4_ds,bin_5_ds,bin_7_ds), 
              gene_type = "arrows", dna_seg_labels = c("bin_4", "bin_5","bin_7"), 
              comparisons = list(blast1,blast2), dna_seg_scale = TRUE, 
              override_color_schemes=FALSE,annotations=list(annot1,NULL,annot2),
              annotation_height=1.7, dna_seg_label_cex = 1,main = "Sulfur assimilation")
#dev.off()
```

![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_gene_synteny_fig1.png)

*NOTE: While we do have control over the colours of the arrows via setting the `bin_n_ds$col` and `bin_n_ds$fill` parameters for each contig (as above), unfortunately there appears to be little flexibility within the `plot_gene_map()` function for setting the colours of the segments joining the arrows (the current options are limited to 'red_blue', 'blue_red', and 'grey').*

Careful analysis would be needed to determine whether this is a genuine rearrangement relative to the rest of the genome, as these are draft genomes and contig orientation can either be forward or reverse. In this case, you can see that genes in bin_5 are in reversed order relative to the other bin contigs, hence, we can manually rotate the contig.

##### Rotate the contig and update the annotation

Rotate the orientation of the contig from bin_5:

```R
blast1 = mutate(blast1, direction = 1)
blast2 = mutate(blast2, direction = 1)

#annot3 <- annotation(x1 = c(bin_4_ds$start[1], bin_4_ds$start[2],bin_4_ds$start[3],bin_4_ds$start[4],bin_4_ds$start[7]),
#                    x2 = c(bin_4_ds$end[1], bin_4_ds$end[2], bin_4_ds$end[3], bin_4_ds$end[6], bin_4_ds$end[7]),
#                    text = c("cysA","cysW","cysU","unknown domain" ,"sbp"),
#                    rot = 0, col = "black")

#edit the color
bin_5_ds$col = "blue"
bin_5_ds$fill = "blue"
```

Regenerate the figure:

```R
#AFTER ROTATING
#pdf("genoplot_rotated.pdf",colormodel = "cmyk",width = 8, height = 4, paper = 'special')
plot_gene_map(dna_segs = list(bin_4_ds,bin_5_ds,bin_7_ds), 
              gene_type = "arrows", dna_seg_labels = c("bin_4", "bin_5","bin_7"), 
              comparisons = list(blast1,blast2), dna_seg_scale = TRUE, 
              override_color_schemes=TRUE,annotations=list(annot1,NULL,annot2),
              annotation_height=1.7, dna_seg_label_cex = 1,xlims=list(NULL, c(Inf, -Inf), NULL),main = "Sulfur assimilation")
#dev.off()
```


![png](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/MGSS2020_DEV/materials/figures/ex15_gene_synteny_fig2.png)


All done! We can see here that compared to bin_5 and bin_7 the following differences are apparent in the gene operon for bin_4:

1. Three genes with unknown domain/function present in bin_4 in between *cysU* and *SBP*
2. *cysW* and *cysA* are highly conserved among the genomes and have higher similarity compared to *cysU* and *SBP*.

---
