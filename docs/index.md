<center>
![image](./theme_images/uoa_ga_nesi_LOGOS.png){width="350"}
</center>
<style>h1 {text-align: center;}</style>
<h1><b>Metagenomics Summer School</b></h1>




- - - 

|<div style="width:350px"> **Day**</div>                                         | **Lesson overview**                           | 
|:---------------------------------------------------|:---------------------------------------|
|[  Day 1  ](""){ .md-button .md-button--primary }   | [1. Intro Session I (pre-summer school): Shell ](./day1/ex1_bash_and_scheduler.md)<br>[2. Intro Session II. : HPC and HPC Job Scheduler](./day1/ex2_1_intro_to_scheduler.md)</br><br>[3. Quality filtering raw reads](./day1/ex2_quality_filtering.md)</br><br>[4. Assembly](./day1/ex3_assembly.md)</br><br>[5. Assembly (part 2)](./day1/ex4_assembly.md)</br><br>[6. Evaluating the assemblies](./day1/ex5_evaluating_assemblies.md)</br>                                       |
|[  Day 2  ](""){ .md-button .md-button--primary }   | [1. Introduction to binning](./day2/ex6_initial_binning.md)<br>[2. Binning (continued)](./day2/ex7_initial_binning.md)</br><br>[3. Bin dereplication](./day2/ex8_bin_dereplication.md)</br><br>[4. Manually refining bins](./day2/ex9_refining_bins.md)</br>                                        |
|[  Day 3  ](""){ .md-button .md-button--primary }   |[1. Identifying viral contigs in metagenomic data](./day3/ex10_viruses.md)<br>[2. Assigning taxonomy to refined prokaryotic bins](./day3/ex11_coverage_and_taxonomy.md)</br><br>[3. Gene prediction](./day3/ex12_gene_prediction.md)</br><br>[4. Gene annotation (part 1)](./day3/ex12_gene_prediction.md)</br><br>[5. Gene annotation (part 2) and coverage calculation](./day3/ex14_gene_annotation_part2.md)</br>                                     |
|[  Day 4  ](""){ .md-button .md-button--primary }   |[1. Gene annotation (part 3)](./day4/ex15_gene_annotation_part3.md)<br>[2. Presentation of data - Intro](./day4/ex16a_data_presentation_Intro.md)</br><br>[3. Presentation of data: Per-sample coverage heatmaps](./day4/ex16b_data_presentation_Coverage.md)</br><br>[4. Presentation of data: Ordinations](./day4/ex16c_OPTIONAL_data_presentation_Ordination.md)</br><br>[5. Presentation of data: KEGG pathway maps](./day4/ex16d_data_presentation_KEGG_pathways.md)</br><br>[6. Presentation of data: Gene synteny](./day4/ex16e_data_presentation_Gene_synteny.md)</br><br>[7. Presentation of data: CAZY annotations heatmap](./day4/ex16f_OPTIONAL_data_presentation_CAZy_annotations.md)</br>                                        |

<br>

??? calendar-days "Timetable 2022"
    

    === "Day 1 : 29<sup>th</sup> Nov"

        |Time|Event|Session leader|
        |:---|:---|:---|
        |9:00 am – 9:30 am|**Introduction**<br>- Welcome<br>- Logging into NeSI|Jian Sheng Boey|
        |9:30 am – 10:30 am|**TASK:** Bash scripting and introduction to Slurm scheduler|Dinindu Senanayake|
        |10:30 am – 10:50 am|**Morning tea break**||
        |10:50 am – 11:10 am|**TASK:** Bash scripting (continued)|Dinindu Senanayake|
        |11:10 am – 11:40 pm|**TALK:** The metagenomics decision tree|Kim Handley|
        |11:40 am – 12:00 pm|**DISCUSSION** Guided discussions on approaches to metagenomics analyses|Kim Handley|
        |12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
        |12:45 pm – 1:45 pm|**TALK:** Quality filtering raw reads<br>**TASK:** Visualisation with *FastQC*<br>**TASK:** Read trimming and adapter removal<br>**TASK:** Diagnosing poor libraries<br>**TALK:** Common issues and best practice<br>**TASK (*Optional*):** Filtering out host DNA|Annie West|
        |1:45 pm – 3:00 pm|**TALK:** Assembly<br>- Choice of assemblers<br>- Considerations for parameters, and when to stop!<br>**TASK:** Prepare data for assembly<br>**TASK:** Exploring assembler options<br>**TASK:** Submitting jobs to NeSI via slurm<br>**TASK:** Run SPAdes and IDBA-UD assembly<br>**TASK (*Optional*):** Submitting variant assemblies to NeSI|Kim Handley|
        |3:00 pm – 3:20 pm|**Afternoon tea break** (*Tea, coffee, and snacks provided*)||
        |3:20 pm – 5:00 pm|**TALK:** Future considerations - co-assembly vs. single assemblies<br>**TASK:** Assembly evaluation<br>**TASK:** Short contig removal|Jian Sheng Boey<br> Kim Handley|


    === "Day 2 : 30<sup>th</sup> Nov"

        |Time|Event|Session leader|
        |:---|:---|:---|
        |9:00 am – 9:15 am|**Introduction**<br>- Overview of yesterday, questions<br>- Overview of today|Annie West|
        |9:15 am – 9:30 am|**RECORDED PRESENTATION:** Genomic adaptations enabling *Acidithiobacillus* distribution across wide ranging hot spring temperatures and pHs (Chanenath Sriaporn)|<br> Kim Handley|
        |9:30 am – 10:30 am|**Binning (part 1)**<br>**TALK:** Overview of binning history<br>- Key parameters and strategies for binning<br>**TASK:** Read mapping|Kim Handley|
        |10:30 am – 10:50 am|**Morning tea break** (*Tea, coffee, and snacks provided*)||
        |10:50 am – 11:20 am|**Binning (part 2)**<br>**TASK:** Multi-binning strategy (*Metabat* and *Maxbin*)|Cesar Facimoto|
        |11:20 am – 12:00 pm|**TALK:** Overview of binning history (*continued*)<br>- Key parameters and strategies for binning|Kim Handley|
        |12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
        |12:45 pm – 2:00 pm|**Binning (part 3)**<br>**TASK:** Bin dereplication via *DAS_Tool*<br>**TASK:** Evaluating bins using *CheckM*|Cesar Facimoto|
        |2:00 pm - 2:30 pm|**Binning (part 4)**<br>**TALK:** Discuss additional dereplication strategies, such as *dRep*<br>**TALK:** How to work with viral and eukaryotic bins<br>**TALK:** Dealing with organisms which possess minimal genomes<br>|Jian Sheng Boey<br> Annie West|
        |2:30 pm - 3:00 pm|**TALK:** Bin refinement<br>- Refinement strategies| |
        |3:00 pm – 3:20 pm|**Afternoon tea break** (*Tea, coffee, and snacks provided*)||
        |3:20 pm – 5:00 pm|**TASK:** Working with *VizBin*<br>|Jian Sheng Boey<br> Annie West|

    === "Day 3 : 1<sup>st</sup> Dec"

        |Time|Event|Session leader|
        |:---|:---|:---|
        |9:00 am – 9:20 am|**Introduction**<br>- Overview of yesterday, questions<br>- Overview of today|Jian Sheng Boey|
        |9:20 am – 10:30 am|**TALK:** Bin taxonomic classification<br>- Bin and species determination<br>**TASK:** Taxonomic classification using GTDB-Tk<br>**TASK:** View phylogenetic trait distribution using *ANNOTREE*|David Waite|
        |10:30 am – 10:50 am|**Morning tea break** (*Tea, coffee, and snacks provided*)||
        |10:50 am – 12:00 pm|**TALK:** Identifying viruses from metagenomic data<br>**TASK:** Identifying viral contigs using *VIBRANT*<br>**RECORDED PRESENTATION:** Genetic exchange in ultra-small Patescibacteria (Emilie Gios)<br>**TASK:** QC of viral contigs using *CheckV*<br>**TASK:** Taxonomic classification of viruses using *vContact2*|Jian Sheng Boey<br>Annie West<br><br>|
        |12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
        |12:45 pm - 1:30 pm|**TALK:** Gene prediction using *prodigal* and other tools (*RNAmer*, *Aragorn*, etc)<br>**TASK:** Predict open reading frames and protein sequences|Jian Sheng Boey <br> Cesar Facimoto|
        |1:30 pm - 2:00 pm|**TALK:** Gene annotation (part 1)<br> - Methods<br>**TASK:** Gene annotation using *DIAMOND* and *HMMER3*<br>**Discussion:** Evaluating the quality of gene assignment<br>|Jian Sheng Boey<br>Cesar Facimoto|
        |2:00 pm – 3:00 pm|**TALK:** Gene annotation (part 2)<br>- Using online resources (e.g. *KEGG, BioCyc, MetaCyc, HydDB, PSORT*)<br>**TASK:** View KEGG annotation in KEGG website|Jian Sheng Boey<br>Cesar Facimoto|
        |3:00 pm – 3:20 pm|**Afternoon tea break**||
        |3:20 pm – 4:30 pm|**TASK:** MAG annotation with *DRAM*<br>**TASK:** Coverage calculation using *Bowtie2*<br>**TASK:** Introduce group project goals<br>**TASK:** Dividing into working groups / get a group name<br>**TASK:** Select a goal from your project|Annie West<br>Kim Handley|
        |4:30 pm – 5:00 pm|**End of day wrap up**|Kim Handley|


    === "Day 4 : 2<sup>nd</sup> Dec"

        |Time|Event|Session leader|
        |:---|:---|:---|
        |9:00 am – 9:15 am|**Introduction**<br>- Overview of yesterday, questions<br>- Overview of today|Jian Sheng Boey|
        |9:15 am – 10:00 am|**Presentation of data**<br>**TALK:** *DRAM* results overview<br>**TASK:** Explore *DRAM* results|Annie West|
        |10:00 am – 10:30 am|**Presentation of data (genome distributions)**<br>**TALK:** Visualising findings (environmental distribution)<br>**TASK:** Coverage heatmap and nMDS ordination |Kim Handley|
        |10:30 am – 10:50 am|**Morning tea break** (*Tea, coffee, and snacks provided*)<br>**TASK:** Workshop survey||
        |10:50 am – 12:00 pm|**Presentation of data (metabolism)**<br>**TALK:** Visualising findings (metabolic maps, genome features, metabolic schematics, and gene trees)<br>**TASK:** KEGG metabolic pathways<br>**TASK:** Gene synteny<br>**TASK:** CAZy heatmaps|Cesar Facimoto<br>Jian Sheng Boey|
        |12:00 pm – 12:45 pm|**Break for lunch** (*lunch not provided*)||
        |12:45 pm – 2:30 pm|**TASK:** Analyse data for group work<br>**TASK:** Prepare group presentation|Kim Handley|
        |2:30 pm – 3:00 pm|**Present and discuss findings**<br>**TASK:** Each group to give an informal presentation of their data|Kim Handley|
        |3:00 pm – 3:20 pm|**Afternoon tea break** (*Tea, coffee, and snacks provided*)||
        |3:20 pm – 3:40 pm|**Present and discuss findings (continued)**<br>**TASK:** Each group to give an informal presentation of their data|Annie West|
        |3:40 pm – 4:00 pm|**End of day wrap up**<br>- Final discussion|Kim Handley<br>Jian Sheng Boey|

<br>


??? example "Appendices"

    |Appendix ID                                         |                        | 
    |:---------------------------------------------------|:---------------------------------------|
    | [Appendix 1](""){ .md-button } | [Dereplicating data from multiple assemblies](./resources/1_APPENDIX_ex8_Dereplication.md)|
    |[Appendix 2](""){ .md-button  } | [Generating input files for "VizBin" from "DAS_Tool" curated bins](./resources/2_APPENDIX_ex9_Generating_input_files_for_VizBin.md) |
    |[Appendix 3](""){ .md-button    }| [Normalise per-sample coverage values by average library size](./resources/3_APPENDIX_ex11_Normalise_coverage_example.md) |
    |[Appendix 4](""){ .md-button    } | [Viral taxonomy prediction via vContact2](./resources/4_APPENDIX_ex11_viral_taxonomy_prediction_via_vContact2.md) |
    |[Appendix 5](""){ .md-button } | [Preparing input files for visualising gene synteny](./resources/5_APPENDIX_ex15_Prepare_gene_synteny_inputs.md) |
<br>

!!! comment-dots "Post-workshop survey"
    
    Thank you for attending Metagenomics Summer School 2022! We would like your feedback on how we have done and what we can improve on. You can provide feedback [here](https://auckland.au1.qualtrics.com/jfe/form/SV_3W4gkA3XZ0qK2do).

<br>

!!! note "License"

    Genomics Aotearoa / New Zealand eScience Infrastructure/ University of Auckland **Metagenomics Summer School** material is licensed under the **GNU General Public License v3.0, 29 June 2007** . ([Follow this link for more information](https://github.com/GenomicsAotearoa/shell-for-bioinformatics/blob/main/LICENSE))

<br>

!!! note "Slides for workshop"

    You can find a copy of the slides presented during the workshop, with published figures removed, in the [slides/](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/slides) folder.

<br>

!!! desktop-download-24 "Snapshots of results to download"

    If you are having trouble downloading files using `scp`, we are providing exemplar output files, which you can download through your browser, [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/docs/resources), or via the following links: 
    
    * [FastQC results](./resources/fastqc_results.zip)
    * [Quast results](./resources/quast_results_sans_reference.zip) and [required references](./resources/quast_references.zip)
    * [Input files for VizBin](./resources/vizbin_files.zip)
    * [Gene annotation tables](./resources/example_annotation_tables.zip)
    * [DRAM output files](./resources/dram_distillation.zip)



