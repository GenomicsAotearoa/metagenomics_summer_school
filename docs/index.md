<center>
![image](./theme_images/uoa_ga_nesi_LOGOS.png){width="350"}
</center>

# <span style="display: block; text-align: center;">**Metagenomics Summer School**</span>

- - - 

| <div style="width:150px">Day</div> | <div style="width:400px">Lesson topic</div> |
| --- | --- |
| [Day 1: Assembly](#){ .md-button .md-button--primary } | [1. (Pre-Summer School) Introduction I: Shell](./day1/ex1_bash_and_scheduler.md)<br>[2. (Pre-Summer School) Introduction II: HPC and job scheduler](./day1/ex2_1_intro_to_scheduler.md)<br>[3. Filter raw reads by quality](./day1/ex2_quality_filtering.md)<br>[4. Assembly I: Assembling contigs](./day1/ex3_assembly.md)<br>[5. Assembly II: Variable parameters](./day1/ex4_assembly.md)<br>[6. Assembly evaluation](./day1/ex5_evaluating_assemblies.md) |
| [Day 2: Binning](#){ .md-button .md-button--primary } | [1. Introduction to binning](./day2/ex6_initial_binning.md)<br>[2. Binning with multiple tools](./day2/ex7_initial_binning.md)<br>[3. Bin dereplication](./day2/ex8_bin_dereplication.md)<br>[4. Manual bin refinement](./day2/ex9_refining_bins.md)<br>[3. Identifying viral contigs in metagenomic data](./day2/ex10.1_viruses.md) |
| [Day 3: Annotation](#){ .md-button .md-button--primary } | [1. Assigning taxonomy to refined prokaryotic bins](./day3/ex11_coverage_and_taxonomy.md)<br>[2. Phylogenomics](./day3/ex11.1_phylogenomics.md)<br>[3. Virus taxonomy](./day3/ex10.2_viruses.md)<br>[4. Gene prediction](./day3/ex12_gene_prediction.md)<br>[5. Gene annotation I: BLAST-like and HMM](./day3/ex13_gene_annotation_part1.md)<br>[6. Gene annotation II: DRAM and coverage calculation](./day3/ex14_gene_annotation_part2.md) |
| [Day 4: Visualisation](#){ .md-button .md-button--primary } | [1. Gene annotation III: DRAM distillation](./day4/ex15_gene_annotation_part3.md)<br>[2. Introduction to data presentation](./day4/ex16a_data_presentation_Intro.md)<br>[3. Coverage heatmaps](./day4/ex16b_data_presentation_Coverage.md)<br>[4. Ordinations](./day4/ex16c_OPTIONAL_data_presentation_Ordination.md)<br>[5. KEGG pathway maps](./day4/ex16d_data_presentation_KEGG_pathways.md)<br>[6. Gene synteny](./day4/ex16e_data_presentation_Gene_synteny.md)<br>[7. CAZy heatmaps](./day4/ex16f_OPTIONAL_data_presentation_CAZy_annotations.md) |

<br>

<style>h1 {text-align: left;}</style>
## :material-calendar: Timetable 2024

=== "Day 1: Tuesday, 3<sup>rd</sup> Sep"

    | <div style="width:100px">Time</div> | <div style="width:350px">Event</div> | <div style="width:125px">Session leader</div> |
    | --- | --- | --- |
    | 09:00 &ndash; 09:25 | **Introductions**<br>Welcome<br>Overview <br>Login to NeSI via [Jupyter Hub](https://jupyter.nesi.org.nz/) and test script | Jian Sheng Boey |
    | 09:25 &ndash; 09:55 | :fontawesome-solid-microphone-lines: **TALK** Metagenomics decision tree<br>:octicons-comment-discussion-16: **DISCUSSION** What are your research questions? | Kim Handley |
    | 09:55 &ndash; 10:30 | **Read QC**<br>:fontawesome-solid-microphone-lines: **TALK** Quality filtering raw reads<br>:fontawesome-solid-laptop-code: **TASK** Visualisation with *FastQC*  and *MultiQC* | Annie West |
    | 10:30 &ndash; 10:50 | :material-tea: *Morning tea* | |
    | 10:50 &ndash; 11:25 | :fontawesome-solid-laptop-code: **TASK** Read trimming & adapter removal<br>:fontawesome-solid-laptop-code: **TASK** Diagnosing poor libraries<br>:fontawesome-solid-microphone-lines: **TALK** Common issues & best practices<br>:fontawesome-solid-laptop-code: **TASK** Filter out host DNA | Annie West |
    | 11:25 &ndash; 12:30 | **Assembly**<br>:fontawesome-solid-microphone-lines: **TALK** Choice of assemblers, parameter considerations, and when to stop!<br>:fontawesome-solid-laptop-code: **TASK** Prepare data for assembly<br>:fontawesome-solid-laptop-code: **TASK** Exploring assembler options<br>:fontawesome-solid-laptop-code: **TASK** Submit assembly jobs to NeSI via Slurm<br>:fontawesome-solid-laptop-code: **TASK** Submit variant assembly jobs | Kim Handley |
    | 12:30 &ndash; 13:30 | :material-silverware-fork-knife: *Lunch* | |
    | 13:30 &ndash; 14:15 | :fontawesome-solid-microphone-lines: **TALK** Other considerations: assembly strategies<br>:fontawesome-solid-microphone-lines: **TALK** Assembly evaluation<br>:fontawesome-solid-laptop-code: **TASK** Evaluate your assemblies<br>:fontawesome-solid-laptop-code: **TASK** Short contig removal | Mike Hoggard |
    | 14:15 &ndash; 15:00 | :fontawesome-solid-microphone-lines: **TALK** Other considerations: rRNA reconstruction<br>:fontawesome-solid-laptop-code: **TASK** Ribosomal RNA reconstruction using *PhyloFlash*<br>:fontawesome-solid-microphone-lines: **TALK** Other considerations: read classification | Annie West<br>Jian Sheng Boey |
    | 15:00 &ndash; 15:20 | :material-tea: *Afternoon tea* | |
    | 15:20 &ndash; 15:50 | :fontawesome-solid-laptop-code: **TASK** Sequence classification using *Kraken2* and *Bracken* | Annie West<br>Jian Sheng Boey |
    | 15:50 &ndash; 16:00 | :fontawesome-solid-microphone-lines: **TALK** Project introduction and description | Kim Handley |

=== "Day 2: Wednesday, 4<sup>th</sup> Sep"

    | <div style="width:100px">Time</div> | <div style="width:350px">Event</div> | <div style="width:125px">Session leader</div> |
    | --- | --- | --- |
    | 09:00 &ndash; 09:10 | **Introductions**<br>Recap of day 1<br>Overview of the day | Annie West |
    | 09:10 &ndash; 10:30 | :fontawesome-solid-microphone-lines: **TALK** Overview of binning history, key parameters<br>:fontawesome-solid-laptop-code: **TASK** Read mapping<br>:fontawesome-solid-laptop-code: **TASK** Multi-binning strategy with MetaBat and MaxBin | Kim Handley |
    | 10:30 &ndash; 10:50 | :material-tea: *Morning tea* | |
    | 10:50 &ndash; 11:45 | :fontawesome-solid-microphone-lines: **TALK** Binning strategies<br>:fontawesome-solid-laptop-code: **TASK** Bin dereplication using *DAS_Tool* | Kim Handley |
    | 11:45 &ndash; 12:10 | :material-presentation-play: **PRESENTATION** Alternative approaches to binning | Amali Thrimawithana |
    | 12:10 &ndash; 13:00 | :material-silverware-fork-knife: *Lunch* | |
    | 13:00 &ndash; 15:00 | :fontawesome-solid-laptop-code: **TASK** Bin evaluation using *CheckM* and *CheckM2*<br>:fontawesome-solid-microphone-lines: **TALK** Bin dereplication across assemblies, working with minimal or non-prokaryotic genomes, and bin refinement strategies<br>:fontawesome-solid-laptop-code: **TASK** Visual bin refinement using *VizBin* | Kim Handley<br>Mike Hoggard |
    | 15:00 &ndash; 15:20 | :material-tea: *Afternoon tea* | |
    | 15:20 &ndash; 15:35 | :fontawesome-solid-laptop-code: **TASK** Checking VizBin results | Annie West |
    | 15:35 &ndash; 16:05 | :fontawesome-solid-microphone-lines: **TALK** Viruses in metagenomic data<br>:fontawesome-solid-laptop-code: **TASK** Identify viral contigs using *VirSorter2* and *CheckV* | Mike Hoggard |
    | 16:05 &ndash; 16:30 | :octicons-comment-discussion-16: **DISCUSSION** How will you approach your research question? | |

=== "Day 3: Thursday, 5<sup>th</sup> Sep"

    | <div style="width:100px">Time</div> | <div style="width:350px">Event</div> | <div style="width:125px">Session leader</div> |
    | --- | --- | --- |
    | 09:00 &ndash; 09:10 | **Introductions**<br>Recap of day 2<br>Recap of binning outputs<br>Overview of the day | Jian Sheng Boey |
    | 09:10 &ndash; 10:20 | :fontawesome-solid-microphone-lines: **TALK** Taxonomic classification: Bin and species determination<br>:fontawesome-solid-laptop-code: **TASK** Taxonomic classification using GTDB-Tk<br>:fontawesome-solid-laptop-code: **TASK** View phylogenetic trait distribution using ANNOTREE | Annie West |
    | 10:20 &ndash; 10:40 | :material-tea: *Morning tea* | |
    | 10:40 &ndash; 11:15 | :fontawesome-solid-microphone-lines: **TALK** Primer to phylogenetic analyses<br>:fontawesome-solid-laptop-code: **TASK** Build a phylogenetic tree using FastTree and visualise tree using iTOL | Jian Sheng Boey |
    | 11:15 &ndash; 12:05 | :fontawesome-solid-laptop-code: **TASK** Exploring results from VirSorter2 and CheckV<br>:fontawesome-solid-microphone-lines: **TALK** Predicting viral taxonomy using vConTACT2<br>:fontawesome-solid-laptop-code: **TASK (Optional)** Visualise vConTACT2 viral gene sharing network in Cytoscape | Mike Hoggard |
    | 12:05 &ndash; 12:50 | :material-silverware-fork-knife: *Lunch* | |
    | 12:50 &ndash; 13:20 | :fontawesome-solid-microphone-lines: **TALK** Tools for predicting genes: Prodigal, RNAmer, Aragorn, etc.<br>:fontawesome-solid-laptop-code: **TASK** Predict ORFs and protein sequences using Prodigal | Jian Sheng Boey |
    | 13:20 &ndash; 14:20 | :fontawesome-solid-microphone-lines: **TALK** Methods for gene annotation<br>:fontawesome-solid-laptop-code: **TASK** Annotate predicted genes using DIAMOND and HMMER3 | Jian Sheng Boey |
    | 14:20 &ndash; 14:50 | :fontawesome-solid-microphone-lines: **TALK** Using online resources (KEGG, BioCyc, MetaCyc, InterPro, PSORT)<br>:fontawesome-solid-laptop-code: **TASK** View KEGG annotations online | Jian Sheng Boey |
    | 14:50 &ndash; 15:05 | :material-tea: *Afternoon tea* | |
    | 15:05 &ndash; 15:30 | :material-presentation-play: **PRESENTATION** Metagenomic analyses *sans* binning | Jess Wallbank |
    | 15:30 &ndash; 15:50 | :fontawesome-solid-laptop-code: **TASK** MAG annotation with DRAM<br>:fontawesome-solid-laptop-code: **TASK** Coverage calculation | Annie West |
    | 15:50 &ndash; 16:00 | :fontawesome-solid-laptop-code: **TASK** Group formation and project topic selection | Kim Handley |

=== "Day 4: Friday, 6<sup>th</sup> Sep"

    | <div style="width:100px">Time</div> | <div style="width:350px">Event</div> | <div style="width:125px">Session leader</div> |
    | --- | --- | --- |
    | 09:00 &ndash; 09:10 | **Introductions**<br>Recap of day 3<br>Overview of the day | Jian Sheng Boey |
    | 09:10 &ndash; 09:55 | :fontawesome-solid-microphone-lines: **TALK** Overview of DRAM results<br>:fontawesome-solid-laptop-code: **TASK** Explore DRAM results | Annie West |
    | 09:55 &ndash; 10:25 | :fontawesome-solid-microphone-lines: **TALK** Visualising environmental distribution<br>:fontawesome-solid-laptop-code: **TASK** Coverage heatmap and nMDS ordination | Mike Hoggard |
    | 10:25 &ndash; 10:45 | :material-tea: *Morning tea* | |
    | 10:45 &ndash; 10:55 | :fontawesome-solid-laptop-code: **TASK** Workshop survey | |
    | 10:55 &ndash; 12:05 | :fontawesome-solid-microphone-lines: **TALK** Visualising genomic and metabolic features<br>:fontawesome-solid-laptop-code: **TASK** KEGG metabolic maps<br>:fontawesome-solid-laptop-code: **TASK** Gene synteny | Jian Sheng Boey<br>Annie West |
    | 12:05 &ndash; 12:55 | :material-silverware-fork-knife: *Lunch* | |
    | 12:45 &ndash; 14:40 | **Group work**<br>:fontawesome-solid-laptop-code: **TASK** Analyse data and prepare presentations | Kim Handley<br>Mike Hoggard<br>Annie West<br>Jian Sheng Boey |
    | 14:40 &ndash; 14:55 | :material-presentation-play: **PRESENTATION** Group presentations | Kim Handley<br>Mike Hoggard<br>Annie West<br>Jian Sheng Boey |
    | 14:55 &ndash; 15:10 | :material-tea: *Afternoon tea* | |
    | 15:10 &ndash; 15:55 | :material-presentation-play: **PRESENTATION** Group presentations | Kim Handley<br>Mike Hoggard<br>Annie West<br>Jian Sheng Boey |
    | 15:55 &ndash; 16:00 | :material-hand-wave: **Wrap up and final discussions** | Jian Sheng Boey |

<br>

??? example "Appendices"

    |Appendix ID                                         |                        | 
    |:---------------------------------------------------|:---------------------------------------|
    | [Appendix 1](""){ .md-button } | [Dereplicating data from multiple assemblies](./resources/1_APPENDIX_ex8_Dereplication.md)|
    |[Appendix 2](""){ .md-button  } | [Generating input files for "VizBin" from "DAS_Tool" curated bins](./resources/2_APPENDIX_ex9_Generating_input_files_for_VizBin.md) |
    |[Appendix 3](""){ .md-button    }| [Normalise per-sample coverage values by average library size](./resources/3_APPENDIX_ex11_Normalise_coverage_example.md) |
    |[Appendix 4](""){ .md-button    } | [Viral taxonomy prediction via vContact2](./resources/4_APPENDIX_ex11_viral_taxonomy_prediction_via_vContact2.md) |
<!--    |[Appendix 5](""){ .md-button } | [Preparing input files for visualising gene synteny](./resources/5_APPENDIX_ex15_Prepare_gene_synteny_inputs.md) | -->

<br>

!!! comment-dots "Post-workshop survey"
    
    Thank you for attending Metagenomics Summer School 2024! We would like your feedback on how we have done and what we can improve on. You can provide feedback [here](https://auckland.au1.qualtrics.com/jfe/form/SV_4JcEoG7GxaugScC).

<br>

!!! copyright "License"

    Genomics Aotearoa / New Zealand eScience Infrastructure / University of Auckland **Metagenomics Summer School** material is licensed under the **GNU General Public License v3.0, 29 June 2007**. ([Follow this link for more information](https://github.com/GenomicsAotearoa/shell-for-bioinformatics/blob/main/LICENSE))

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
