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

!!! note "Timetable 2022"
    The timetable of daily events can be found [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/Timetable%20for%20MGSS%202022.md)

??? example "Appendices"

    |Appendix ID                                         |                        | 
    |:---------------------------------------------------|:---------------------------------------|
    | [Appendix 1](""){ .md-button } | [Dereplicating data from multiple assemblies](./resources/1_APPENDIX_ex8_Dereplication.md)|
    |[Appendix 2](""){ .md-button  } | [Generating input files for "VizBin" from "DAS_Tool" curated bins](./resources/2_APPENDIX_ex9_Generating_input_files_for_VizBin.md) |
    |[Appendix 3](""){ .md-button    }| [Normalise per-sample coverage values by average library size](./resources/3_APPENDIX_ex11_Normalise_coverage_example.md) |
    |[Appendix 4](""){ .md-button    } | [Viral taxonomy prediction via vContact2](./resources/4_APPENDIX_ex11_viral_taxonomy_prediction_via_vContact2.md) |
    |[Appendix 5](""){ .md-button } | [How to generate the blast files provided by Gene synteny](./resources/5_APPENDIX_ex15_gene_synteny_Generate_blast_files.md) |
    |[Appendix 6](""){ .md-button } | [Extract Gene of Interest](./resources/6_APPENDIX_ex15_gene_synteny_grab_GOI.md)|


<br>
<br>

!!! note "License"

    Genomics Aotearoa / New Zealand eScience Infrastructure/ University of Auckland **Metagenomics Summer School** material is licensed under the **GNU General Public License v3.0, 29 June 2007** . ([Follow this link for more information](https://github.com/GenomicsAotearoa/shell-for-bioinformatics/blob/main/LICENSE))

!!! note "Slides for workshop"

    You can find a copy of the slides presented during the workshop, with published figures removed, in the [slides/](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/slides) folder.


!!! desktop-download-24 "Snapshots of results to download"

    If you are having trouble downloading files using `scp`, we are providing exemplar output files, which you can download through your browser, [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/docs/resources), or via the following links: 
    
    * [FastQC results](./resources/fastqc_results.zip)
    * [Quast results](./resources/quast_results_sans_reference.zip) and [required references](./resources/quast_references.zip)
    * [Input files for VizBin](./resources/vizbin_files.zip)
    * [Gene annotation tables](./resources/example_annotation_tables.zip)
    * [DRAM output files](./resources/dram_distillation.zip)



