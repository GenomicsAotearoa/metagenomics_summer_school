<center>

# **Metagenomics Summer School**

</center>




<center>

|<div style="width:290px"> **Day**</div>                                         | **Lesson overview**                           | 
|:---------------------------------------------------|:---------------------------------------|
|[  Day 1  ](""){ .md-button .md-button--primary }   | [1. Introduction to shell and HPC Job Scheduler](./day1/ex1_bash_and_scheduler.md)<br>[2. Quality filtering raw reads](./day1/ex2_quality_filtering.md)</br><br>[3. Assembly](./day1/ex3_assembly.md)</br><br>[4. Assembly 2](./day1/ex4_assembly.md)</br><br>[5. Evaluating the assemblies](./day1/ex5_evaluating_assemblies.md)</br>                                       |
|[  Day 2  ](""){ .md-button .md-button--primary }   | [1. Introduction to binning](./day2/ex6_initial_binning.md)<br>[2. Introduction to binning pt 2](./day2/ex7_initial_binning.md)</br><br>[3. Bin dereplication](./day2/ex8_bin_dereplication.md)</br><br>[4. Manually refining bins](./day2/ex9_refining_bins.md)</br>                                        |
|[  Day 3  ](""){ .md-button .md-button--primary }   |[1. Identifying viral contrigs](./day3/ex10_viruses.md)<br>[2. Per-sample coverage and assigning taxonomy](./day3/ex11_coverage_and_taxonomy.md)</br><br>[3. Gene Prediction](./day3/ex12_gene_prediction.md)</br><br>[4. Gene annotation - Pt 1](./day3/ex12_gene_prediction.md)</br><br>[5. Gene Annotation - Pt 2](./day3/ex14_gene_annotation_part2.md)</br>                                     |
|[  Day 4  ](""){ .md-button .md-button--primary }   |[1. Gene Annotation - Pt 3](./day4/ex15_gene_annotation_part3.md)<br>[2. Presentation of data - Intro](./day4/ex16a_data_presentation_Intro.md)</br><br>[3. Presentation of data : Per-sample coverage heatmaps](./day4/ex16b_data_presentation_Coverage.md)</br><br>[4. Presentation of data :Ordinations](./day4/ex16c_OPTIONAL_data_presentation_Ordination.md)</br><br>[5. Presentation od data : KEGG pathway maps](./day4/ex16d_data_presentation_KEGG_pathways.md)</br><br>[6. Presentation of data : Gene synteny](./day4/ex16e_data_presentation_Gene_synteny.md)</br><br>[7. Presentation of data : CAZY annotations heatmap](./day4/ex16f_OPTIONAL_data_presentation_CAZy_annotations.md)</br>                                        |

</center>



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


!!! tip "Snapshots of results to download"

    If you are having trouble downloading files using `scp`, we are providing exemplar output files, which you can download through your browser, [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/tree/master/docs/resources), or via the following links: 
    
    * [FastQC results](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/docs/resources/fastqc_results.zip)
    * [Quast results](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/docs/resources/quast_results.zip)
    * [Input files for VizBin](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/docs/resources/vizbin_files.zip)
    * [Tables of gene annotations](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/docs/resources/example_annotation_tables.zip)
    * [DRAM output files](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/docs/resources/DRAM_results.zip)