# Introduction to binning

### Objectives

* Remove short contigs from the data set
* Obtain coverage profiles for assembled contigs via read mapping
* Create initial bins using **MetaBAT** and **MaxBin**

---

### Remove short contigs from the data set

Ideally, we do not want to be creating bins from all of the assembled contigs, as there is often a long tail of contigs which are only several *k*-mers long. These have little biological meaning, as they are too short for robust gene annotation, and they can introduct a significant degree of noise in the clustering algorithms used for binning. We therefore identify a suitable threshold for a minimum length of contigs to be considered for binning.

We have already done this in the previous exercise](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/day2/ex5_evaluating_assemblies.md) so we could either use the existing filtering at 1,000 bp in length, or move to something stricter. Most binning tools have a default cut-off for minimum contig size - **MetaBAT** uses a default minimum of 2,500 bp, and recommends at least 1,500. By contrast, **MaxBin** sets the minimum length at 1,000 bp.

---

### Obtain coverage profiles for assembled contigs via read mapping



---

### Create initial bins using *MetaBAT* and *MaxBin*

