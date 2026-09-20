# Modules and R libraries needed

!!! note
    Last updated: 21 Sep 2026

## Modules in material
Comment next to modules indicates if there is a later version currently on NESI

```
module load FastQC/0.12.1
module load MultiQC/1.13-gimkl-2022a-Python-3.10.5 	#nesi latest MultiQC/1.24.1-foss-2023a-Python-3.11.6
module load Trimmomatic/0.39-Java-1.8.0_144
module load BBMap/39.01-GCC-11.3.0
module load SPAdes/4.0.0-foss-2023a-Python-3.11.6
module load seqmagick/0.8.4-gimkl-2020a-Python-3.8.2
module load Miniconda3/23.10.0-1
module load USEARCH/11.0.667-i86linux32 	#nesi latest USEARCH/11.0.667-i86linux64
module load Kraken2/2.1.3-GCC-11.3.0  	#nesi latest Kraken2/2.1.6-GCC-12.3.0
module load Bracken/2.7-GCC-11.3.0
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0 	#nesi latest SAMtools/1.23.1-GCC-12.3.0
module load MetaBAT/2.15-GCC-11.3.0 	#nesi latest MetaBAT/2.17-GCC-12.3.0
module load MaxBin/2.2.7-GCC-11.3.0-Perl-5.34.1
module load DAS_Tool/1.1.5-gimkl-2022a-R-4.2.1
module load CheckM/1.2.3-foss-2023a-Python-3.11.6
module load CheckM2/1.0.1-Miniconda3
module load GTDB-Tk/2.4.0-foss-2023a-Python-3.11.6 	#nesi latest GTDB-Tk/2.7.1-foss-2023a-Python-3.11.6
module load FastTree/2.1.11-GCCcore-9.2.0
module load prodigal/2.6.3-GCCcore-7.4.0
module load Metaxa2/2.2.3-gimkl-2022a
module load DRAM/1.3.5-Miniconda3
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5 # different version used in drep section to earlier
module load drep/2.3.2-gimkl-2018b-Python-3.7.3 	#nesi latest drep/3.4.2-gimkl-2022a-Python-3.10.5
module load seqtk/1.4-GCC-11.3.0 	#nesi latest seqtk/1.5-GCC-12.3.0
module load BLAST/2.16.0-GCC-12.3.0 	#nesi latest BLAST/2.17.0-GCC-15.2.0
```



## Optional modules: these are used only in supp material or in dropdowns/optional sections in main material

```
module load Python/3.10.5-gimkl-2022a 	#nesi latest Python/3.14.4-foss-2026
module load QUAST/5.0.2-gimkl-2018b 	#nesi latest QUAST/5.2.0-gimkl-2022a
module load IDBA-UD/1.1.3-GCC-11.3.0
module load pigz/2.7
module load DIAMOND/2.1.6-GCC-11.3.0 	#nesi latest DIAMOND/2.2.1-GCC-15.2.0
module load HMMER/3.3.2-GCC-11.3.0 	#nesi latest HMMER/3.4-GCC-15.2.0
module load SignalP/6.0g-gimkl-2022a-Python-3.10.5
module load CONCOCT/1.0.0-gimkl-2018b-Python-2.7.16 	#nesi latest CONCOCT/1.1.0-gimkl-2020a-Python-3.8.2
module load Bowtie2/2.4.5-GCC-11.3.0 	# different version to earlier in main material
module load Python/3.8.2-gimkl-2020a 	# different version again to earlier
module load SAMtools/1.15.1-GCC-11.3.0 	# different version to earlier in main material
module load R/4.2.1-gimkl-2022a 	#nesi latest R/4.6.0-foss-2026
module load prodigal-gv/2.9.0-GCC-11.3.0 	#nesi latest prodigal-gv/2.11.0-GCC-15.2.0
module load Apptainer/1.2.5  	#nesi latest Apptainer/1.3.1
module load MCL/14.137-gimkl-2020a  	#nesi latest MCL/22.282
module load DIAMOND/2.1.9-GCC-11.3.0 	# different version again to earlier
module load Python/3.11.6-foss-2023a 	# third different version 
```

## R libraries needed

```r
install.packages('ade4')
install.packages('genoPlotR')
install.packages('pheatmap')
install.packages('gplots')
install.packages('vegan')
install.packages('tidyverse')

if (!require(BiocManager)) {
  install.packages("BiocManager")
  BiocManager::install("pathview", update = FALSE)
}

if (!require(BiocManager)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST", update = FALSE)
}
```

## Optional R libraries 
```r
install.packages("scales")
```