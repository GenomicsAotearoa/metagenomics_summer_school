# Bin coverage 2: Coverage calculation


## Calculating per-sample coverage stats for prokaryotic bins

One of the first questions we often ask when studying the ecology of a system is: What are the pattens of abundance and distribution of taxa across the different samples? With bins of metagenome-assembled genome (MAG) data, we can investigate this by mapping the quality-filtered unassembled reads back to the refined bins to then generate coverage profiles. Genomes in higher abundance in a sample will contribute more genomic sequence to the metagenome, and so the average depth of sequencing coverage for each of the different genomes provides a proxy for abundance in each sample. 

We can do this using the tool `CoverM`. This exercise will take place in the `10.gene_annotation_and_coverage/` folder which contains our results on the mock metagenome data. `CoverM` will retrieve the per-contig coverages from the *.bam* files calculated earlier and calculate a mean coverage for each MAG (`dastool_bins`). 

!!! terminal-2 "Running `CoverM` to get per-sample coverage of each bin as a whole"

    ```bash
    # Load CoverM
    module purge
    module load CoverM/0.7.0-GCC-12.3.0

    # Run CoverM
    coverm genome --min-covered-fraction 0 \
        --methods mean variance covered_fraction relative_abundance \
        --genome-fasta-directory dastool_bins \
        --bam-files ../5.binning/*.bam \
        --output-file coverm_mags.tsv
    ```

The coverage table will be generated as coverm_mags.tsv. The key columns of interest are the `sampleX Mean` (coverage), `sampleX Covered Fraction`,  and `sampleX Relative Abundance (%)`. For a full explanation of all available statistical methods cf. [https://github.com/wwood/CoverM#calculation-methods](https://github.com/wwood/CoverM#calculation-methods)

If we were interested in the individual per-contig coverages (and not their per-MAG averages), we would use `CoverM` as follows:

!!! terminal-2 "Running `CoverM` to get per-contig coverages"

    ```bash
    # Load CoverM
    module purge
    module load CoverM/0.7.0-GCC-12.3.0

    # Run CoverM
    coverm contig \
        --methods mean variance \
        --bam-files ../5.binning/*.bam \
        --output-file coverm_contigs.tsv
    ```

---
