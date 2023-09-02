# APPENDIX (ex15): Prepare input for gene synteny visualisation

In order to produce gene synteny plots using `genoPlotR` as outlined in [Gene synteny](../day4/ex16e_data_presentation_Gene_synteny.md), we need to know the annotations and relative nucleotide positions of our genes of interest. We can use our annotations file generated via [homology and domain searches](../day3/ex13_gene_annotation_part1.md) (here, we have parsed and aggregated these annotations using an [in-house custom script](https://github.com/GenomicsAotearoa/environmental_metagenomics/tree/master/analysis_tools/annotationaggregator_v0.1)) or [DRAM](../day3/ex14_gene_annotation_part2.md) in order to obtain relevant genes and their labels. Gene nucleotide positions relative to assembled contigs/scaffolds can be obtained from `prodigal` outputs. 

Additionally, we need BLAST outputs for comparing between genes along their contigs. For this, we rely on outputs from pairwise `tBLASTx` (translates a nucleotide database then searches it using a translated nucleotide query) to perform sequential comparisons across different bins.

For this example, we use our aggregated annotations (provided in `10.gene_annotation_and_coverage/example_annotation_tables`) and prodigal outputs to generate input files for visualising the synteny of genes that encode the [sulfate/thiosulfate ABC transporter](https://metacyc.org/META/NEW-IMAGE?type=ENZYME&object=ABC-70-CPLX) involved in [assimilating extracellular sulfate](https://www.kegg.jp/module/M00616). We will begin by subsetting our genes of interest based on KEGG orthology numbers and relevant annotation labels, followed by a pairwise `tBLASTx` of our contigs.

## 1. Navigating the working directory

We will be working in `10.gene_annotation_and_coverage/` where you will find the `example_annotation_tables/`, `filtered_bins/` and `predictions/` sub-directories already prepared for you. 

!!! terminal "code"
    
    ```bash
    # Navigate to working directory
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.gene_annotation_and_coverage

    # View relevant sub-directories
    ls example_annotation_tables/
    ```

!!! circle-check "Terminal output"

    ```
    bin_0.annotation.aa  bin_3.annotation.aa  bin_6.annotation.aa  bin_9.annotation.aa
    bin_1.annotation.aa  bin_4.annotation.aa  bin_7.annotation.aa
    bin_2.annotation.aa  bin_5.annotation.aa  bin_8.annotation.aa
    ```

!!! terminal "code"

    ```
    ls predictions/
    ```

!!! circle-check "Terminal output"

    ```
    bin_0.filtered.genes.faa              bin_5.filtered.genes.faa
    bin_0.filtered.genes.fna              bin_5.filtered.genes.fna
    bin_0.filtered.genes.gbk              bin_5.filtered.genes.gbk
    bin_0.filtered.genes.no_metadata.faa  bin_5.filtered.genes.no_metadata.faa
    bin_0.filtered.genes.no_metadata.fna  bin_5.filtered.genes.no_metadata.fna
    bin_1.filtered.genes.faa              bin_6.filtered.genes.faa
    bin_1.filtered.genes.fna              bin_6.filtered.genes.fna
    bin_1.filtered.genes.gbk              bin_6.filtered.genes.gbk
    bin_1.filtered.genes.no_metadata.faa  bin_6.filtered.genes.no_metadata.faa
    bin_1.filtered.genes.no_metadata.fna  bin_6.filtered.genes.no_metadata.fna
    bin_2.filtered.genes.faa              bin_7.filtered.genes.faa
    bin_2.filtered.genes.fna              bin_7.filtered.genes.fna
    bin_2.filtered.genes.gbk              bin_7.filtered.genes.gbk
    bin_2.filtered.genes.no_metadata.faa  bin_7.filtered.genes.no_metadata.faa
    bin_2.filtered.genes.no_metadata.fna  bin_7.filtered.genes.no_metadata.fna
    bin_3.filtered.genes.faa              bin_8.filtered.genes.faa
    bin_3.filtered.genes.fna              bin_8.filtered.genes.fna
    bin_3.filtered.genes.gbk              bin_8.filtered.genes.gbk
    bin_3.filtered.genes.no_metadata.faa  bin_8.filtered.genes.no_metadata.faa
    bin_3.filtered.genes.no_metadata.fna  bin_8.filtered.genes.no_metadata.fna
    bin_4.filtered.genes.faa              bin_9.filtered.genes.faa
    bin_4.filtered.genes.fna              bin_9.filtered.genes.fna
    bin_4.filtered.genes.gbk              bin_9.filtered.genes.gbk
    bin_4.filtered.genes.no_metadata.faa  bin_9.filtered.genes.no_metadata.faa
    bin_4.filtered.genes.no_metadata.fna  bin_9.filtered.genes.no_metadata.fna
    ```

## 2. Obtain annotations for genes of interest

### 2.1 Subset annotations

Based on the KEGG module for [sulfate-sulfur assimilation](https://www.kegg.jp/module/M00616), we need annotations that match `K02048`, `K23163`, `K02046`, `K02047`, and `K02045`. Out of all the columns, we are only interested in the gene ID and the KEGG annotations. We know that column 1 is the gene ID, so lets find the column index that contains the KEGG annotations.

!!! terminal "code"
    
    ```bash
    grep -P "K0204[5-8]|K23163" example_annotation_tables/*.aa \
      | awk -F '\t' '
          {
            for (i = 1; i <= NF; i++) {
              if (match($i, "K[0-9]{5}"))
                printf("Column %d contains KO of interest.\n", i)
            }
          }
        ' \
      | uniq

    # Column 26 contains KO of interest.
    ```

Deconstruct the code block above:

| Code | Action |
| :--  | :--    |
| `grep -P "K0204[5-8]\|K23163" *.aa` | Finds the rows that contain the pattern for our KO of interest in all annotation files within this directory. |
| `awk -F '\t' '` | Starts the `awk` programme to parse the output of `grep` and tells it that the outputs have tab-delimited columns. |
| `for (i = 1; i <= NF; i++) {` | Starts the C-style `for` loop within `awk` by assigning it an initial value of 1 (`i = 1`), continue the loop as long as the value of `i` is lower than the total number of rows (`i <= NF`), and increment by 1 as the loop progresses (`i++`). |
| `if (match($i, "K[0-9]{5}"))` | Starts an `if` loop and evaluates if a column contains string that matches the pattern `K[0-9]{5}` (which is the regex for a KO number). This is evaluated per row. |
| `printf("Column %d contains KO of interest.\n", i)` | If there is a match in the row, print out the resulting column index. |
| `uniq` | Returns only unique results. |

The code above says that all KO annotations are in column 26. This is true for all annotation files. Lets move on to subset to the KO numbers we want and cut out the columns we need. We will do this per bin.

!!! terminal "code"
    
    ```bash
    for bin_number in {0..9}; do
    
      grep -P "K0204[5-8]|K23163" example_annotation_tables/bin_${bin_number}.annotation.aa \
        | cut -f 1,26 > bin_${bin_number}.goi.aa
    
      if [ ! -s bin_${bin_number}.goi.aa ]; then
        rm bin_${bin_number}.goi.aa
      fi
    
    done
    ```

??? info "Deconstructing the code block"

    The above code block loops through each bin annotation file to subset relevant rows from the main annotation files (`grep -P "K0204[5-8]|K23163" example_annotation_tables/bin_${bin_number}.annotation.aa`) then selects the first and 26<sup>th</sup> columns (`cut -f 1,26`). 
    
    It also evaluates if the output is an empty file (in the case where the KO of interest is not found in the bin annotation) (`if [ ! -s bin_${bin_number}.goi.aa ]`) and removes it if it is. 
    
We have 4 files after running the above code.

!!! terminal "code"
    
    ```bash
    ls *.aa
    ```

!!! circle-check "Terminal output"

    ```
    bin_3.goi.aa  bin_5.goi.aa  bin_8.goi.aa  bin_9.goi.aa
    ```

### 2.2 Check for contiguous and sequential outputs

We need to check that the genes are contiguous (all on the same contig) and sequential (identify gaps between genes and fill them). Remember that prodigal headers (for which were propagated through our annotation workflows) always have `contigID_geneID` where `geneID` is always relative to the order it was found on the contig's nucleotide sequence. With that in mind, lets check for gene contiguity first.

!!! terminal "code"
    
    ```bash
    for aa_file in *goi.aa; do
      bin=$(basename ${aa_file} .goi.aa)
      contig=$(cut -f 1 ${aa_file} | sed -E 's/_[0-9]+$//g' | sort -u | wc -l)
    
      if [ $contig -gt 1 ]; then
        printf "%s has %d contig(s).\n" $bin $contig
      fi
      
    done
    ```

!!! circle-check "Terminal output"

    ```
    bin_5 has 2 contig(s).
    ```

??? note "Deconstructing the code block"
    
    The code above loops through all annotation subsets and finds the number of contigs in each file. It does this by using `cut -f 1` which selects the first columns, then `sed -E 's/_[0-9]+$//g'` removes the trailing gene numbers of the gene ID, followed by a dereplication of contigs via `sort -u` and finally counts the number of lines of the results using `wc -l`.

    if the bin annotation subset has more than 1 contig (`if [ $contig -gt 1 ]; then`), then it reports/prints the bin and number of contigs found (`printf "%s has %d contig(s).\n" $bin $contig`). 

The code block above shows that we have non-contiguous genes in bin 5. We will need to remove the extra gene(s) originating from another contig. Lets check which contig we need to remove.

!!! terminal "code"
    
    ```bash
    cat bin_5.goi.aa
    ```

!!! circle-check "Terminal output"

    ```    
    bin_5_NODE_52_length_158668_cov_0.373555_136    K02048: cysP; sulfate ABC transporter substrate-binding protein
    bin_5_NODE_95_length_91726_cov_0.379302_67      K02048: sbp; sulfate-binding protein
    bin_5_NODE_95_length_91726_cov_0.379302_68      K02046: cysT; sulfate transporter CysT
    bin_5_NODE_95_length_91726_cov_0.379302_69      K02047: cysW; sulfate transporter CysW
    bin_5_NODE_95_length_91726_cov_0.379302_70      K02045: cysA; sulfate.thiosulfate ABC transporter ATP-binding protein CysA
    ```

We find that contig `bin_5_NODE_52_length_158668_cov_0.373555` is not contiguous with the set of genes in contig `bin_5_NODE_95_length_91726_cov_0.379302`. We will make a note of that. Now lets check for continuity of genes.

!!! terminal "code"
    
    ```bash
    for aa_file in *goi.aa; do
      bin=$(basename ${aa_file} .goi.aa)
      gene_number=$(cut -f 1 ${aa_file} | sed -E 's/.*_([0-9]+$)/\1/g' | uniq)
      echo "$bin: $(echo ${gene_number[@]})"
    done
    ```
!!! circle-check "Terminal output"

    ```
    bin_3: 128 132 133 134
    bin_5: 136 67 68 69 70
    bin_8: 55 56 57 58
    bin_9: 147 148 149 150
    ```

??? info "Deconstructing the code block"

    We loop through each of the annotation subset to select the first column (`cut -f 1`), extract the trailing gene order/number (`sed -E 's/.*_([0-9]+$)/\1/g'`) and returns unique entries (`uniq`). The result is stored as an array named `gene_number`.

    It also reports the results of each bin and the array of gene numbers.

Results from the above code block shows that bin 3 has a gap between genes. We need to make sure we get the genes in between gene number 128 and 132.

### 2.3 Create subsets of annotation files with correct genes

Based on previous checks, we find that we need to:

- remove contig `bin_5_NODE_52_length_158668_cov_0.373555` from bin 5 when we create the final annotation file.
- fill in gaps in genes for bin 3.

!!! terminal "code"
    
    ```bash
    # Remove non-contiguous gene from bin_5.goi.aa
    sed -i '/bin_5_NODE_52_length_158668_cov_0.373555/d' bin_5.goi.aa
    
    # Obtain genes between 129 and 131 in contig bin_3_NODE_53_length_158395_cov_1.135272 and add to annotations
    cat example_annotation_tables/bin_3.annotation.aa \
      | grep -E "bin_3_NODE_53_length_158395_cov_1.135272_" \ # Search for required contig
      | grep -E "_129|_130|_131" \                            # Search for required gene numbers
      | cut -f 1,26 >> bin_3.goi.aa                           # Select columns 1 and 26
    ```

### 2.4 Clean up annotation tables
Now we have gene annotations that are contiguous and continuous. Lets clean up the tables by:

1. Sorting entries by gene order
1. Separating KO from other annotation information
1. Adding headers
1. Replace anything that doesn't have KO numbers with an annotation so we do not get empty fields

!!! terminal "code"
    
    ```bash
    for aa_file in *.aa; do
      bin=$(basename ${aa_file} .goi.aa)
      sort -u -k 1 ${aa_file} \
        | sed -e 's/: /\t/g' \
        | sed '1i\Query gene\tKO\tAnnotation' \
        | awk -F '\t' -v OFS='\t' '
            {
              if (!match($2, "K")) {
                $3=$2
              }
              
              print
            }
          ' > ${bin}_cys.txt
    done
    ```

??? info "Deconstructing the code block"
    
    We loop through each annotation subset, then:

    | Code | Action |
    | :--  | :--    |
    | `sort -u -k 1` | Dereplicate and sort entries based on the first column |
    | `sed '1i\Query gene\tKO\tAnnotation'` | Add header |
    | `awk -F '\t' -v OFS='\t'` | Initiate the `awk` programme by setting the input and output field separator as tabs |
    | `if (!match($2, "K")) {$3=$2}` | If entries in the second column that **do not** have a "K", use value in the second column as values for the third column |
    | `print` | Return/print out all lines |

## 3. Create comparison tables

After obtaining relevant annotations for our genes of interest, we need to align the sequences against each other in order to identify sequence homology and directionality across our genes of interest. This can be achieved via pairwise sequence alignment using BLAST (or BLAST-like algorithms). Here, we will use `tBLASTx` to obtain requisite comparison tables. Here, we will do this using the command line. However, you can also use online BLAST as [outlined here](./5_APPENDIX_ex15_gene_synteny_Generate_blast_files.md). This is a 2 step process:

1. Obtain nucleotide sequences for relevant contigs in each bin
1. Run (sequentially) a pairwise `tBLASTx` between each bin

??? note "Which BLAST algorithms to use?"
    
    For visualising gene synteny using `genoPlotR`, one can use outputs of `BLASTn` or `tBLASTx`. The determination of which algorithm to use depends on the phylogenetic relationship between MAGs. For closely related genomes, you can use `BLASTn` (or its variants) given that the genes will likely be conserved at the nucleotide sequence level. However, if your genomes are phylogenetically distant, you will need to use `tBLASTx`. This allows us to compare gene homology at the amino acid sequence level and retain nucleotide position information.

### 3.1 Subset contigs per bin

We begin by getting contig headers from our annotation files then using that to search and subset binned contigs. To subset sequences in FASTA format, we will be using `seqtk`

!!! terminal "code"
    
    ```bash
    module purge
    module load seqtk/1.3-gimkl-2018b
    
    for gene_file in *_cys.txt; do
      bin=$(basename ${gene_file} _cys.txt)
      cut -f 1 ${gene_file} | tail -n+2 | sed -e 's/_[0-9]*$//g' | uniq > ${bin}_cys.contigID
      seqtk subseq filtered_bins/${bin}.filtered.fna ${bin}_cys.contigID > ${bin}_cys_contig.fna
    done
    ```

### 3.3 Run `tBLASTx`

We do not need to perform pairwise comparisons for all bin combinations (you can, if you want to). We will run tBLASTx in the order which we want to observe the genes in: `bin_3`, `bin_5`, `bin_8`, then `bin_9`. To do this, we set an array for bin IDs then run them through sequentially.

!!! terminal "code"
    
    ```bash
    # Set array
    bin_array=(bin_3 bin_5 bin_8 bin_9)
    
    # Initiate index
    i=0
    
    # Set maximum array index
    max_i=$(( ${#bin_array[@]}-1 ))
    
    # Run tBLASTx, comparing bins sequentially in the order of the array
    while [ $i -lt $max_i ]; do
      qbin=${bin_array[$i]}
      sbin=${bin_array[$(($i+1))]}
      printf "Comparing %s against %s\n" $qbin $sbin
      tblastx -query ${qbin}_cys_contig.fna -subject ${sbin}_cys_contig.fna \
              -out blast_${qbin}_${sbin}.txt -outfmt "7 std ppos frames"
      ((i++)) # Increment i after each comparison
    done
    ```

---

The above workflow should generate the required files for the annotation subset (`<binID>_cys.txt`) and pairwise `BLAST` comparisons (`blast_<query binID>_<subject binID>.txt`). Along with these outputs, copy the `prodigal` predictions (either `*.faa` or `*.fna`; do not use `no_metadata` files) from relevant bins into your working directory and follow the steps outlined in [Presentation of data: gene synteny](../day4/ex16e_data_presentation_Gene_synteny.md) to generate synteny plots.

