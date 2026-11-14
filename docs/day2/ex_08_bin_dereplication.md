# Binning II: Bin dereplication

!!! info "Objectives"

    - [Bin dereplication using `DAS_Tool` - Creating input tables](#bin-dereplication-using-das_tool-creating-input-tables)
    - [Bin dereplication using `DAS_Tool`](#bin-dereplication-using-das_tool)
    - [Dereplicating bins across multiple assemblies](#dereplicating-bins-across-multiple-assemblies)

---

## Bin dereplication using `DAS_Tool` - Creating input tables

As we discussed in the previous exercise, we have now generated three sets of bins from the same single assembly. With this mock data set we can see that `MetaBAT` recovered 12 bins, while `MaxBin` and `MetaDecoder` recovered 10. Note that we are aiming to recover prokaryote genomes using these binning tools, and 10 bacterial and archaeal genomes were used in the creation of this mock community. In our case here, these tools have likely recovered 10 bins each of the same genomes. The additional two bins identified by `MetaBAT` may be the result of noise introduced into the binning process by the viral contigs included in the data. Furthermore, it is not clear which tool has done a better job of recruiting contigs to each bin - we very rarely expect to see the complete genome recovered from these kinds of data, so while it is probably the case that an equivalent bin is present in the `MetaBAT`, `MaxBin` and `MetaDecoder` outputs, they will likely be of differing quality.

`DAS_Tool` is a program designed to analyse the bins in each of our binning sets and determine where these equivalent triplets exist and return the 'best' one. `DAS_Tool` does not use the actual bins, but a set of text files that link contigs to their corresponding bins in each of the bin sets. We can produce these files using `bash`.

For this exercise, we will continue working in the `4.binning/` directory.

### Create contig/bin tables - `MetaBAT`

For each of our binning tools, we need to extract the contigs assigned to each bin and create a single file that reports these as

!!! success ""

    ```
    Contig[tab]Bin
    ```

This can be done with a bit of `bash` scripting. There's quite a bit going on here, so we'll provide the full command, and then a step-by-step explanation of what's happening.

!!! warning "Remember to update `<YOUR FOLDER>` to your own folder"

!!! terminal-2 "Navigate to working directory"

    ```bash
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/4.binning/
    ```

!!! terminal-2 "Create `metabat_refs.tsv`"

    ```bash linenums='1'
    > metabat_refs.tsv
    for bin_path in metabat/*.fa; do
        bin_name=$(basename ${bin_path} .fa)
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> metabat_refs.tsv
    done
    ```

You can check the contents of this file using `less` or `head`:

!!! terminal-2 "Inspect `metabat_association.txt`"

    ```bash
    head -n 5 metabat_refs.tsv
    ```

!!! circle-check "Terminal output"

    ```
    NODE_3_length_835083_cov_0.372452   metabat.10
    NODE_5_length_586102_cov_0.372632   metabat.10
    NODE_20_length_292112_cov_0.373169  metabat.10
    NODE_24_length_265402_cov_0.371105  metabat.10
    NODE_33_length_247431_cov_0.373159  metabat.10
    ```

We will now walk through the content of the "Create `metabat_refs.tsv`" code block above, breaking apart each individual step.


!!! terminal-2 "Line 1: Create/ensure empty output file"

    ```bash
    > metabat_refs.tsv
    ```

We don't want to append to the contents of a potentially existing file: As literately nothing is on the left side of the redirection operator `>`, an empty file is created or the content of an existing one deleted.

!!! terminal-2 "Line 2: Initiate a `for` loop"

    ```bash
    for bin_path in metabat/*.fa; do
    ```

This is the initial loop, returning each file in the `metabat/` folder that ends with the file extension *.fa*.

!!! terminal-2 "Line 3: Remove path prefix"

    ```bash
        bin_name=$(basename ${bin_path} .fa)
    ```

As we have previously seen, the `basename` command removes the path information from the input variable `bin_path` (i.e. - `metabat/metabat.1.fa` becomes `metabat.1.fa`) and assigns it to a new variable, `bin_name`. As an optional additional parameter, we can pass extra pieces of text to be removed from the variable, in this case the *.fa* extension.

!!! terminal-2 "Line 4-5: Subset information to create and update `metabat_refs.tsv`"

    ```bash
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> metabat_refs.tsv
        done
    ```

We will break this line down step-by-step:

| <div style="width:230px">Command</div> | Description |
| :--- | :--- |
| `grep ">" ${bin_path}` | The `grep` command searches the input file `bin_path` for lines containing the `>` character, which is the *FASTA* demarcation for a sequence name. |
| `sed 's/>//g'` | The first `sed` command replaces the `>` character with the empty character `''`, as we do not need this character in the final file. |
| `sed "s/$/\t${bin_name}/g"` | We now use `sed` again, this time to replace the `$` character. In many command line tools and software environments (including `R` and `python`) the characters `^` and `$` are used as shortcuts for the beginning and ending of a line, respectively. By using the character in this way, we are telling `sed` to replace the end-of-line with the text `\t${bin_name}`. `sed` will parse this text to mean the tab character followed by the content of the `bin_name` variable. The nature of `sed` is that it will automatically insert a new end-of-line. |
| `>> metabat_refs.tsv` | This is similar to the stdout redirection we have previously used, but the double use of the `>` character means that we are appending our text to the end of the file `metabat_refs.tsv`. Because we are looping through several files in this exercise, if we were to use the single `>` character then on each new *fastA* file read, the content of `metabat_refs.tsv` would be replaced.<br><br>It is important to note that because we are **_appending_** to the file, not **_replacing_** the contents, if you make a mistake in the command and need to re-run it, you will need to explicitly delete the `metabat_refs.tsv` file using `rm`, otherwise your new (correct) output will be pasted to the end of your old (incorrect) output.|

### Create contig/bin tables - `MaxBin`, `MetaDecoder`

The process for creating the `MaxBin` and `MetaDecoder` tables is basically the same: Because we have already normalised all bin filename extensions to *.fa*, we only have to change the names for the input directory and the output file.

!!! terminal-2 "Create `maxbin_refs.tsv`"

    ```bash linenums="1"
    > maxbin_refs.tsv
    for bin_path in maxbin/*.fa; do
        bin_name=$(basename ${bin_path} .fa)
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> maxbin_refs.tsv
    done
    ```

!!! terminal-2 "Create `mdecoder_refs.tsv`"

    ```bash linenums="1"
    > mdecoder_refs.tsv
    for bin_path in metadecoder/*.fa; do
        bin_name=$(basename ${bin_path} .fa)
        grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> mdecoder_refs.tsv
    done
    ```


!!! warning "Warning for unbinned contigs"

    Both `MetaBAT` and `MaxBin` have the option to output unbinned contigs after binning completes. We have not used that parameter here, but if you do choose to enable it you will end up with another *FASTA* file in your output directory which you will need to avoid in the loops for creating `DAS_Tool` tables.

---

## Bin dereplication using `DAS_Tool`

We are now ready to run `DAS_Tool`. This can be done from the command line, as it does not take a particularly long time to run for this data set.

!!! terminal-2 "Run `DAS_Tool`" 

    ```bash linenums="1"
    # Remove modules to ensure a clean environment
    module purge

    # Load DAS Tool
    module load DAS_Tool/1.1.5-gimkl-2022a-R-4.2.1

    # Create DAS_Tool output directory
    mkdir -p dastool_out/

    # Run DAS_Tool
    DAS_Tool --threads 2  --write_bins \
             --search_engine diamond \
             --bins metabat_refs.tsv,maxbin_refs.tsv,mdecoder_refs.tsv \
             --labels MetaBAT,MaxBin,MetaDecoder \
             --contigs spades_assembly/spades_assembly.m1000.fa \
             --outputbasename dastool_out/
    ```

!!! circle-check "Terminal output"

    ```
    DAS Tool 1.1.5 
    Analyzing assembly 
    Predicting genes 
    Annotating single copy genes using diamond 
    Dereplicating, aggregating, and scoring bins 
    Writing bins
    ```

As usual, we will break down the parameters:

|<div style="width:200px">Parameter</div>|Description|
|:---|:---|
|`--threads`|Number of threads to use|
|`--write_bins`|An argument telling `DAS_Tool` whether or not to write out a new set of bins<br>This is recommended, because `DAS_Tool` can create slices of old bins based on marker composition (see [the paper](https://www.nature.com/articles/s41564-018-0171-1) for details)|
|`--search_engine diamond`|Specify whether to use `usearch`, `diamond`, or `BLAST` as the alignment tool for comparing gene sequences|
|`--bins`|A comma-separated list of the contig/bin files we wish to process|
|`--labels`|A comma-separated list of the binning tools used|
|`--contigs`|Path to the assembly used in binning|
|`--outputbasename`|Output directory for all files|

When `DAS_Tool` has completed, we will have a final set of bins located in the folder path `dastool_out/_DASTool_bins`. Have a look at the output and see which bins made it to the final selection. Did a single binning tool pick the best bins, or are the results a split between `MetaBAT`, `MaxBin` and `MetaDecoder`?

---

## Dereplicating bins across multiple assemblies

In this workshop, we have generated a set of putative MAGs by binning scaffolds taken from a *single co-assembly* of all four samples. Alternatively, we may have chosen to generate multiple assemblies (for example, mini-co-assemblies for each sample group, or individual assemblies for each sample). In this case, it would be necessary to work through the binning process for each assembly, and then conduct an additional dereplication step *across* the multiple assemblies to generate a single set of dereplicated bins for all assemblies. 

