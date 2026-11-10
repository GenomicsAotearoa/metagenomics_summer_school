# BONUS: Varying assembly parameters

!!! info "Objectives"

    * [Changing parameters for `SPAdes` assembly](#examining-the-effect-of-changing-assembly-parameters-for-spades)
    * [Changing parameters for `IDBA-UD` assembly](#examining-the-effect-of-changing-assembly-parameters-for-idba-ud)

All work for this exercise will occur in the `2.assembly/` directory.

---

## Examining the effect of changing assembly parameters for `SPAdes`

For this exercise, there is no real structure. Make a few copies of your initial slurm scripts and tweak a few of the assembly parameters. Make sure to also change the assembly output file (e.g. from `spades_assembly` to `spades_assembly_var1`). We will compare the effects of these changes in the next lesson.

Make a few copies of your `SPAdes` slurm script like so;

!!! terminal "code"

    ```bash
    cp spades_assembly.sl spades_assembly_var1.sl
    ```

Change a few of the parameters for run time. Some potential options include

1. Change the *k*-mer sizes to either a different specification, or change to the `auto` option
1. Disable error correction
1. Assemble without the `--meta` flag
1. Employ a coverage cut-off for assembling

Submit two or three jobs per variation.

!!! warning "Outputs per variant job"

    For all variations of your assemblies, please remember to modify `-o` argument from `-o spades_assembly/` in your slurm scripts to another name (e.g. `-o spades_assembly_var_kmer/`)

---

## Examining the effect of changing assembly parameters for `IDBA-UD`

Make variants of your `IDBA-UD` assembly script and change some parameters. Some potential options include

1. Change the minimum/maximum *k*-mer sizes, or the *k*-mer step size
1. Change the alignment similarity parameter
1. Adjust the prefix length for the *k*-mer sub-table 

Submit two or three jobs per variation.

!!! warning "Outputs per variant job"

    For all variations of your assemblies, please remember to modify `--out` argument from `--out idbaud_assembly/` in your slurm scripts to another name (e.g. `--out idbaud_assembly_var_kmer/`)


---

Now try evaluating how well the new assembly/assemblies did by adapting code in the previous task on assembly evaluation (i.e. using seqmagick and BBMap’s stats.sh command). Note that the name of the assembly input directory will be different.
