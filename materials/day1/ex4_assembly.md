# Assembly (part 2)

### Objectives

* Examine the effect of changing parameters for assembly

All work for this exercise will occur in the `3.assembly/` directory.

---

### Examine the effect of changing assembly parameters

For this exercise, there is no real structure. Make a few copies of your initial slurm scripts and tweak a few of the asembly parameters. You will have a chance tomorrow to compare the effects of these changes tomorrow.

#### *SPAdes* parameters

Make a few copies of your `SPAdes` slurm script like so;

```bash
cp assembly_spades.sl assembly_spades_var1.sl
```

Change a few of the parameters for run time. Some potential options include

1. Change the *k*-mer sizes to either a different specification, or change to the `auto` option
1. Disable error correction
1. Assemble without the `--meta` flag
1. Employ a coverage cutoff for assembling

#### *IDBA-UD* parameters

Make variants of your `IDBA-UD` assembly script and change some parameters. Some potential options include

1. Change the minimum/maximum *k*-mer sizes, or the *k*-mer step size
1. Change the alignment similarity parameter
1. Adjust the prefix length for the *k*-mer sub-table 

Submit two or three jobs per variation.

---