# Dereplicating a data set with *dRep*

### Objectives

* Understand the common issues with using `dRep` and `CheckM`
* Use `CheckM` and `dRep` together to dereplicate a set of genomes

---

### Using *dRep* and *CheckM*

Before we begin to use `dRep`, it is important to understand the workflow that it applies to a data set. The basic idea of `dRep` is that genomes or MAGs are processed as follows:

1. Genomes are scored for completeness and contamination estimates using `CheckM`
1. Genomes are assigned to primary clusters using a quick and rough average nucleotide identidy (ANI) calculation
1. Clusters of genomes sharing greater than 90% ANI are grouped together and ANI is calculated using a more sensitive method
1. Where groups of genomes sharing >99% ANI are found, the best (determined by completeness and contamination statistics) is picked as a representative of the cluster

When run on its own, `dRep` will automatically try to run `CheckM` in the background. There are two problems with this approach, namely:

1. `dRep` is written in version 3.6 of the `python` language, and the version of `CheckM` avaiable on NeSI is written in version 2.7. These are not compatible with each other
1. There are two parameters in `CheckM` which speed up the workflow through multithreading, but `dRep` only has access to one of them

For these reasons, when working on NeSI we run `CheckM` on our data set first, and then pass the results directly into `dRep`, avoiding the need for `dRep` to try to call `CheckM`.

---

### Use *CheckM* and *dRep* together to dereplicate a set of MAGs

For this exercise, we will be working with a different set of MAGs to the mock community, as there is not enough strain-level variation in the mock metagenome for `dRep` to actually remove any MAGs.

We will write a single slurm script to run all necessary commands, then analyse the content.

```bash
#!/bin/bash -e
#SBATCH -A nesi02659
#SBATCH -J checkm_drep
#SBATCH --res SummerSchool
#SBATCH --time 2:00:00
#SBATCH --mem 80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -e checkm_drep.err
#SBATCH -o checkm_drep.out

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/10.drep_example/

# Step 1
module load CheckM/1.0.13-gimkl-2018b-Python-2.7.16
checkm lineage_wf -t 10 --pplacer_threads 10 -x fa --tab_table -f checkm.txt \
                  input_bins/ checkm_out/

# Step 2
echo "genome,completeness,contamination" > dRep.genomeInfo
cut -f1,12,13 checkm.txt | sed 's/\t/.fa\t/' | sed 's/\t/,/g' | \
    tail -n+2 >> dRep.genomeInfo

# Step 3
module purge
module load drep/2.3.2-gimkl-2017a-Python-3.6.3 MUMmer/3.23-gimkl-2017a

dRep dereplicate --genomeInfo dRep.genomeInfo -g input_bins/*.fa -p 10 drep_output/
```

Walking through this script, step by step, we are performing the following tasks:

##### Step 1

This should look familiar to you. Here we are simply loading `CheckM`, then running it over a set of MAGs to get the completeness and contamination estimates for our data.

##### Step 2

When running `dRep`, we have the option to either let `dRep` execute `CheckM` in the background, or we can pass a comma-separated file of the MAG name and its statistics. Unfortunately, `dRep` does not take the `CheckM` output itself, so we must use some shell commands to reformat the data. To achieve this, we use the following steps:

```bash
echo "genome,completeness,contamination" > dRep.genomeInfo
```

This line creates the header row for the `dRep` file, which we are calling `dRep.genomeInfo`.

```bash
cut -f1,12,13 checkm.txt | sed 's/\t/.fa\t/' | sed 's/\t/,/g' | tail -n+2 >> dRep.genomeInfo
```

This line cuts the columns 1, 12, and 13 from the `CheckM` output table, which correspond to the MAG name and completeness/contamination estimates. We then redirect these columns using the `|` character and use them as input in a `sed` command. Because `CheckM` reports our MAG names without their *fastA* file extension, but `dRep` requires the extension to be present in the MAG name, we add the trailing `.fa` with our `sed` command. This also gives us an opportunity to replace the tab-delimiting character from the `CheckM` output with the comma character that `dRep` uses for marking columns in the table. We then pass the output into a second `sed` command to replace the tab between columns 12 and 13 with a comma.

We then use another redirect (`|`) to pass the resulting text stream to the `tail` command. The way we are calling `tail` here will return every row in the text stream except for the first, which means that we are getting all the MAG rows but not their column names. We remove these names because ``dRep` uses a different naming convention for specifying columns.

We append the MAG statistics to the end of the `dRep.genomeInfo` file, whih contains the correct column names for `dRep`.

##### Step 3

Here we simply load the modules required for `dRep`, then execute the command. Because of the compatibility issues between the `python` version required by `CheckM` and `dRep`, we use the `module purge` command to unload all current modules before loading `dRep`. This removes the `CheckM` library, and its `python2.7` dependency, allowing the `dRep` and `python3.6` to load correctly.

The parameters for `dRep` are as follows:

|Parameter|Function|
|:---|:---|
|**dereplicate**|Activate the *dereplicate* workflow from `dRep`|
|**--genomeInfo ...**|Skip quality checking via `CheckM`, instead use the values in the table provided|
|**-g ...**|List of MAGs to dereplicate, passed by wildcard|
|**-p ...**|Number of processors to use|
|**drep_output/**|Output folder for all outputs|

When `dRep` finishes running, there are a few useful outputs to examine:

```bash
drep_output/dereplicated_genomes/   # The representative set of MAGs
drep_output/figures/                # Dendrograms to visualise the clustering of genomes
drep_output/data_tables/            # The primary and secondary clustering of the MAGs, and scoring information
```

---