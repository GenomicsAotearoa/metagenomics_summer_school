# APPENDIX (ex11) : viral taxonomy prediction via *vContact2*


**1. Predict genes via `prodigal`**

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy
```

Example slurm script:

!!! terminal "code"

    ```bash
    #!/bin/bash -e

    #SBATCH --account       nesi02659
    #SBATCH --job-name      prodigal
    #SBATCH --time          00:05:00
    #SBATCH --mem           1GB
    #SBATCH --cpus-per-task 2
    #SBATCH --error         prodigal.err
    #SBATCH --output        prodigal.out


    # Load dependencies
    module purge
    module load prodigal/2.6.3-GCC-11.3.0

    # Set up working directories
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy

    mkdir -p viral_taxonomy

    # Run main analyses 
    srun prodigal -p meta -q \
    -i checkv_combined.fna \
    -a viral_taxonomy/checkv_combined.faa 
    ```

**2. Generate required mapping file for `vContact2`**

Use `vContact2`'s `vcontact2_gene2genome` script to generate the required mapping file from the output of `prodigal`.

*NOTE: update `/path/to/conda/envs/vContact2/bin` in the below script to the appropraite path.*

!!! terminal "code"

    ```bash
    # activate vcontact2 conda environment
    module purge
    module load Miniconda3
    source activate vContact2

    # Load dependencies
    export PATH="/path/to/conda/envs/vContact2/bin:$PATH"
    module load DIAMOND/0.9.32-GCC-9.2.0
    module load MCL/14.137-gimkl-2020a

    # run vcontact2_gene2genome
    vcontact2_gene2genome -p viral_taxonomy/checkv_combined.faa -o viral_taxonomy/viral_genomes_g2g.csv -s 'Prodigal-FAA'

    # deactivate conda environment
    conda deactivate
    ```

**3. Run `vContact2`**

Example slurm script:

!!! terminal "code"

    ```bash
    #!/bin/bash -e
    
    #SBATCH --account       uow03498
    #SBATCH --job-name      vcontact2_test
    #SBATCH --time          02:00:00
    #SBATCH --mem           20GB
    #SBATCH --cpus-per-task 10
    #SBATCH --error         vcontact2_%j.err
    #SBATCH --output        vcontact2_%j.out
    
    
    # activate vcontact2 conda environment
    module purge
    module unload XALT
    module load Singularity
    
    container=/opt/nesi/containers
    
    # Bind external filesystem paths to container image
    export SINGULARITY_BIND="$PWD,$container,/nesi/nobackup/nesi02659/MGSS_U"
    
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/viral_taxonomy/
    
    
    # Run vcontact2
    singularity run $container/vcontact2.simg vcontact2 \
    -t $SLURM_CPUS_PER_TASK \
    --raw-proteins checkv_combined.faa \
    --rel-mode Diamond \
    --proteins-fp viral_genomes_g2g.csv \
    --db 'ProkaryoticViralRefSeq201-Merged' \
    --c1-bin /opt/conda/bin/cluster_one-1.0.jar \
    --output-dir vConTACT2_Results
    ```

**4. Predict taxonomy of viral contigs based on ouput of `vContact2`**

`vContact2` doesn't actually *assign* taxonomy to your input viral contigs. It instead provides an output outlining which reference viral genomes your viral contigs clustered with (if they clustered with any at all). Based on how closely they clustered with any reference genome(s), you can then use this to *predict* the likely taxonomy of the contig. 

From the `vContact2` online docs:

!!! quote ""
    One important note is that the taxonomic information is not included for user sequences. This means that each user will need to find their genome(s) of interest and check to see if reference genomes are located in the same VC. If the user genome is within the same VC subcluster as a reference genome, then there's a very high probability that the user genome is part of the same genus. If the user genome is in the same VC but not the same subcluster as a reference, then it's highly likely the two genomes are related at roughly genus-subfamily level. If there are no reference genomes in the same VC or VC subcluster, then it's likely that they are not related at the genus level at all.

The summary output of `vContact2` is the file `vConTACT2_Results/genome_by_genome_overview.csv`. As the comment above notes, one approach would be to search this file for particular contigs of interest, and see if any reference genomes fall into the same viral cluster (VC), using this reference to predict the taxonomy of the contig of interest.

The following `python` script is effectively an automated version of this for all input contigs (*Note: this script has not been widely tested, and so should be used with some degree of caution*). This script groups together contigs (and reference genomes) that fall into each VC, and then for each, outputs a list of all taxonomies (at the ranks of 'Order', 'Family', and 'Genus', separately) that were found in that cluster. The predictions (i.e. the list of all taxonomies found in the same VC) for each rank and each contig is output to the table `tax_predict_table.txt`. 

!!! note "Note"
    The taxonomies are deliberately enclosed in square brackets (`[ ]`) to highlight the fact that these are **predictions**, rather than definitive taxonomy **assignments**.*

For future reference, a copy of this script is available for download [here](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/scripts/vcontact2_tax_predict.py)

!!! terminal "code"

    ```bash
    module purge
    module load Python/3.8.2-gimkl-2020a
    
    cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/8.coverage_and_taxonomy/
    
    ./vcontact2_tax_predict.py \
    -i viral_taxonomy/vConTACT2_Results/genome_by_genome_overview.csv \
    -o viral_taxonomy/
    ```

---
