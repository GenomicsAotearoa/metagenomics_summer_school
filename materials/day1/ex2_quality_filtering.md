# Quality filtering raw reads

### Objectives

* Visualising raw reads with `FastQC` on NeSI
* Read trimming and adapter removal with `trimmomatic`
* Diagnosing poor libraries
* Understand common issues and best practices

---

### Visualising raw reads

`FastQC` is an extremely popular tool for checking your sequencing libraries, as the visual interface makes it easy to identify the following issues:

1. Adapter/barcode sequences
1. Low quality regions of sequence
1. Quality drop-off towards the end of read-pair sequence

To activate `FastQC` on NeSI, you need to first load the module using the command

```bash
module load FastQC/0.11.7
```

These exercises will take place in the `2.fastqc/` folder. First, navigate to this folder. Copy the command below into your terminal (logged in to NeSI), replacing `<YOUR FOLDER>`, and then running the command.

```bash
cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/2.fastqc/
```

`FastQC` can then be run either interactively (i.e. with a GUI) or from the command line. The difference is simply whether or not you provide any input files to `FastQC` when it loads.

```bash
# Load interactively

fastqc mock_R1.good.fastq.gz mock_R2.good.fastq.gz
```

Now download the resulting files using `scp` and view them on your computer. Note that `FastQC` does not load the forward and reverse pairs of a library in the same window, so you need to be mindful of how your samples relate to each other. 

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig1_overview.PNG)

At a first glance, we can see the follow statistics:

1. The data is stored in Sanger / Illumina 1.9 encoding. *This will be important to remember when read trimming.*
1. There are 100,000 reads in the file
1. The maximum sequence length is 251 base pairs. *This is good to check, since the data were generated using Illumina
2x250 bp sequencing.*

Have a quick look through the left hand columns. As you can see, the data has passed most of the basic parameters.

**Per-base sequence quality**

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig2_quality.PNG)

**Per base sequence content**

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig3_content.PNG)

**Adapter Content**

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig4_adapter.PNG)

**Per sequence GC content**

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig5_gc.PNG)

The only aspect of the data that `FastQC` is flagging as potentially problematic is the GC% content of the data set. This is a common observation, as we are dealing with a mixed community and organisms and it is therefore unlikely that there will be a perfect normal distribution around an average value. For example, a community comprised of low- and high-GC organisms would manifest a bimodal distribution of peaks which would be a problematic outcome in terms of the expectations of `FastQC`, but completely consistent with the biology of the system.

Lets take a look at a library with significant errors. Load the sequence file `mock_R1.adapter_decay.fastq` and compare the results with the `mock_R1.good.fastq` file.

Which of the previous fields we examined are now flagged as problematic? How does this compare with your expectation? Are there any which should be flagged which are not?

**Per-base sequence quality**

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig6_taildecay.PNG)

---

### Read trimming and adapter removal with *trimmomatic*

There are a multitude of programs which can be used to quality trim sequence data and remove adapter sequence. For this exercise we are going to use `trimmomatic`, but this should in no way be interpreted as an endorsement of `trimmomatic` over equivalent tools like `BBMap`, `sickle`, `cutadapt` or any other.

For a first run with `trimmomatic`, type the following commands into your console:

```bash
module load Trimmomatic/0.39-Java-1.8.0_144

trimmomatic PE -threads 10 -phred33 \
               mock_R1.adapter_decay.fastq.gz mock_R2.adapter_decay.fastq.gz \
               mock_R1.qc.fastq.gz mock_s1.qc.fastq.gz mock_R2.qc.fastq.gz mock_s2.qc.fastq.gz \
               HEADCROP:10 SLIDINGWINDOW:4:30 MINLEN:80
```

There is a lot going on in this command, so here is a breakdown of the parameters in the command above

|Parameter|Type|Description|
|:---|:---:|:---|
|PE|*positional*|Specifies whether we are analysing single- or paired-end reads|
|-threads 10|*keyword*|Specifies the number of threads to use when processing|
|-phred33|*keyword*|Specifies the fastq encoding used|
|mock_R1.adapter_decay.fastq.gz / mock_R2.adapter_decay.fastq.gz|*positional*|The paired forward and reverse reads to trim|
|mock_R1.qc.fastq.gz|*positional*|The file to write forward reads which passed quality trimming, if their reverse partner also passed|
|mock_s1.qc.fastq.gz|*positional*|The file to write forward reads which passed quality trimming, if their reverse partner failed (orphan reads)|
|mock_R2.qc.fastq.gz / mock_s2.qc.fastq.gz|*positional*|The reverse-sequence equivalent of above|
|HEADCROP:10|*positional*|Adapter trimming command. Remove the first 10 positions in the sequence|
|SLIDINGWINDOW:4:30|*positional*|Quality filtering command. Analyses each sequence in a 4 base pair sliding window, truncating if the average quality drops below Q30|
|MINLEN:80|*positional*|Length filtering command, discard sequences that are shorter than 80 base pairs after trimming|

This signifcantly improves the output.

##### Quality of reads

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig7_qualtrim.PNG)

##### Nucleotide distribution

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex3_fig8_headtrim.PNG)

---

### Considerations when working with *trimmomatic*

**Order of operations**

The basic format for a `trimmomatic` command is

```bash
trimmomatic PE <keyword flags> <sequence input> <sequence output> <trimming parameters>
```

The trimming parameters are processed in the order you specify them. This is a deliberate behaviour, but can have some unexpected consequences for new users.

For example, consider these two scenarios:

```bash
trimmomatic PE <keyword flags> <sequence input> <sequence output> SLIDINGWINDOW:4:30 MINLEN:80

trimmomatic PE <keyword flags> <sequence input> <sequence output> MINLEN:80 SLIDINGWINDOW:4:30
```

In the first run we would not expect any sequence shorter than 80 base pairs to exist in the output files, but we might encounter them in the second command. This is because in the second command we remove sequences shorter than 80 base pairs, **_then_** perform quality trimming. If a sequence is trimmed to a length shorter than 80 base pairs **_after_** trimming the `MINLEN` filtering does not execute a second time. In the first instance, we are excluding performing trimming **_before_** size selection, so any reads that start longer than 80 base pairs, but are trimmed to under 80 base pairs during quality trimming will be caught in the `MINLEN` run.

A more subtle consequence of this behaviour is the interplay between adapter removal and quality filtering. In the command above, we try to identify adapters **_before_** quality trimming. Why do you think this is?

---

### *Optional:* Working with the ILLUMINACLIP command

`Trimmomatic` also provides the `ILLUMINACLIP` command, which can be used to pass a *fastA* file of adapter and barcode sequences to be found and removed from your data. The format for the `ILLUMINACLIP` parameter can be quite confusing to work with and so we generally favour a `HEADCROP` where possitble. If you wish to work with the `ILLUMINACLIP` command, then you can see its use here in the `trimmomatic` [manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

```bash
ILLUMINACLIP:iua.fna:1:25:7
             |       | |  |
             |       | |  |
             |       | |  Simple clip threshold
             |       | Palindrome clip threshold
             |       Seed mismatches
             File of expected sequences
```

There is always some subjectivity in how sensitive you want your adapter (and barcode) searching to be. If the settings are too strict you might end up discarding real sequence data that only partially overlaps with the Illumina adapters. If your settings are not strict enough then you might leave partial adapters in the sequence. Where possible, we favour the use of simple positional trimming.

---

### Diagnosing poor libraries

Whether a library is 'poor' quality or not can be a bit subjective. These are some aspects of the library that you should be looking for when evaluating `FastQC`:

1. Does the sequencing length match what you ordered from the facility?
1. If the sequences are shorter than expected, is adapter read-through a concern?
1. What does the sequence quality look like for the whole length of the run? Are there any expected/unexpected regions of quality degradation?
1. Are adapters and/or barcodes removed? (look at the *Per base sequence content* to diagnose this)
1. Is there unexpected sequence duplication? (this can occur when low-input library preparations are used)
1. Are over-represented *k*-mers present? (this can be a sign of adapter and barcode contamination)

---

### Understand common issues and best practices

1. Do I need to remove (rare) adapters?
1. I don’t know if adapters have been removed or not
1. How do I identify and remove adapter read-through
1. Identifying incomplete barcode/adapter removal
1. Over aggressive trimming
1. GC skew is outside of expected range

---
