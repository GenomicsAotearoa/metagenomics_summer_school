# Introduction to the command line

---

### Objectives

* Navigating your file system
* Copying, Moving, Renaming and Removing files
* Examining file contents
* Redirection, manipulation and extraction

---

### Navigating your file system

Find the current location by running the command `pwd` which stands for *print working directory*. At any given time, **current working directory** is your current default directory.

```bash
pwd
# /home/UserName/
```

We can see what files and subdirectories are in this directory by running `ls`, which stands for "listing".

```bash
ls
```

Navigating to the `MGSS_Intro/` directory can be done with the`cd` command which stands for *change directory*.

```bash
cd MGSS_Intro/
```

Run the `ls` command to list the content of the current directory. Check whether there are two *.fastq* files.

The `mkdir` command (*make directory*) is used to make a directory. Enter `mkdir` followed by a space, then the directory name you want to create

```bash
mkdir backup/
```

---

### Copying, Moving, Renaming and Removing files

Make a second copy of `Test_1.fastq` and rename it as `Test_1_backup.fastq`.  Then move that file to `backup/` directory.

```bash
cp Test_1.fastq Test_1_backup.fastq

mv Test_1_backup.fastq backup
```

Navigate to the `backup/` directory and use the `mv` command to rename and move `Test_1_backup.fastq` as `Test_1_copy.fastq` to the directory immediately above.

```bash
cd backup/

mv Test_1_backup.fastq ../Test_1_copy.fastq
```

Return to the directory immediately above, check whether the `Test_1_copy.fastq` was moved and renamed as instructed and remove it by using the `rm` command.

```bash
cd ..

rm Test_1_copy.fastq
```

See whether you can remove the `backup/` directory by using the `rm` command as well. 

```bash
rm backup/
# rm : can not remove 'backup/': Is a directory
```

By default, `rm` will not delete directories. This can be done by using the `-r` (recursive) option.

```bash
rm -r backup
```

---

### Examining file contents

There are a number of ways to examine the content of a file. `cat` and `less` are two commonly used programs for a quick look. Check the content of `Test_1.fastq` by using these commands. Take a note of the differences. 

```
cat Test_1.fastq
# less Test_1.fastq
```

A few useful shortcuts for navigating in `less`

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex1_less_shortcuts.jpg) 

There are ways to take a look at parts of a file. For example the `head` and `tail` commands will scan the beginning and end of a file, respectively. 

```bash
head Test_1.fastq

tail Test_1.fastq
```

Adding `-n` option to either of these commands will print the first or last *n* lines of a file.

```bash
head -n 1 Test_1.fastq
# @SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
```

---

### Redirection, manipulation and extraction

Although using `cat` and `less` commands will allow us to view the content of the whole file, most of the time we are in search of particular characters (strings) of interest, rather than the full content of the file. One of the most commonly used command-line utilities to search for strings is `grep`. Let's use this command to search for the string NNNNNNNNNN in `Test_2.fastq` file.

```bash
grep NNNNNNNNNN Test_2.fastq
```

Retrieve and discuss the output you get when `grep` was executed with the `-B1` and `-A1` flags.

```bash
grep -B1 -A2 NNNNNNNNNN Test_2.fastq
```

In both occasions, outputs were printed to the terminal where they can not be reproduced without the execution of the same command. In order for "string" of interest to be used for other operations, this has to be "redirected" (captured and written into a file). The command for redirecting output to a file is `>`. Redirecting the string of bad reads that was searched using the `grep` command to a file named `bad_reads.txt` can be done with

```bash
grep -B1 -A2 NNNNNNNNNN Test_2.fastq > bad_reads.txt
```

Use the `wc` command to count the number of words, lines and characters in the `bad_reads.txt` file.

```bash
wc bad_reads.txt
```

Add `-l` flag to `wc` command and compare the number with the above output

```bash
wc -l bad_reads.txt
```

In an instance where the same operation has to be applied to multiple input files and the outputs are to be redirected to the same output file, it is important to make sure that the new output is not over-writing the previous output. This can be avoided with the use of `>>` (append redirect) command which will append the new output to the end of the file, rather than overwriting it. 

```bash
grep -B1 -A2 NNNNNNNNNN Test_1.fastq >> bad_reads.txt
```

Executing the same operation on multiple files with the same (or different) file extension can be done with the use of *wildcards* which are symbols or special characters that represent other characters. For an example. Using `*` wildcard, we can run the previous `grep` command on both files at the same time

```bash
grep -B1 -A2 NNNNNNNNNN *.fastq >> bad_reads.txt
```

The only objective of the above redirection example is to search for a string in a file(s), write the output to a file and then count the number of lines in the output. Generating output files for short routine tasks like that will end up generating an excessive number of files with little value. The `|` (pipe) command is a commonly used method to apply an operation for an ouput without creating intermediate files. It takes the output generated by one command and uses that output as input to another command. 

```bash
grep -B1 -A2 NNNNNNNNNN Test_2.fastq | wc -l
```

In an instance where we aren't after a particular string, but want sections from each line of file(s) to be extracted and written into an output file for further operations, `cut` command is a great utility. It can be used to cut parts of a line by **byte position, character and field**. Basically the cut command slices a line and extracts the text.
