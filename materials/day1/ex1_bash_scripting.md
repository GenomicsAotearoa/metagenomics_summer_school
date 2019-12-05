# Introduction to the command line

---

### Objectives

* Navigating your file system
* Copying, Moving, Renaming and Removing files
* Examining file contents
* Redirection, manipulation and extraction

---

### Navigating your file system

Find the current location by running the command `pwd` which stands for *print working directory*. At any give time, **current working directory** is current default directory.

```bash
pwd
# /home/UserName/
```

We can see what files and subdirectories are in this directory by running `ls`, which stands for "listing".

```bash
ls
```

Navigating to the `MGSS_Intro/` directory can be done with `cd` command which stand for *change directory*.

```bash
cd MGSS_Intro/
```

Run the `ls` command to list the content of current directory. Check whether there are two *.fastq* files.

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

Navigate to `backup/` directory and use `mv` command to rename and move `Test_1_backup.fastq` as `Test_1_copy.fastq` to the dirctory immediately above .

```bash
cd backup/

mv Test_1_backup.fastq ../Test_1_copy.fastq
```

Return the directory immediately above, check whether the `Test_1_copy.fastq` was moved and renamed as instructed and remove it by using the `rm` command.

```bash
cd ..

rm Test_1_copy.fastq
```

See whether you can remove the `backup/` directory by using `rm` command as well. 

```bash
rm backup/
# rm : can not remove 'backup/': Is a directory
```

By default, `rm` will not delete directories. This can be done by using `-r` (recursive) option.

```bash
rm -r backup
```

---

### Examining file contents

There are number of ways to examine the content of a file. `cat` and `less` are two commonly used programs for a quick look. Check the content of `Test_1.fastq` by using these commands. Take a note on the differences. 

```
cat Test_1.fastq
# less Test_1.fastq
```

A few useful shortcuts for navigating in `less`

![](https://github.com/GenomicsAotearoa/metagenomics_summer_school/blob/master/materials/figures/ex1_less_shortcuts.jpg) 

There are ways to take a look at part of a file. For an example `head` and `tail` command will scan the beginning and end of a file, respectively. 

```bash
head Test_1.fastq

tail Test_1.fastq
```

Adding `-n` option to either of these command will print the first or last *n* lines of a file.

```bash
head -n 1 Test_1.fastq
# @SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
```

---

### Redirection, manipulation and extraction

Although using `cat` and `less` commands will allow us to view the content of the whole file, most of the time we are in search of particualr characters (strings) of interest than the full content of the file. One of the most commonly used command-line utility to search for strigs is `grep`. Let's use this command to search for the string NNNNNNNNNN in `Test_2.fastq` file.

```bash
grep NNNNNNNNNN Test_2.fastq
```

Retrive and discuss the output you get when `grep` was executed with `-B1` and `-A1` flags.

```bash
grep -B1 -A2 NNNNNNNNNN Test_2.fastq
```

In both occassions, outputs were printed on the terminal where they can not be reproduced without the execution of the same command. In order for "string" of interest to be used for other operations, this has to be "redirected" (captured and written into a file). The command for redirecting output to a file is `>`. Redirecting the string of bad reads that was searched using the `grep` command to a file named `bad_reads.txt` can be done with

```bash
grep -B1 -A2 NNNNNNNNNN Test_2.fastq > bad_reads.txt
```

Use the `wc` command to count he number of words, lines and characters in the `bad_reads.txt` file.

```bash
wc bad_reads.txt
```

Add `-l` flag to `wc` command and compare the number with the above output

```bash
wc -l bad_reads.txt
```

In an instance where the same operatioin has to applied for multiple input files and the outputs are to be redirected to the same output file, it is important to make sure the the new output is not over-writitng the previous output. This can be avoided with use of `>>` (append redirect) command which will append the new output to the end of the file, rather than overwriting it. 

```bash
grep -B1 -A2 NNNNNNNNNN Test_1.fastq >> bad_reads.txt
```

Executing the same operation on multiple files with the same file extention (or different) can be done the use of *wildcards* which are symbols or special characters that represent other characters. For an example. Using `*` wildcard, we can run the previous `grep` command to both files at once

```bash
grep -B1 -A2 NNNNNNNNNN *.fastq >> bad_reads.txt
```

The only objective of above redirection example is to search for a string in a file(s), write the output to a file and then count the number of lines of the  output. Generating output files for short routine taks like that will end up generating an excessive number of files with little value. The `|` (pipe) command is a commonly used method to apply an operation for an ouput without creating intermediate files. It takes the output generated by one command and uses that output as input to antoher command. 

```bash
grep -B1 -A2 NNNNNNNNNN Test_2.fastq | wc -l
```

In an instance where we aren't after a particular string, but want  sections from each lines of file(s) to be extracted and written into an output file for further operations, `cut` command is a great utility. It can be used to cut parts of a line by **byte position, character and field**. Basically the cut command slices a line and extracts the text.




### Cut
Print selected parts of lines from each FILE to standard output
- When invoking cut, use the -b, -c, or -f option, but only one of them.

```bash
cut --help
```

**Create a file called names.txt and add random names** (first and last), numbers e.g. Ngoni Faya 19
```bash
cut -d " " -f 1 names.txt
cut -d " " -f 1-3 names.txt
cut -d " " -f 1,3 names.txt
```
The tab character is the default delimiter that cut uses to determine what constitutes a field. So, if your file's fields are already delimited by tabs, you don't need to specify a different delimiter character.



### basename
Basename is a function in UNIX that is helpful for removing a uniform part of a name from a list of files. In this case, we will use basename to remove the .fastq extension from the files that we’ve been working with.
```bash
basename Test_1.fastq .fastq
```

### sed
sed is a stream editor. A stream editor is used to perform basic text transformations on an input stream (a file, or input from a pipeline)

```bash
sed --help
```

