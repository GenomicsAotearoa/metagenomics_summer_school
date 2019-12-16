# Useful command line shortcuts

---

### Overview

* Bash
* Nano
* SLURM commands
* SBATCH flags

---

### Bash commands

|Command|Function|Example|
|:---|:---|:---|
|`cat`|Read and print contents of a file|`cat my_file.txt`|
|`grep`|Search a file for lines containing text of interest|`grep "words" my_file.txt`|
|`sed`|Find and replace text in a file|`sed 's/words/pictures/' my_file.txt`|
|`head`|Print the first *n* lines from a file|`head -n 10 my_file.txt`|
|`tail`|Print the last *n* lines from a file|`tail -n 10 my_file.tx`t|
|`less`|Partially open a text file for viewing, readings lines as the user scrolls|`less my_file.txt`|
|`cut`|Retrieve particular columns of a file.<br>Uses tab as the delimiting character by default, can change with `-d` parameter|cut -f1,2,5 -d ',' my_file.txt`|
|`paste`|Opposite of `cut`, stitch together files in a row-wise manner|`paste my_file.txt my_other_file.txt`|

---

### Nano shortcuts

|Command|Function|
|:---|:---|
|<kbd>Ctrl</kbd> + <kbd>X</kbd>|Close the file|
|<kbd>Ctrl</kbd> + <kbd>O</kbd>|Save the file|
|<kbd>Ctrl</kbd> + <kbd>W</kbd>|Search for a string in the file|
|<kbd>Ctrl</kbd> + <kbd>Space</kbd>|Go forward one word in a line|
|<kbd>Alt</kbd> + <kbd>Space</kbd>|Go backwards one word in a line|
|<kbd>Ctrl</kbd> + <kbd>A</kbd>|Go to the beginning of a line|
|<kbd>Ctrl</kbd> + <kbd>E</kbd>|Go to the end of a line|
|<kbd>Ctrl</kbd> + <kbd>W</kbd>, <kbd>Ctrl</kbd> + <kbd>V</kbd>|Go to the last line of the file|
|<kbd>Ctrl</kbd> + <kbd>K</kbd>|Cut the current line and save it to the clipboard|
|<kbd>Ctrl</kbd> + <kbd>U</kbd>|Paste (**_u_**n-cut) the text in the clipboard|
|<kbd>Ctrl</kbd> + <kbd>Shift</kbd> + <kbd>6</kbd>|Select text|

---

### SLURM commands

|Command|Function|Example|
|:---|:---|:---|
|`sbatch`|Submits a SLURM script to the queue manager|`sbatch job.sl`|
|`squeue`|Display jobs in the queue<br>Display jobs belonging to user *usr9999*<br>Display jobs on the *long* partition|`squeue`<br>`squeue -u usr9999`<br>`squeue -p long`|
|`sacct`|Displays all the jobs run by you that day<br>Display job 123456789|`sacct`<br>`sacct -j 123456789`|
|`scancel`|Cancel a queued or running job|`scancel 123456789`<br>`scancel -u usr9999`|
|`sshare`|Shows the [Fair Share](https://support.nesi.org.nz/hc/en-gb/articles/360000743536-Fair-Share) scores for all projects of which *usr9999* is a member|`sshare -U usr9999`|
|`sinfo`|Shows the current state of the SLURM partitions|`sinfo`|
|`sview`|Launches a GUI for monitoring SLURM jobs|`sview`|

---
