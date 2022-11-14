# Working with job scheduler

<center>![image](../theme_images/scheduler_image.png){width="500"}</center>

An HPC system might have thousands of nodes and thousands of users. How do we decide who gets what and when? How do we ensure that a task is run with the resources it needs? This job is handled by a special piece of software called the scheduler. On an HPC system, the scheduler manages which jobs run where and when. In brief, scheduler is a 


* Mechanism to control access by many users to shared computing resources
* Queuing / scheduling system for users’ jobs
* Manages the reservation of resources and job execution on these resources 
* Allows users to “fire and forget” large, long calculations or many jobs (“production runs”)

!!! info "Why do we need a scheduler ?"

    * To ensure the machine is utilised as fully as possible
    * To ensure all users get a fair chance to use compute resources (demand usually exceeds supply)
    * To track usage - for accounting and budget control
    * To mediate access to other resources e.g. software licences

    **Commonly used schedulers**

    * Slurm
    * PBS , Torque
    * Grid Engine
    * LSF – IBM Systems

    <center>![image](../theme_images/slurm_comms2compute.png){width="800"}</center>

    <small>Researchers can not communicate directly to  Compute nodes from the login node. Only way to establish a connection OR send scripts to compute nodes is to use scheduler as the carrier/manager</small>


## Life cycle of a slurm job

<br>
<center>![image](../theme_images/batch_system_flow%20.png){width="800"}</center>
<br>



??? info "Commonly used Slurm commands"

    | Command        | Function                                                                                             |
    |:---------------|:------------------------------------------------------------------------------------------------------|
    | `sbatch`       | Submit non-interactive (batch) jobs to the scheduler                                                 |
    | `squeue`       | List jobs in the queue                                                                               |
    | `scancel`      | Cancel a job                                                                                         |
    | `sacct`        | Display accounting data for all jobs and job steps in the Slurm job accounting log or Slurm database|
    | `srun`         | Slurm directive for parallel computing                                                                      |
    | `sinfo`        | Query the current state of nodes                                                                     |
    | `salloc`       | Submit interactive jobs to the scheduler                                                             |


About

## Anatomy of a slurm script and submitting first slurm job 🧐

As with most other scheduler systems, job submission scripts in Slurm consist of a header section with the shell specification and options to the submission command (`sbatch` in this case) followed by the body of the script that actually runs the commands you want. In the header section, options to `sbatch` should be prepended with `#SBATCH`.

<br>
![image](../theme_images/anatomy_of_a_slurm_script.png){width="700"}
<br>

!!! quote ""

    Commented lines are ignored by the bash interpreter, but they are not ignored by slurm. The `#SBATCH` parameters are read by slurm when we submit the job. When the job starts, the bash interpreter will ignore all lines starting with `#`. This is very similar to the shebang mentioned earlier, when you run your script, the system looks at the `#!`, then uses the program at the subsequent path to interpret the script, in our case `/bin/bash` (the program `bash` found in the */bin* directory

---

??? info "Slurm variables"

    | header          | use                                 | description                                          |
    |:--------------- |:------------------------------------|:-----------------------------------------------------|
    |--job-name 	  | `#SBATCH --job-name=MyJob` 	        |The name that will appear when using squeue or sacct. |
    |--account 	      | `#SBATCH --account=nesi12345` 	    |The account your core hours will be 'charged' to.     |
    |--time 	      | `#SBATCH --time=DD-HH:MM:SS` 	    |Job max walltime.                                     |
    |--mem 	          | `#SBATCH --mem=512MB` 	            |Memory required per node.                             |
    |--cpus-per-task  | `#SBATCH --cpus-per-task=10` 	    |Will request 10 logical CPUs per task.                |
    |--output 	      | `#SBATCH --output=%j_output.out` 	|Path and name of standard output file. `%j` will be replaced by the job ID.         |
    |--mail-user 	  | `#SBATCH --mail-user=me23@gmail.com`|address to send mail notifications.                   |
    |--mail-type 	  | `#SBATCH --mail-type=ALL` 	        |Will send a mail notification at BEGIN END FAIL.      |
    |                 | `#SBATCH --mail-type=TIME_LIMIT_80` |Will send message at 80% walltime.                    |


??? question "Exercise"
    Copy the contents of the `BLAST/` folder to your current directory, using the following command

    ```bash
    cp -r /nesi/nobackup/nesi02659/SLURM/BLAST ./
    ```

    We will then navigate into this directory with the `cd` command, and inspect the text of the file `blast-test.sh` using `less` or `nano`.

    ```bash
    cd BLAST/

    less blast-test.sh
    ```

    Evaluate the contents of the `blast-test.sh` script. Take a note of the basic slurm variables, path variables, etc. We will revisit these in the afternoon, when you create your own slurm scripts.

    Submit the script to the job queue as below.

    ```
    sbatch blast-test.sh
    ```