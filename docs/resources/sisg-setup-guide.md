# SISG Setup Guide: NeSI OnDemand and Slurm

This guide shows you how to log in to NeSI OnDemand, set up SSH access to the HPC login node, and submit a Slurm job to the SISG reservation.

For full details, see the NeSI support docs: <https://docs.nesi.org.nz/Interactive_Computing/OnDemand/how_to_guide/>

| Item | Value |
|---|---|
| Project code | `nesi02659` |
| Reservation | `MPI_mc062` |
| Login node | `login03.hpc.nesi.org.nz` |

---

## 1. First Access

### a. Quick Access

Go to <https://ondemand.nesi.org.nz/pun/sys/shell/ssh/login.hpc.nesi.org.nz> and follow all steps required to log in using your credentials.

If this does not work follow the step-by-step guide:

### b. Step-by-step Access Guide

#### 1. Log in to OnDemand

1. Go to <https://ondemand.nesi.org.nz/>.
2. Select your institution and sign in.
3. Enter the 6-digit code from your authenticator app. On your first login, OnDemand asks you to set this up.

See [Log in](https://docs.nesi.org.nz/Interactive_Computing/OnDemand/how_to_guide/) in the NeSI docs for screenshots.

#### 2. Open a terminal

On the OnDemand home page, select **Applications >> Mahuika Shell Access**.

---

## 2. Set up your SSH key

Do this once. In the terminal, run the commands below.

```bash
cd .ssh
ssh-keygen -t rsa -b 4096 -C "$USER@login03.hpc.nesi.org.nz"
cat id_rsa.pub >> authorized_keys
```

Press **Enter** at each `ssh-keygen` prompt to accept the defaults.

For example, if your user ID is `matt.bixley`, the second command is:

```bash
ssh-keygen -t rsa -b 4096 -C "matt.bixley@login03.hpc.nesi.org.nz"
```

---

## 3. Submit a Slurm job to the reservation

Create a file called `sisg-test.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=sisg-test
#SBATCH --account=nesi02659
#SBATCH --reservation=MPI_mc062
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out

echo "Running on $(hostname)"
```

Submit it:

```bash
sbatch sisg-test.sl
```

The reservation is "magnetic": jobs under `nesi02659` are sent to it automatically. Keep the `--reservation=MPI_mc062` line anyway so it's clear where your job runs.

To check the reservation is active:

```bash
scontrol show reservation MPI_mc062
```

---

## 4. Check your job

| Task | Command |
|---|---|
| See your jobs | `squeue --me` |
| Cancel a job | `scancel <jobid>` |
| Read the output | `cat sisg-test-<jobid>.out` |

For more on Slurm, see the NeSI support docs: <https://docs.nesi.org.nz/>
