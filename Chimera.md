# MFF HPC metacenter

Clusters are now part of the MFF HPC metacenter, which is [SLURM in a multi-cluster setup](https://slurm.schedmd.com/multi_cluster.html).
The clusters are managed by the SLURM resource manager.

### Table of contents
 - [Chimera cluster](#chimera-cluster)
 - [KSI/DSE clusters](#ksidse-clusters)
 - [Graphical interface for cluster monitoring](#graphical-interface-for-cluster-monitoring)
 - [Authentication](#authentication)
 - [Paths](#paths)
 - [SLURM Crash Course](#slurm-crash-course)
     - [Useful examples](#useful-examples)
	 - [Job resource usage](#job-resource-usage)
 - [Description of available clusters](#description-of-available-clusters)
      - [Parlab cluster specification](#parlab-cluster-specification)
	  - [Gnulab cluster specification](#gnulab-cluster-specification)
	  - [Chimera cluster specification](#chimera-cluster-specification)
    		- [Using MPI](#mpi)
		- [Checkpointing](#checkpointing)
		- [Mathematica](#mathematica)
 - [Charliecloud](#charliecloud)
 - [Apptainer](#apptainer)
 - [Choosing an interface](#choosing-an-interface)
 - [Installed software on clusters](#installed-software-on-clusters)
 - [Accounts management](#accounts-management) 


## Chimera cluster

Faculty cluster of MFF UK. Mostly composed of computers from physics departments. Used for both research and teaching.

## KSI/DSE clusters

We have two computational clusters **parlab** and **gpulab**.
They are used for exercises, homework assignments, and research, where you need high-performance computing.

> Do not execute any of your code on the front-end nodes directly, use worker nodes instead!

All unknown terms (front-end node, worker node) will be explained later. Front-end nodes do not have installed any development software.

## Graphical interface for cluster monitoring

The cluster usage and load can be monitored conveniently in real time using the `sview` tool, see the figure below. In the left vertical panel the little boxes represent individual nodes and their colour their current load and availability. Clicking on the box will open a window containing full info about the node's hardware and its features that can be used in SLURM to select nodes with a particular hardware to run your job.

Toggling between info for any of the three clusters is achieved by using the drop-down menu located on the top right part of the sview window.

Sview can be used to monitor running jobs, view the partitions and node info. For example node-oriented information can be added to the default set of tabs as a new tab using the tab `Visible Tabs`. The resulting view is exemplified on the figure below.

To run `sview` you need to enable X11 forwarding by adding the `-X` option to your ssh command which you use to connect to the cluster.

![Fig.1](sview_example.png "SVIEW tool for cluster monitoring")

Only on Chimera there is also the command line tool `freenodes` that can be used to view activity on the cluster. It was developed for the use on cluster *Snehurka* - mathematical section of MFF UK. For details please visit the help page [there](https://cluster.karlin.mff.cuni.cz/freenodes/). It displays both node-oriented and partition-oriented load statistics for the cluster. It is particularly useful for a quick assesment of which partitions are most available.

Our cluster is little bigger and so is the output of this command. Here is example of the output:

![Fig.2](freenodes.png "Example of output of freenodes on chimera cluser")

## Authentication

Authentication is now done using the MFF LDAP server, which is a selection of appropriate accounts from [CAS UK](https://idp.cuni.cz/cas/login).
Use the same username and password to log in as you use to log in to SIS. To run jobs on chimera, you need associated account - you can login to [JupyterHUB](https://hpc.troja.mff.cuni.cz:8000/) and start session, or if you want/need specific account association send an email that includes including your CAS username and affiliation in the message to jiri.eliasek@matfyz.cuni.cz or your account coordinator.

> You need to change the password at least once at [CAS UK](https://idp.cuni.cz/cas/login) for this to work. The default password received with a fresh MFF account does not work for reasons outside of our control.

## Paths

The home directory is now automatically created during the first login.
The path to the home directory is `/home/<CASlogin>`.
You can ask the administrators to create a personal directory for big data or a directory for shared project files in `/work` directory.
The `/home` and `/work` directories are on the local disk array for each cluster.


To facilitate data sharing between clusters in the metacenter, the `/home` and `/work` directories are shared with other clusters
using the tunneled NFS protocol and are available on fixed paths according to the cluster location,
e.g. the `/home` directory of the Chimera cluster located in Troja has the path `/troja/home` on all clusters (including the Chimera cluster).

> Use these shared directories with caution, as they could overwhelm the network connection for that location.

**Copy your frequently used data to the local cluster directories (`/home` or `/work`) where the computation will take place.**

#### Temporary directories on chimera cluster
There are multiple choices for temporary files (typically needed for the duration of your job) for you jobs:
 - `$TMPDIR` - It is set up automatically by slurm on tmpfs (in RAM). Files there are cleaned up after job is finished. **Data there are counted to the memory limit of the job.**
    - You can use this for setup of your job - i.e. copy / compile binaries and run your job there and copy results to persistent place after finishing - this automatically cleans your environment and is easy solution for cleaning preempted jobs that should be restarted from the start
	- Fast temporary space if you need save small amount of data during running task

 - `/scratch` -  This is a fast network-attached storage shared across all nodes. Only nodes with Infiniband take full advantage of the speed.
 If your temporary files are large and do not fit memory (`$TMPDIR`) alongside your running task, please create your own directory in `/scratch/tmp` and cleanup after yourself.
    - You should use this, if your temporary data are large - i.e. you would have to request too much memory or it will not fit into memory with your task.
	- You can use it as a temporary storage between multiple jobs, but if you want to use this for longer term storage, plese contact administrator.

 - `/scratch_local` - *This option is not recommended.* 
 Few nodes have also local scratch.  It is on nodes with `local_scratch` feature. You can also request these by specifying minimal amount of space for local scratch using sbatch option `--tmp=<size>[units]`. It is mounted to `/scratch_local` - **The directory always exist, but is mounted only on nodes with `local_scratch` feature.** You have to clean up after yourself and be aware that slurm does not check how much of the scratch is already in use, only maximum capacity specified in configuration so the amount of space requested using `--tmp` option is not guaranteed. 
     - You should not use these, unless you know what you are doing. Only some older nodes have local scratch and its usage have been mostly moved to `/scratch`.

## SLURM Crash Course

All informations about SLURM can be found on its
[SLURM documentation](https://slurm.schedmd.com/documentation.html)
page or on [SLURM tutorials](https://slurm.schedmd.com/tutorials.html) page.
Anyway, we have provided a short description of SLURM and of our clusters.

### Terminology

A **cluster** is a bunch of **nodes**.
Nodes are grouped together to **partitions**.
Partitions may overlap, ie. one node can be in more partitions.
**Feature** is a string describing a feature of a node, e.g. `avx2` for a node capable of executing AVX2 instructions.
Each node has two sets of features: current features and available features.
Usually, they are same.
But in some cases, the node is capable of changing current features set on demand.

SLURM manages **resources** and schedules jobs according to available resources.
The most important resources are CPUs and memory, but it is able to manage other generic resources (GRES), like GPUs, etc.
Moreover, the time is resource as well.

A **user** is identified by his/her login.
**Account** is a billing entity (well, we won't charge you for using our clusters).
Each user must have assigned an account. Moreover, user can be assigned to more accounts and use them depending on what he/she is doing.
Accounts can only be allowed access to certain partitions.

A user can launch a **job**.
Job has a state, has some reserved and assigned resources, and returns an error code after completion.
Job consists from **steps** (usually 1 step).
Each step is executed by a bunch of **tasks** (usually 1 task), where resources are equally assigned to the tasks of a step.
Jobs are inserted to a **scheduling queue**, where you can find them.

Partition has a priority.
A job submitted to a partition with higher **priority** can preempt an another job submitted to a partition with lower priority.
Our clusters use only two kind of preemption:

- **SUSPEND** - the preempted job is suspended and releases CPUs. Unfortunately, it does not release memory or GPUs.
- **REQUEUE** - the preempted job is killed and requeued to the scheduling queue, waiting for available resources.
All resources including GPUs are released.
Jobs submitted to the partition should use some form of checkpointing, otherwise all work can be lost, when the job is killed during preemption.

### Important commands

There are many commands (see
[SLURM man pages](https://slurm.schedmd.com/man_index.html) or [SLURM command summary](https://slurm.schedmd.com/pdfs/summary.pdf)
). The most important commands are:


- **srun** - reserves and assigns resources and runs an interactive job, use it only for short or debug jobs
- **sbatch** - reserves and assigns resources and submit a batch as a job, use it for long jobs (>1 hour)
- **salloc** - only reserves and does not assign resources, use srun or sbatch for assigning resources subsequently
- **sinfo** - get info about cluster state
- **scancel** - cancels a job
- **squeue** - view info about scheduling queue
- **sbcast** - distributes a file to the nodes allocated to a job
- **sstat** - display info about a job
- **seff** - displays statistics of resources for finnished job

> On **Chimera** there is also script *smypartition* - display your accounts and partitions you can use (ALL means any account can be used)

Job submission commands (srun, sbatch, salloc) have a common set of important options:

| Shrt | Long                 | Description |
| ---- | -------------------- |------------ |
| `-A` | `--acount=` 	      | Charge resources used by the job to the specified account. If not specified, user's default account is charged |
| `-B` | `--extra-node-info=` | Select only nodes with at least specified number of sockets, cores per socket, and threads per core. This option does not specify the resource allocation, its just constraint. |
| `-C` | `--constraint=`      | Select only nodes with matching features |
| `-c` | `--cpus-per-task=`   | Number of CPUs per 1 task |
| `-e` | `--error=`           | Standard error stream is redirected to the specified file |
| `-G` | `--gpus=`            | Specifies the number of GPUs required for the job in the form `[type:]count`. It is a shortcut for --gres=gpu:type:count. |
|      | `--gres=`            | Specifies comma delimited list of GRES. Each entry on the list is in form `name[[:type]:count]` |
| `-i` | `--input=`           | Standard input stream is redirected from the specified file |
| `-J` | `--job-name=`        | Job name |
| `-L` | `--licenses=`        | Specifies comma delimited list of licenses allocated to the job. Each entry on the list is in form `name[:count]` |
| `-M` | `--clusters=`        | Clusters to issue commands to (comma separated list or `all`) |
| `-m` | `--distribution=`    | Select distribution method for tasks and resources. For more info see documentation |
|      | `--mem=`             | Specify the real memory required per node, default unit is MB. To allocate all available memory on the node use `--mem 0` |
|      | `--mem-per-cpu=`     | Specify the memory required per allocated CPU |
|      | `--mem-bind=`        | Specify the memory binding. For more info see documentation |
| `-N` | `--nodes=`           | A minimum of allocated nodes. Default is 1 |
| `-n` | `--ntasks=`          | Number of tasks. Default is 1 |
| `-o` | `--output=`          | Standard output stream is redirected to the specified file |
| `-p` | `--partition=`       | Request a specific partition for the resource allocation. If not specified, default partition os chosen |
| `-t` | `--time=`            | Set a limit on the total run time of a job. If not specified, default time for a selected partition is used. Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours". |
|      | `--threads-per-core=`| If set to 1, will turn off multithreading which is usefull for many applications; simmilar to `--hint=nomultithread` |
| `-a` | `--array=`           | Specify indexes for job array; runs multiple of the same jobs (with the reqirements set for each job). Only difference is the environment variable SLURM_ARRAY_JOB_ID is set to one of the indexes |

### Useful examples

`srun --pty bash`

Stars interactive job (with bash as shell).

> Interactive jobs do not work across clusters, please log into the correct frontnode to run interactive jobs.

`salloc`

Starts a shell in a minimal environment (resources). Useful for making small builds or debugging in a restricted environment.

`salloc -N2 -n4`

Starts a shell in a minimal environment with 4 tasks on 2 nodes from which you can run srun and see all outputs. Useful as an interactive debugging tool for MPI jobs.

`sbatch -p gpu-long --gpus=1 mygpucode.sh`

Starts a batch job which will have one NVIDIA V100 card available.

`sinfo -o "%P %L %l"`

Prints info about default and maximum job time for partitions in a cluster.

`sinfo -o "%n %b %f"`

Prints info about current and available feature set for nodes in a cluster.

`sinfo -o "%.5P %.5a %.10l %.16F %.10p %N"`

Prints partition summary with priority.

`sinfo -Nel`

Prints nodes list with resources for every partition.

`srun -p mpi-homo-short -n 128 -N 8 --mem-per-cpu=2G mympijob`

Starts a MPI job with 128 tasks/ranks spanning over 8 nodes in the homogeneous partition assigning 1 CPU and 2 GB RAM to each task.

#### Examples of script headers

These are examples of script header with *Slurm* options, that can be run by `sbatch script.sh`. You can use this format also  with options mentioned above.

#### Standard single-node single-thread job

This is the simple example of single thread jobs. You need only specify how long it takes (otherwise it useses default value, depending on the partition) and memory consumption. If you run an array jobs, these requirements are for each job in the array.
```
#!/bin/bash
#SBATCH --partition=ffa           # partition you want to run job in
#SBATCH --time=00:15:00           # walltime for the job in format (days-)hours:minutes:seconds
#SBATCH --mem=24000               # memory resource
#SBATCH --job-name="job_name"     # change to your job name
#SBATCH --output=output.txt       # stdout and stderr output file
#SBATCH --mail-user=user@example.com # send email when job changes state to email address user@example.com


# some real work here
```


##### Standard multi-node job

```
#!/bin/bash
#SBATCH --partition=ffa           # partition you want to run job in
#SBATCH --time=00:15:00           # walltime for the job in format (days-)hours:minutes:seconds
#SBATCH --nodes=2                 # number of nodes (can be only 1)
#SBATCH --ntasks-per-node=4       # processes per node
#SBATCH --mem=24000               # memory resource per node
#SBATCH --job-name="job_name"     # change to your job name
#SBATCH --output=output.txt       # stdout and stderr output file
#SBATCH --mail-user=user@example.com # send email when job changes state to email address user@example.com
#SBATCH --exclusive               # Use whole node
#SBATCH --signal SIGTERM@60       # send signal SIGTERM 60s before ending job (preemption, time limit)


# some real work here
```

This will run job in ffa partition with length of 15 minutes on 2 nodes, 4 tasks per node with 24000MB of RAM on each node - and nothing else will be run on these nodes.
Standard output and error will be redirected into output.txt. Slurm will send SIGTEM to the job 60s before walltime ends, or in case of *gang* preemption. And email will be sent to the user@example.com if job changes begins, ends, or fails. (this can be changed by `--mail-type` option.

##### Job-array: submitting many jobs at the same time

```
#!/bin/bash
#SBATCH --partition=ffa           # partition you want to run job in
#SBATCH -A project                # use account project instead of default
#SBATCH --time=12:00:00         # walltime for the job in format (days-)hours:minutes:seconds
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=1       # processes per node
#SBATCH --cpus-per-task=8         # cpus per tasks
#SBATCH --mem-per-cpu=2G          # memory resource per cpu
#SBATCH --array 0-127             # duplicate the job 128 times - slurm sets SLURM_ARRAY_TASK_ID environment variable which corresponds to the index each job in the array
#SBATCH -d afterok:10:12?after:13 # run this job after jobs with ID 10 and 12 ends successfully and job 13 starts executing (this line is optional and only an example of course)
#SBATCH --job-name="job_name"     # change to your job name
#SBATCH --error=error\_%a.txt         # redirects stderr to error\_${SLURM_ARRAY_TASK_ID}.txt
#SBATCH --output=output\_%a.txt       # redirects stdout to output\_${SLURM_ARRAY_TASK_ID}.txt
#SBATCH --mail-user=user@example.com # send email when job changes state to email address user@example.com


# the environment variable $SLURM_ARRAY_JOB_ID is set to the ID of the first job in the array
# some real work here dependent on $SLURM_ARRAY_TASK_ID
```

Runs 128 jobs, every with 8 CPUs and 2GB of RAM per CPU. All is run in partition ffa and with account "project". Max time is 12 hours for every job. Jobs start after job 10 and 12 ends successfully and job 13 start to run.
Here the standard error and standard output are redirected to separate files.
Emails in this case are sent about array as a whole.

### Recommendations

- Prefer using `sbatch` for long jobs.
Moreover `sbatch` allows setting of SLURM parameters in the executed shell-script,
you don't need to write them always on the command-line.
Moreover, jobs executed by `sbatch` can be requeued, when cancelled (e.g. for priority reasons).

- Don't forget to request GPU using `--gpus=1`, when running on gpu-xxx partitions.

- Do not use `srun --pty bash` followed by executing a long computation from the command-line.
When the job finishes, you will block resources, until it either timeouts or you will exit the bash.
Again, prefer using `sbatch` for long jobs.

- Set mail-user and mail-type option in `sbatch` shell using `#SBATCH --mail-user=mymail@isp.edu` and `#SBATCH --mail-type=END,FAIL`.
SLURM will send you an email, when the job finishes. There are many other mail-type variants.

#### Job resource usage

It is useful to check if your jobs are using assigned resources. It this section we will notice some linux and slurm tools to help you find how much resources your jobs are using and what may indicate there is a problem. 

- you can use `time ./program` (or `/usr/bin/time -v ./program` for more details) to get some statistics about used cpu time (and more) that the executable `./program` used. 
  - this is helpful to discover problems when your program spends too much time in `system`.
- you can check running jobs by simply ssh to the node your job is running on and 
 - use basic tools like `top` or `htop` - you can use option `-u$USER` to filter your processes.
 - you can use `sstat job_id` to get information about running jobs - you can `sstat job_id | less -S` for better viewing experience; some useful fields:
   - *MaxRSS* - maximum resident set size (memory footprint) of all tasks in job
	 - if this is close to memory allocated for task in slurm, there might be problems - processes killed for trying to use more memory, caching problems, ...
   - MinCPU - Minimum (system + user) CPU time of all tasks in job.
	 - if this is considerably smaller then job (wall time) * #CPUs, there is probably some problem (often waiting for I/O)


 - after the job finishes you can ask use tools like `seff job_id` to get resource usage; or you can use `sacct -j job_id` to get even more information about job
    - the `seff` command presents statistics in user friendly way
	- CPU efficiency around 50% is usually expected on nodes with SMT
	- if CPU efficiency is low (< 30%), there is probably some issue that can be fixed
	- if CPU efficiency is close to 100% and you are running on node with SMT this may indicate some issues with running too many threads and using tools like `time` would be helpful to check if there is a problem

 >If the CPU efficiency is low, there is usually some kind of problem - often it is stuck on I/O operations - if this happens when running larger batch of jobs, that read large data, please consider to run smaller number of the jobs simultaneously (consider using arrays with % operator) - this also often hinders other users access to storage

Examples of `seff` usage:
 - here is example of job with CPU Efficiency close to 100%
```
$ seff  1111
Job ID: 1111
Cluster: chimera
User/Group: user1/user1
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:09:13
CPU Efficiency: 96.68% of 00:09:32 core-walltime
Job Wall-clock time: 00:02:23
Memory Utilized: 14.44 MB
Memory Efficiency: 0.35% of 4.00 GB
```
 - here is example of job with low CPU Efficiency; here the job was indeed limited by I/O operation:
```
$ seff  123 
Job ID: 123
Cluster: chimera
User/Group: user2/user2
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 2
CPU Utilized: 00:00:47
CPU Efficiency: 1.83% of 00:42:54 core-walltime
Job Wall-clock time: 00:21:27
Memory Utilized: 6.01 GB
Memory Efficiency: 42.92% of 14.00 GB
```
 

## Description of available clusters

**Chimera** cluster can be accessed using SSH. Front-end server **hpc.troja.mff.cuni.cz**

> fingerprints:
>
> 3072 SHA256:C+WGqtnspLv9GQDMkwyDF3UHr6k2M3oOzruxwD5DSrI localhost (RSA)
>
> 256 SHA256:F/rlEJM7MEF4n/NBpGnrjgzIxZizwRIbRw3o7dEtm5Y localhost (ECDSA)
>
> 256 SHA256:5PkELOKy8uy8s/PnwqMwe9rW9A/CUVLmouqlZBmPRYY localhost (ED25519)
>
> 3072 MD5:14:d3:c5:56:63:54:57:12:90:2b:3f:fe:bc:28:7f:ad localhost (RSA)
>
> 256 MD5:f6:8e:21:77:da:8b:33:2d:5a:a7:1b:cf:da:20:21:84 localhost (ECDSA)
>
> 256 MD5:a6:ca:5f:3e:ca:e5:31:fc:83:b3:b4:db:75:ed:20:e3 localhost (ED25519)

As mentioned above, KSI/DSE have two clusters **parlab** and **gpulab**.
Access to the cluster is always through the front-end server using SSH on port 42222.
Front-end servers have the same name as the cluster,

i.e. **parlab.ms.mff.cuni.cz**
> <span style="color:red">new fingerprints (system reinstalled)</span>
>
> 256 SHA256:xp8egf//QRHd09JzUjA+MAar5HAuku1UMtORBSgnrRQ [parlab.ms.mff.cuni.cz]:42222 (ECDSA)
>
> 3072 SHA256:ZQ+h0UOulcrc9m6FYnfG71UEAXd8u1VPayPu5KdsvzQ [parlab.ms.mff.cuni.cz]:42222 (RSA)
>
> 256 SHA256:QPF6hDJJ6sXSMtLMeVbU6F90qYSFn6519I3RavvcFoE [parlab.ms.mff.cuni.cz]:42222 (ED25519)

> <span style="color:blue">old fingerprints (not valid anymore, just for checking)</span>
>
> 256 SHA256:oUhUEdCKxIxVphDzbiCcln/xwSFDJUKAXVEc/L2ltS0 [parlab.ms.mff.cuni.cz]:42222 (ECDSA)
>
> 3072 SHA256:WadA7b4labjw5VI+S8WoEAxmniHUYAJNVMck4tOZY0s [parlab.ms.mff.cuni.cz]:42222 (RSA)
>
> 256 SHA256:2TcysbzcJ600W43X6xHc4tlF76khxNkPfu5J19YYLPg [parlab.ms.mff.cuni.cz]:42222 (ED25519)

and **gpulab.ms.mff.cuni.cz**
> fingerprints

> 3072 SHA256:I7RfaXDgCK7u1Ucfy8JqNYAdd85IxbTVbm+WE2YYZ/0 [gpulab.ms.mff.cuni.cz]:42222 (RSA)

> 256 SHA256:spPrDC+cAvJ//jOzXApk4QuBMdvIwJ0h3iFZj6v/mCY [gpulab.ms.mff.cuni.cz]:42222 (ECDSA)

> 256 SHA256:63H1C+6PzFhmm/Sr+8Mhh6F/emtqG2mL2wHZjQoKz+U [gpulab.ms.mff.cuni.cz]:42222 (ED25519)

#### Course users

Student logins will be added to the appropriate course account after the first week of classes.

#### Other users

Use your CAS login. Lower case username is required for ssh login to work correctly.

#### For everyone

Faculty clusters use one external read-only MFF LDAP server, so you will not be able to change your password, use CAS instead.
Each course taught on our clusters has its account.
Any research group or project have their account as well.
Logins are assigned to the corresponding account depending on visiting relevant courses or working in a research/project group.

Faculty clusters have access to the same disk array using NFS.
You may find your home mounted on `/home`.
Moreover, users or projects can have an additional space mounted on `/work`.

You may use an environment variable **$TMPDIR** set to a private local temporary directory for a job.
It is created on every node allocated to the job before the job starts and it is completely removed after the job finishes.
The temporary directory can by used as a scratchpad.
 - On chimera cluster it resides in RAM and there are multiple [scratch](#temporary-directories-on-chimera-cluster) options. **Be aware that it is counted to your memory resource.**

 - On Parlab/Gpulab `$TMPDIR` is on a local SSD RAID. 
 - In both cases it is faster to access data here then accessing data on a remote storage.
 - On the other hand, **the space is limited, be careful!**


> **NOTE:** *Cores* in next tables are per *socket*; also *HT* are per *core* - you can also use `sinfo -N -O NodeList,SocketCoreThread,Memory:.10,Features:.50 | uniq` to get information about sockets, cores, SMT, memory, and features.

### Parlab cluster specification

All nodes are interconnected by InfiniBand FDR (56 Gb/s) for high-performance messaging using MPI.
Moreover, they are interconnected by 10 GbE for all other traffic.
The front-end server is connected by 10 GbE to the external world.
The latest version of OpenMPI is installed on all nodes.

Parlab nodes

| Node names | CPU                  | Sockets | Cores | HT | RAM    | GRES | Additional info |
| ---------- | -------------------- | ------- | ----- | -- | ------ | ---- | --------------- |
| w[201-208] | Intel Xeon Gold 6130 | 2       | 16    | 2  | 128 GB | | |

Parlab partitions

| Name             | Nodes      | Priority | Timelimit | Preemption | Allowed accounts               | Intended use |
| ---------------- | ---------- | -------- | --------- | ---------- | ------------------------------ | ------------ |
| ffa              | all        | 10       | 12 hours  | REQUEUE    | any                                                 | default, free-for-all, any job |
| mpi-homo-ffa     | w[201-208] | 30       | 12 hours  | REQUEUE    | any                                                 | free-for-all, MPI jobs on homogeneous nodes |
| mpi-homo-long    | w[201-208] | 100      | 7 days    | REQUEUE    | ksi, kdss, kdsstudent, nprg042s, nprg058s, nprg054s | debugging on newer CPUs, MPI debugging, long jobs |
| mpi-homo-short   | w[201-208] | 150      | 2 hours   | REQUEUE    | ksi, kdss, kdsstudent, nprg042s, nprg058s, nprg054s | executing short jobs on 2-socket system, MPI jobs on homogeneous nodes |

### Gpulab cluster specification

All nodes are interconnected by 10 GbE. The front-end server is connected by 10 GbE to the external world.

Gpulab nodes (integrated in chimera with their own environment)

| Node names   | CPU                    | Sockets | Cores | HT | RAM     | GRES | Additional info |
| ------------ | ---------------------- | ------- | ----- | -- | ------- | ---- | --------------- |
| volta01      | Intel Xeon Silver 4110 | 2       | 8     | 2  | 256 GB  | gpu:V100:2, vram:16G | 1x NVIDIA V100 PCIe 16 GB, latest V100 compatible CUDA (12.9.1) |
| volta02      | Intel Xeon Silver 4110 | 2       | 8     | 2  | 192 GB  | gpu:V100:2, vram:16G | 2x NVIDIA V100 PCIe 16 GB, latest V100 compatible CUDA (12.9.1) |
| volta03      | Intel Xeon Silver 4110 | 2       | 8     | 2  | 192 GB  | gpu:V100:2, vram:16G | 2x NVIDIA V100 PCIe 16 GB, latest V100 compatible CUDA (12.9.1) |
| volta04      | Intel Xeon Gold 5218   | 2       | 16    | 2  | 384 GB  | gpu:V100:1, vram:32G | 1x NVIDIA V100 SXM2 32 GB, latest V100 compatible CUDA (12.9.1) |
| volta05      | Intel Xeon Gold 5218   | 2       | 16    | 2  | 384 GB  | gpu:V100:4, vram:32G | 4x NVIDIA V100 SXM2 32 GB, latest V100 compatible CUDA (12.9.1) |
| ampere01     | AMD EPYC 9554P         | 1       | 64    | 2  | 384 GB  | gpu:L40:4, vram:48G  | 4x NVIDIA L40 PCIe 48 GB, latest CUDA  |
| ampere02     | AMD EPYC 9374F         | 2       | 32    | 2  | 768 GB  | gpu:A100:2, vram:80G | 2x NVIDIA A100 PCIe 80 GB, latest CUDA |
| ampere03     | AMD EPYC 9554          | 1       | 64    | 2  | 384 GB  | gpu:L40:1, vram:40G  | 1x NVIDIA L40 PCIe 40 GB, latest CUDA  |
| hopper01     | AMD EPYC 9554          | 2       | 48    | 2  | 1024 GB | gpu:H100:2, vram:80G | 2x NVIDIA H100 PCIe 80 GB, latest CUDA |
| bw01         | AMD EPYC 9454          | 1       | 48    | 1  | 768 GB  | gpu:H100:2, vram:80G | 2x NVIDIA RTX PRO Blackwell Server Edition PCIe 96 GB, latest CUDA |

Gpulab partitions (inside chimera)

| Name             | Nodes                                 | Priority  | Timelimit | Preemption | Allowed accounts               | Intended use |
| ---------------- | ------------------------------------- | --------- | --------- | ---------- | ------------------------------ | ------------ |
| ffa-gpulab       | all                                   | 10        | 12 hours  | REQUEUE    | any                            | default, free-for-all, any job |
| gpu-ffa-gpulab   | volta[01-05],ampere[01-03],hopper[01] | 20        | 12 hours  | REQUEUE    | any                            | free-for-all, any GPU job |
| gpu-long         | volta[01-05],ampere[01-02],hopper[01] | 100       | 7 days    | REQUEUE    | kdss, kdsstudent, ksviprj      | long GPU jobs |
| gpu-short        | volta[01-05],ampere[01-02],hopper[01] | 150       | 2 hours   | REQUEUE    | kdss, kdsstudent, ksviprj      | debugging GPU tasks, executing short GPU jobs |
| gpu-long-ksi     | volta[01-05],ampere[01-02],hopper[01] | 200       | 7 days    | REQUEUE    | ksi, ksiprj, ksibio            | long GPU jobs for KSI department |
| gpu-short-ksi    | volta[01-05],ampere[01-02],hopper[01] | 250       | 2 hours   | REQUEUE    | ksi, ksiprj, ksibio            | debugging GPU tasks, executing short GPU jobs |
| gpu-short-teach  | volta[01-05],ampere[01-02],hopper[01] | 300       | 2 hours   | REQUEUE    | nprg042s, nprg042t, nprg058s, nprg058t | parallel computing courses |
| gpu-bio          | ampere[01-03],hopper[01]              | 800       | 7 days    | REQUEUE    | ksibio                         | GPU jobs for Bioinformatics study program |

### Chimera cluster specification

All nodes are interconnected by > 1Gb Ethernet, nodes with feature IB are also interconnected with HDR infiniband (and also have 10Gb Ethernet).

Chimera nodes


| Node names      | CPU                            | Sockets | Cores | HT | RAM    | GRES | Additional info   |
| --------------- | ------------------------------ | ------- | ----- | -- | ------:| ---- | ----------------- |
| deepth1[3-4]    | AMD EPYC 7282                  | 2       | 16    | 2  | 128 GB |      | IB                |
| utf-asus[1-3]   | AMD EPYC 7302                  | 2       | 16    | 2  | 256 GB |      | local_scratch, IB |
| utf-hm1         | AMD EPYC 7763                  | 2       | 64    | 2  | 2 TB   |      | local_scratch, swap, IB |
| node[0-3,6]     | Intel(R) Xeon(R) E5-2630 v4    | 2       | 10    | 2  | 256 GB |      | local_scratch     |
| ucjf-asus1  | AMD EPYC 7532                      | 2       | 32    | 2  | 512 GB |      | IB                |
| ucjf-asus1  | AMD EPYC 7532                      | 2       | 32    | 2  | 512 GB |      | local_scratch, IB |
| ucjf-asusb[1-2] | AMD EPYC 7352                  | 2       | 24    | 2  | 256 GB |      | local_scratch, IB |
| onos2           | Intel(R) Xeon(R) CPU E5-2670 v3| 2       | 12    | 2  | 128 GB |      |                   |
| onos[3-4]       | Intel(R) Xeon(R) CPU E5-2670 v3| 2       | 12    | 2  | 256 GB |      |                   |
| c10t4[0-1]      | Intel(R) Xeon(R) CPU E5-2640 v3| 2       | 8     | 2  | 64 GB  |      |                   |
| c11t4[2-3]      | Intel(R) Xeon(R) CPU E5-2640 v3| 2       | 8     | 2  | 64 GB  |      |                   |
| c14t3[2-3]      | Intel(R) Xeon(R) CPU E5-2640 v3| 2       | 8     | 2  | 64 GB  |      |                   |
| c23t46          | Intel(R) Xeon(R) CPU E5-2630 v3| 2       | 10    | 2  | 64 GB  |      |                   |
| c74t57          | AMD EPYC 7451                  | 2       | 24    | 1  | 128 GB |gpu:quadro_m2000 vram:4G| NVIDIA Quadro M2000 4 GB |
| c71t55,c70t54,c72t56,c73t38 | AMD EPYC 7451                  | 2       | 24    | 1  | 128 GB |      |                   |
| c73t38          | AMD EPYC 7451                  | 2       | 24    | 1  | 128 GB |      |                   |
| c42t47          | Intel(R) Xeon(R) CPU E5-2650 v4| 2       | 12    | 2  | 64 GB  |      |                   |
| c41t19          | Intel(R) Xeon(R) CPU E5-2650 v4| 2       | 12    | 1  | 64 GB  |      |                   |
| c40t17,c43t48   | Intel(R) Xeon(R) CPU E5-2630 v4| 2       | 10    | 1  | 64 GB  |      |                   |
| c27t15,c28t14   | Intel(R) Xeon(R) CPU E5-2630 v4| 2       | 10    | 2  | 64 GB  |      |                   |
| c67t02          | Intel(R) Xeon(R) CPU E5-2630 v4| 2       | 10    | 1  | 64 GB  |gpu:quadro_m2000 vram:4G| NVIDIA Quadro M2000 4 GB |
| c64t03,c39t09   | Intel(R) Xeon(R) CPU E5-2630 v4| 2       | 10    | 2  | 64 GB  |      |                   |
| c65t04,c66t05   | Intel(R) Xeon(R) CPU E5-2650 v4| 2       | 12    | 2  | 64 GB  |      |                   |
| c77t06,c78t10,c75t34 | AMD EPYC 7532             | 2       | 32    | 1  | 256 GB |      |                   |
| c93t33          | AMD EPYC 7543                  | 2       | 32    | 1  | 512 GB |      | 16 numas          |
| c94t26          | AMD EPYC 7543                  | 2       | 32    | 1  | 256 GB |      | 16 numas          |
| c08t08          | Intel(R) Xeon(R) CPU E5-2650 v2| 2       |  8    | 1  | 64 GB  |      | no AVX2, FMA      |
| mff-a2-[01-02]  | AMD EPYC 7543                  | 2       | 32    | 2  | 512 GB |gpu:T600:1 vram:4G | IB, NVIDIA T600 4 GB |
| *-a2-[03-18]    | AMD EPYC 7543                  | 2       | 32    | 2  | 512 GB |      | IB                |
| mff-a2-19       | AMD EPYC 7543                  | 2       | 32    | 2  | 512 GB |gpu:L40:1 mps:100 vram:48G| IB, NVIDIA Tesla L40 48 GB |
| ucjf-a4-01      | AMD EPYC 9354                  | 2       | 32    | 2  | 384 GB |gpu:L40:1 mps:100 vram:48G| IB, NVIDIA Tesla L40 48 GB |
| r[31-35]        | Intel(R) Xeon(R) Gold 6140 CPU | 2       | 18    | 2  | 128 GB |      | IB(EDR)           |
| r[36-39,41-47]  | Intel(R) Xeon(R) Gold 6240 CPU | 2       | 18    | 2  | 128 GB |      | IB(EDR)           |
| r40             | Intel(R) Xeon(R) Gold 6240 CPU | 4       | 18    | 2  | 512 GB |      | IB(EDR)           |
| r50             | AMD EPYC 9654                  | 2       | 96    | 2  | 1.5 TB |gpu:MI210:2 vram:64G| IB(EDR), AMD MI210 |
| fuuk-1[2-4]     | AMD EPYC 9654                  | 1       | 96    | 2  | 1.5 TB |gpu:L40:4 vram:48G| IB(EDR), 4x NVIDIA TESLA L40s 48 GB |

Chimera node features

|Feature| description |
|--------|-|
| `ib`   | node has infiniband|
|`local_scratch`| nde has local scratch|
| `avx2`, `sse42`, `fma`,  | CPU support this instruction set| 
|`avx512_old`| CPUs supporting `avx512f`, `avx512dq`, `avx512cd`, `avx512bw`, `avx512vl`, `avx512_vnni` instruction sets (older Intel CPUs)|
| `avx512` | supports `avx512_old` and `avx512ifma`, `avx512_bf16`, `avx512vbmi`, `avx512_vbmi2`, `avx512_bitalg`,`avx512_vpopcntdq` instruction sets (new AMD / Intel CPUs)|
|`#C`| Procesor has # CPU cores|
|`#T`| Core has # SMT threads |
|`#N`| Machine has # numa nodes |
|`#S`| Machine has # sockets |



Chimera partitions

| Name             | Nodes                                | Priority  | Timelimit | DefaultTime | Preemption   | Oversubs | Allowed accounts               | Intended use |
| ---------------- | ------------------------------------ | --------- | --------- | ----------- | ------------ | -------- | ------------------------------ | ------------ |
| ffa              | all without debug and subset of auuk, ucjf, mff-a2-[01-02] | 5         | 12 hours  | 1 hour        | no           | FORCE:1  | any                            | free-for-all, any job |
| ffa-short        | all without debug and subset of auuk, ucjf, mff-a2-[01-02] | 5         | 2 hours   | 30 minutes    | no           | FORCE:1  | any                            | default, free-for-all, any job|
| ffa-long         | mff-a2-19                            | 5         | 07 days   | 1 day       | no           | FORCE:1  | any                            | free-for-all, any, QoS only 1/2 of the node |
| ffa-preempt*     | all without mff-a2-[01-02]                                 | 5         | 36 hours  | 1 hour        | REQUEUE      | NO       | any                            | free-for-all, any job, with preemption |
| ffa-checkpoint*  | all without subset of ucjf , mff-a2-[01-02]                | 5         | 36 hours  | 1 hour        | REQUEUE      | NO       | any                            | free-for-all, jobs with checkpointing, with preemption, grace time 60s |
| auuk             | subset of c[00-99]t[00-99], auuk-a2-[07-15,17-18]          | 40        | 8 hours   | 4 hours       | NO           | FORCE:1  | auuk                           | jobs for AUUK |
| auuk-preempt*    | c[00-99]t[00-99], auuk-a2-[07-15,17-18]                    | 20        | 12 hours  | 8 hours       | REQUEUE      | FORCE:1  | auuk                           | jobs for AUUK with preemption |
| utf              | utf-* node[0-1,6], r4[1-7]           | 40        | infinite  | 1 day       | no           | FORCE:1  | utf                            | jobs for UTF  |
| utf-pejcha       | r-[31-40]                            | 40        | infinite  | 1 day       | no           | FORCE:1  | utf-pejcha                     | jobs for Pejcha's group  |
| fuuk             | deepth1[3-4]                         | 40        | infinite  | 1 day       | no           | FORCE:1  | fuuk                           | jobs for FUUK |
| fuuk-barvik      | fuuk-1[2-4]                          | 40        | infinite  | 1 day       | no           | FORCE:1  | fuuk-barvik                    | jobs for Barvik's group |
| ucjf             | ucjf-*                               | 40        | infinite  | 1 day       | no           | FORCE:1  | ucjf                           | jobs for UCJF |
| math             | r[41-47,50]                          | 40        | infinite  | 1 day       | no           | FORCE:1  | math                           | jobs for mathematical section |
| debug            | onos[3-4]                            | 100       | 12 hour   | 1 hour      | no           | FORCE:1  | any                            | for interactive jobs only, QoS debug (max CPUs:24 and memory:200GB) |
| checkpoint       | onos[3-4], c14t3[2-3]                | 50        | 10 minutes| 5 minutes   | no           | FORCE:1  | any                            | for starting jobs with checkpoints |
| gpu-ffa          | mff-a2-19, ucjf-a4-01, fuuk-1[2-4]   | 10        | 12 hours  | 1 hour      | no           | FORCE:1  | any                            | only jobs that require gpu |
| edu              | mff-a2-[01-02]                       | 20        |  3 hours  | 2  hours    | no           | FORECE:1 | any                            | for education, only one job per user

> *preemptable partitions have lower fair share billing (TRESBillingWeights="CPU=") ; ffa-preempt 0.1; ffa-checkpoint 0.15; auuk-preempt 0.25



#### Special features

There are modules (LMod) accessible via command `module`:
`module avail` lists modules
`module load module_name` loads module to environment.

Preemption is set up on partitions indicated by preemption set to `REQUEUE` or `GANG,SUSPEND` (**ffa-preempt**, **ucjf**, **auuk**).
Job running in lower priority partitions with preemption set up can be preempted by jobs on partitions with
higher priority (with Oversubs set to FORCE:1) if there are no free resources left for the new job.

If **Preemption** is set to REQUEUE Slurm sends signal (SIGCONT unless specified differently in sbatch,...) to low priority job and
kills it one minute later (if it not ends itself). If job is killed this way it will be requeued.

If **Preemption** is set to GANG,SUSPEND
the low priority job is suspended and resumed after the high priority job ends.

Details about preemtpion can be found [here](https://slurm.schedmd.com/preempt.html).

#### MPI

###### Intel MPI
For now there is one tested MPI implementation installed - Intel oneAPI. To use it you need to load module oneapi/mpi like this:
`module load oneapi/mpi`
This is needed for compiling code and for running MPI job (in script submitted by sbatch before running mpi application).
```diff
- You need to run MPI jobs with `srun` and not with `mpirun` or `mpiexec` provided by Intel.
```
`mpirun` and `mpiexec` are not binding CPUs to tasks correctly.

This provides wrapers for fortran na C compilers `mpiifort`, `mpiicc`, `mpif90`, etc. For details please see Intel [documentaton](https://www.intel.com/content/www/us/en/develop/documentation/mpi-developer-guide-linux/top.html).

###### OpenMPI
There are compiled versions of OpenMPI.
 - If unsure you should use module `openmpi-5.0.8`.
   - The default `openmpi-5.0.8` (or `openmpi-5.0.8/ib`) is set up for running on nodes with InfiniBand - this can (and should) be enforced by `-C ib` in SLURM options. 
   - If you want it to run program on node(s) without InfiniBand, please use `openmpi-5.0.8/tcp`. (If you are running over multiple nodes, and any of them are without Infiniband, you have to use this version.)
   - These modules will correctly setup `srun` options and environment variables for runtime. It should not matter which one you are using at compile time, unless the building of software is running MPI code. But the result should be independent on `ib` / `tcp`.
   - There is also package `openmpi-5.0.8-base` which sets up only paths to OpenMPI environment.

- You can use older module `openmpi-4.1.5`. 
    - It is set up similarly: default for InfiniBand nodes (`openmpi-4.1.5/ib`), nodes without InfiniBand `openmpi-4.1.5/tcp` and  `openmpi-4.1.5-base`.
    

#### Fair Share
On Chimera cluster is enabled Fair Share feature.
This means that users that used cluster less in last weeks have higher priority than users that used cluster more. For now, this only account for allocated CPUs.  The half-life of the fairshare is now set to 14 days - after 14 day the job contribution to "Fair share" is halved.
Overall priority calculation also depends on other factors - size of job, how long it is in queue and more. Details can be found [here](https://slurm.schedmd.com/priority_multifactor.html). For detailed parameter on cluster priority settings please check `scontrol show config | grep Priority`.

### Default resources
Default resources for all partitions are now set to:
- 1 node
- 1 cpu
- 1 GB RAM/CPU

Time is set between 30 minutes and one day depending on the partition, please check `scontrol -o show part | cut -d' ' -f1,8` for details on front-end of cluster you want to run jobs at.

### Checkpointing

- The `DMTCP` system is available on the cluster. An example of use is in the directory `user_checkpointing_demo/low_priority_job_with_dmtcp/` of this repository.
- it is done with `source scheckpoint` - there are some options available - see `source scheckpoint -h` for details

### Mathematica
- only for MFF and CERGE
- is available on chimera cluster by using `module load Mathematica`
- for interactive GUI environment you should:
	- allocate resources you need - ie. `salloc -n1`; you can specify other things, like partition - the result sould look something like:

			$salloc -n1 -p ffa-preempt
			salloc: Granted job allocation 9324110
        	salloc: Nodes ucjf-asusb2 are ready for job
			$

	- leave this open (in this case we got node `ucjf-asusb2`)
	- login to cluster with X11 forwarding `ssh -Y -l username hpc.troja.mff.cuni.cz`
	- login to allocated node with X11 forwarding `ssh -Y NODE` (in our example the `NODE` is `ucjf-asusb2`)
	- run mathematica - `module load Mathematica; mathematica`


## Charliecloud (gpulab nodes)

Charliecloud provides user-defined software stacks (UDSS) for HPC.
It allows you to run nearly any software stack (like TensorFlow) on the cluster even it is not system-wide installed and available.
All informations about Charliecloud can be found on its [Charliecloud documentation](https://hpc.github.io/charliecloud/) page.


### Basic workflow

This workflow is valid for Charliecloud version 0.34.

1. #### Get or create Docker image

We don't use Docker anymore, as it poses significant security hole. Docker was replaced by [Podman](https://podman.io/).
Podman is installed on all worker nodes and the command `podman` is aliased by `docker`.
You can use scripts using `docker` command without changing them.

You can either pull already prepared Docker image (e.g. for TensorFlow) or you may create your own one from `Dockerfile`.
You will make this step only once for the given UDSS.
Of course, you must restart the whole workflow, if there is a new version of the UDSS.

If you are pulling a Docker image, use

`ch-image pull dockertag`

If you are building your own Docker image from `Dockerfile`, use

`ch-image build -t dockertag .`

2. #### Create tar/directory image from Docker image

You must convert prepared Docker image to a directory structure.
You will make this step only once for the given UDSS.

Use command

`ch-convert -i ch-image -o dir dockertag imgdir`

which exports Docker image with *dockertag* tag to a directory *imgdir*,
which will be used for executing the container.

3. #### Import CUDA libraries

This step is required only for UDSS with CUDA requirement (like TensorFlow with GPU support).
If your UDSS does not require CUDA, skip this step.
You will make this step only once for the given UDSS, but you should repeat this step, when new host drivers/CUDA are installed.
It works only with directory structure.
All commands for this step must be run on volta[01-05] or ampere[01-02] nodes.

Execute on gpulab

`srun -p gpu-short --gpus=1 ch-fromhost --nvidia imgdir`

which copies some necessary CUDA files from the host to your image directory structure.

4. #### Execute created UDSS by SLURM

This step is executed many times as necessary on any node of parlab and gpulab clusters.

For an interactive job run

`srun <slurm params> ch-run <charlie options> imgdir <my img command>`

You will probably use SLURM batch job mode more often, as the length of computation is usually several hours or days.
In this case use

`ch-run <charlie options> imgdir <my img command>`

in your shell script passed to the `sbatch` SLURM command.

Beware, that by default, your home is not bound to the image.
You may bind home or any additional directories by using options `--bind=/some/dir`
(which will appear as `/some/dir` in your UDSS environment) or by `--bind=/source/dir:/dest/dir`.
All Charliecloud images have free mountpoints `/mnt/[0-9]`.

### Advanced techniques, troubleshooting, and notes

#### Builders

Charliecloud brings a term **builder**. A builder is anything capable of producing a Linux filesystem tree
from either a prepared set of container images or from some container description.

Currently Charliecloud supports several builders. We have enabled only two of them: Docker and ch-image.

`ch-image` command/builder builds an image in unprivileged mode from `Dockerfile`.
This command doesn't need Docker, it can be run anywhere on workers.

When using `ch-image` command, you have to use parameter `-i ch-image` in step 2.

#### CUDA import errors

It can happen, executing `ch-fromhost` from the 3rd step will produce some errors, which looks something like

	/sbin/ldconfig.real: Can't stat /usr/local/nvidia/lib: No such file or directory
	/sbin/ldconfig.real: Can't stat /usr/local/nvidia/lib64: No such file or directory
	/sbin/ldconfig.real: Can't stat /usr/local/lib/x86_64-linux-gnu: No such file or directory

Ignore safely these errors, they do no harm to you.

This warning/notice was contributed by Vít Kabele.

#### CUDA is not working inside your UDSS

If you suspect that CUDA is not working inside your container, run the `nvidia-smi` command from the container command line.
If `nvidia-smi` prints the CUDA version correctly, then CUDA is functional. However, if it prints "ERR", CUDA does not work.

In this case, follow the checklist:
- Did you correctly import the CUDA libraries? See step 3 of the basic workflow.
- Is `libcuda.so` (or `libcuda.so.1`) loadable? Check `LD_LIBRARY_PATH` environment variable inside your container.
If not set to the CUDA library directory, set it to the correct path, e.g. `export LD_LIBRARY_PATH=/usr/local/cuda/lib64`.
Be careful, if the variable is already set to some additional paths.

## Apptainer

Apptainer (formerly *Singularity*) simplifies the creation and execution of containers, ensuring software components are encapsulated for portability and reproducibility.
Please see [Apptainer documentation](https://apptainer.org/docs/user/1.2/) for details.

If you want to run containers with `--fakeroot` option on the cluster, please contact the administrator to be added to `/etc/subuid` and `/etc/subgid` on headnode.


## Choosing an interface

The [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) allows AI clients such as Codex, Gemini, Claude Code, and OpenCode to use tools provided by external processes. Chimera currently has two MCP-based AI-on-cluster projects, but they serve different purposes:

| Project | Purpose | Execution model |
| ------- | ------- | --------------- |
| [ChimeraMCP](https://gitlab.mff.cuni.cz/rse/projects/ChimeraMCP) | Guarded access to Slurm, jobs, logs, user-scoped files, and archive operations | The local AI client starts typed MCP servers over SSH; submitted compute runs through Slurm |
| [ChimeraSandbox](https://gitlab.mff.cuni.cz/rse/projects/apptainermcp) | A persistent Apptainer environment with shell and file tools for CPU/GPU workloads | A Slurm job holds an Apptainer instance on a compute node and the AI client attaches to it |

### ChimeraMCP

ChimeraMCP lets MCP-capable AI agents operate on Chimera through guarded Slurm and archive tools. Instead of giving a model raw shell access, you connect it to typed MCP servers that run on Chimera under your Unix account.

The companion [ChimeraMCP OpenCode HPC Agent](https://gitlab.mff.cuni.cz/rse/projects/chimeramcp-opencodeagent) repository provides ready-made OpenCode agent, subagent, and configuration templates for using ChimeraMCP from OpenCode.

ChimeraMCP currently ships two MCP servers through the administrator-managed `chimera-mcp` module:

- `chimera-slurm`, started by `chimera-slurm-mcp`, for Slurm queue inspection, job submission, job logs, accounting, quotas, software module lookup, and user-scoped workspace file operations.
- `chimera-filecompress`, started by `chimera-filecompress-mcp`, for archive planning, creation, and extraction for `zip`, `tar`, `tar.gz`, `tar.bz2`, and `tar.xz`.

#### Warning

An MCP-connected model can act on your behalf through the tools you expose to it. With ChimeraMCP, that can include reading and writing files in the configured work area, submitting jobs, canceling or updating your jobs when those tools are enabled, moving files, creating or extracting archives, and reading job logs.

Anything returned by the server can also be shown to the connected AI service, including file paths, log contents, job metadata, error messages, and snippets of your data. Do not connect ChimeraMCP to a model or provider you do not trust, and do not ask the agent to handle secrets, passwords, tokens, private keys, or sensitive data you are not allowed to share.

ChimeraMCP adds guardrails, but it is not a security sandbox. You remain responsible for your files, jobs, resource usage, and compliance with Chimera site policy.

#### Prerequisites

Normal users do not install ChimeraMCP locally and do not need to clone the ChimeraMCP repository. The standard setup is the shared Chimera module:

```bash
module load chimera-mcp
```

You need:

- A Chimera account and SSH access, usually `USER@hpc.troja.mff.cuni.cz`.
- Access to the shared `chimera-mcp` module on Chimera.
- A writable Chimera work area. In the standard configuration, the user-facing data root is `/work/$USER`, derived from `/work/MCP_scripts` and `/work/MCP_logs`.
- A local MCP-capable client, such as Codex, Gemini, Claude Code, or OpenCode.

Your local client starts each MCP server over SSH. Replace `USER@hpc.troja.mff.cuni.cz` with your Chimera login. Configure key-based SSH authentication before using `BatchMode=yes`, because MCP clients cannot answer an interactive password prompt from an MCP subprocess.

#### Check the installation

Log in to Chimera and confirm that the module exposes both commands:

```bash
ssh USER@hpc.troja.mff.cuni.cz
module load chimera-mcp
which chimera-slurm-mcp
which chimera-filecompress-mcp
```

#### Connect your model

Configure both servers as separate stdio MCP servers. Keep the remote `bash -lc 'module load chimera-mcp && exec ...'` command as one SSH argument.

<details>
<summary><strong>Codex</strong>: <code>~/.codex/config.toml</code></summary>

```toml
[mcp_servers.chimera-slurm]
command = "ssh"
args = [
  "-T",
  "-o", "BatchMode=yes",
  "-o", "LogLevel=ERROR",
  "USER@hpc.troja.mff.cuni.cz",
  "bash -lc 'module load chimera-mcp && exec chimera-slurm-mcp'"
]

[mcp_servers.chimera-filecompress]
command = "ssh"
args = [
  "-T",
  "-o", "BatchMode=yes",
  "-o", "LogLevel=ERROR",
  "USER@hpc.troja.mff.cuni.cz",
  "bash -lc 'module load chimera-mcp && exec chimera-filecompress-mcp'"
]
```

</details>

<details>
<summary><strong>Gemini</strong>: <code>~/.gemini/settings.json</code></summary>

```json
{
  "mcpServers": {
    "chimera-slurm": {
      "command": "ssh",
      "args": [
        "-T",
        "-o", "BatchMode=yes",
        "-o", "LogLevel=ERROR",
        "USER@hpc.troja.mff.cuni.cz",
        "bash -lc 'module load chimera-mcp && exec chimera-slurm-mcp'"
      ],
      "timeout": 600000,
      "trust": true
    },
    "chimera-filecompress": {
      "command": "ssh",
      "args": [
        "-T",
        "-o", "BatchMode=yes",
        "-o", "LogLevel=ERROR",
        "USER@hpc.troja.mff.cuni.cz",
        "bash -lc 'module load chimera-mcp && exec chimera-filecompress-mcp'"
      ],
      "timeout": 600000,
      "trust": true
    }
  }
}
```

</details>

<details>
<summary><strong>Claude Code</strong>: <code>claude mcp add</code></summary>

```bash
claude mcp add --scope user --transport stdio chimera-slurm -- \
  ssh -T -o BatchMode=yes -o LogLevel=ERROR USER@hpc.troja.mff.cuni.cz \
  "bash -lc 'module load chimera-mcp && exec chimera-slurm-mcp'"

claude mcp add --scope user --transport stdio chimera-filecompress -- \
  ssh -T -o BatchMode=yes -o LogLevel=ERROR USER@hpc.troja.mff.cuni.cz \
  "bash -lc 'module load chimera-mcp && exec chimera-filecompress-mcp'"
```

Check registration with:

```bash
claude mcp list
```

</details>

<details>
<summary><strong>OpenCode</strong>: <code>~/.config/opencode/opencode.jsonc</code></summary>

Merge the following `mcp` entries into `~/.config/opencode/opencode.jsonc` or a project `opencode.jsonc` and replace `USER`:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "chimera-slurm": {
      "type": "local",
      "command": [
        "ssh",
        "-T",
        "-o", "BatchMode=yes",
        "-o", "LogLevel=ERROR",
        "USER@hpc.troja.mff.cuni.cz",
        "bash -lc 'module load chimera-mcp && exec chimera-slurm-mcp'"
      ],
      "enabled": true,
      "timeout": 600000
    },
    "chimera-filecompress": {
      "type": "local",
      "command": [
        "ssh",
        "-T",
        "-o", "BatchMode=yes",
        "-o", "LogLevel=ERROR",
        "USER@hpc.troja.mff.cuni.cz",
        "bash -lc 'module load chimera-mcp && exec chimera-filecompress-mcp'"
      ],
      "enabled": true,
      "timeout": 600000
    }
  }
}
```

Check registration with:

```bash
opencode mcp list
```

OpenCode users can also use the ready-made agent templates in the [ChimeraMCP OpenCode HPC Agent](https://gitlab.mff.cuni.cz/rse/projects/chimeramcp-opencodeagent) project.

</details>

Restart your AI client after changing MCP configuration. If the remote shell cannot find `module`, your site may need to source its Lmod initialization script before `module load chimera-mcp`. See the [ChimeraMCP documentation](https://gitlab.mff.cuni.cz/rse/projects/ChimeraMCP/-/tree/main/docs) for the deployment and troubleshooting reference.

#### Optional agent skills

ChimeraMCP also provides optional [agent skill bundles](AI_on_cluster/ChimeraMCP/skills/README.md). Skills are reusable instruction packages for AI clients; they are not MCP servers, they do not install ChimeraMCP, and they do not replace the MCP server configuration above. Use them when your client supports skills and you want the model to load Chimera-specific workflow guidance on demand.

Included skills:

- `chimera-slurm-mcp`: live and managed Slurm workflows through Chimera Slurm MCP, including queue inspection, owned-job logs, accounting, failure diagnosis, workspace files, and guarded submissions.
- `chimera-filecompress-mcp`: archive planning, creation, extraction, protected-root handling, and Slurm routing for large archive work through Chimera FileCompress MCP.
- `chimera-cluster`: general Chimera guidance outside MCP-managed workflows, such as SSH, JupyterHub, storage, modules, manual Slurm, GPUs, and Apptainer.
- `chimera-mcp-server-authoring`: developer guidance for adding or reshaping ChimeraMCP servers so they follow the repository conventions.

Copy the whole skill directory, not only `SKILL.md`, because bundled references and assets are part of the skill behavior. From a ChimeraMCP checkout:

```bash
mkdir -p ~/.codex/skills ~/.gemini/skills ~/.claude/skills
cp -R skills/chimera-cluster skills/chimera-mcp-server-authoring skills/chimera-slurm-mcp skills/chimera-filecompress-mcp ~/.codex/skills/
cp -R skills/chimera-cluster skills/chimera-mcp-server-authoring skills/chimera-slurm-mcp skills/chimera-filecompress-mcp ~/.gemini/skills/
cp -R skills/chimera-cluster skills/chimera-mcp-server-authoring skills/chimera-slurm-mcp skills/chimera-filecompress-mcp ~/.claude/skills/
```

Restart the AI client after installing or updating skills so it reloads the skill metadata.

#### What ChimeraMCP can do

All ChimeraMCP tools return the same response envelope:

```json
{"ok": true, "data": "..."}
```

or:

```json
{"ok": false, "error_code": "...", "message": "..."}
```

Use `chimera-slurm` to inspect current cluster state and software modules:

- current partitions and node state
- cached Lmod module search and module details
- user Slurm associations and disk quotas

Use `chimera-slurm` to inspect your jobs and the live cluster queue:

- `slurm_list_jobs` lists jobs owned by the current Unix user.
- `slurm_list_all_jobs` gives a read-only live view of the whole queue.
- array tools summarize and inspect array job tasks.

Owned-job detail tools can read job metadata, accounting, efficiency, stdout and stderr logs, and structured failure diagnoses.

Use `chimera-slurm` to submit work through guarded paths:

- built-in templates with `slurm_list_templates` and `slurm_submit_template_job`
- stored custom scripts
- project scripts under the configured user data root
- managed DMTCP checkpoint/restart for custom and project scripts
- advanced user-managed DMTCP tools for scripts that intentionally contain DMTCP logic

The active security level may hide custom, project, DMTCP, or control tools. Template submission does not support DMTCP in v1.

When enabled by policy, `chimera-slurm` can cancel or update jobs owned by the server user. It does not allow controlling other users' jobs.

The Slurm server can also operate on files and directories under the configured user data root, normally `/work/$USER` in the standard deployment. These tools are intended for job inputs, small generated notes, scripts, logs, and workflow artifacts. They do not provide unrestricted filesystem access.

Use `chimera-filecompress` for archive work:

- plan archive creation or extraction before doing it
- create `zip`, `tar`, `tar.gz`, `tar.bz2`, and `tar.xz` archives
- safely extract supported archive formats
- reject plain relative paths, unsafe archive members, and protected-root access
- auto-rename outputs rather than overwrite existing paths
- route large archive operations through Slurm using the `FileCompress_job` template

FileCompress accepts absolute paths and `~/` paths. Filesystem permissions still apply.

#### Safe use

Start with read-only requests when connecting a new model. Before submitting, canceling, updating, extracting, or moving data, ask the agent to show the exact planned action, target paths, job inputs, Slurm overrides, and expected outputs.

Good first prompts:

```text
Use Chimera Slurm MCP to list my running and pending jobs.
```

```text
Use Chimera Slurm MCP to summarize current partitions, quotas, and my Slurm associations.
```

```text
Use Chimera FileCompress MCP to plan creating ~/project.tar.gz from ~/project.
```

The complete server reference, configuration rules, and security model are maintained in the [ChimeraMCP documentation](https://gitlab.mff.cuni.cz/rse/projects/ChimeraMCP/-/tree/main/docs).

### ChimeraSandbox (experimental)

ChimeraSandbox runs a small MCP server inside a persistent Apptainer instance on a Slurm compute node. It gives an OpenCode agent four tools:

| Tool | Purpose |
| ---- | ------- |
| `bash` | Run shell commands inside the container |
| `read_file` | Read a size-limited file |
| `write_file` | Create or overwrite a file |
| `list_dir` | List a directory |

The submitted project directory is mounted at `/workspace`, which is the intended location for inputs, environments, checkpoints, logs, and results. Container state and background processes survive OpenCode reconnects while the Slurm job remains alive. Only data stored on a persistent host filesystem, including the mounted `/workspace`, survives the end of the allocation.

This workflow is experimental. Unlike ChimeraMCP, its `bash` and file tools do not enforce guarded Slurm, path, or owned-job policies. Commands run with the permissions of your cluster account, and `/workspace` contains real user data. Review requested operations accordingly.

#### Add ChimeraSandbox to OpenCode

Start the ChimeraSandbox Slurm job and Apptainer instance yourself, using the instructions in the [ChimeraSandbox repository](https://gitlab.mff.cuni.cz/rse/projects/apptainermcp). The GitLab project path is still `rse/projects/apptainermcp.git`; the user-facing project name is ChimeraSandbox.

Build the image on a compute node, not on the login node:

```bash
git clone git@gitlab.mff.cuni.cz:rse/projects/apptainermcp.git ChimeraSandbox
cd ChimeraSandbox
./scripts/build_submit.sh --output ~/chimera_sandbox.sif
```

Then start a sandbox allocation. On the current Chimera configuration, GPU jobs use the `gpu-ffa` partition; use `sinfo` if you need to verify partition names live:

```bash
mkdir -p ~/chimera-sandbox-test
./scripts/submit.sh \
  --mount ~/chimera-sandbox-test \
  --image ~/chimera_sandbox.sif \
  --partition gpu-ffa \
  --gres mps:5 \
  --time 01:00:00 \
  --user USER \
  --login hpc.troja.mff.cuni.cz
```

When the job is running, note its Slurm job ID. You can get it from the submission output or with:

```bash
squeue --me
```

Then add a `chimera-sandbox` entry to the `mcp` object in your global or project OpenCode configuration. Replace `USER` with your Chimera login and `JOBID` with the running Slurm job ID:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "chimera-sandbox": {
      "type": "local",
      "enabled": true,
      "command": [
        "ssh", "-T", "USER@hpc.troja.mff.cuni.cz",
        "srun", "--jobid=JOBID", "--overlap",
        "apptainer", "exec", "instance://chimerasandbox",
        "python", "-u", "/opt/mcp/server.py"
      ]
    }
  }
}
```

Restart OpenCode and check the connection with `opencode mcp list`. A useful first request is to run `hostname`, `nvidia-smi`, `pwd`, and `python --version` through the sandbox `bash` tool. We verified this flow on Chimera with a running `chimera-sandbox` MCP server, `/workspace` bind mount, Python/R/mamba/git, and an NVIDIA L40 visible inside the container.

Closing or restarting OpenCode does **not** release the Slurm allocation. The requested CPUs, memory, and GPU remain reserved until the job reaches its wall-time or is canceled. When finished, remove or disable the stale OpenCode entry and release the resources:

```bash
scancel JOBID
```


## Installed software on clusters

### KSI/DSE clusters

- **Charliecloud** - version 0.38


### Chimera

- **apptainer** - version 1.2.2-1.el9

There are also some preinstalled software packages on Chimera cluster which are accessible via the `module` command.
For now there is:
- **Molpro** - version 2012.1. You need a licence token to run this software. Please contact one of K.Houfek, Z.Masin, J.Eliasek to obtain one. Once you have it the license token should be placed into `~/.molpro/token`. Alternatively, run Molpro like this `molpro molpro.inp -k "CONTENT_OF_THE_LICENSE_TOKEN" \`.
- **Intel Parallel Studio XE 2019**
- **Intel oneAPI**
- **UKRmol+ codes** - this code comprises of UKRmol-in (v 3.1.2) and UKRmol-out (v 3.1). For details contact one of Z.Masin, K.Houfek.
- **DMTCP** - Version 3.0. It is used with the checkpoint scripts.
- **Mathematica** - version 13.3


## Accounts management

Accounts are associated with resouces / project and can be used to access specific partitions. If you want to access these partitions and you don't have access to the account you should, please contact coordinator for the account or administator.

coordinators for chimera (physics section):

| account | coordinator | 
|---------|-------------|
|auuk|doc. RNDr. Ladislav Šubr, Ph.D.|
|utf|Jiří Eliášek|

coordinators for chimera (mathematics section):

| account | coordinator | 
|---------|-------------|
|math|Jaroslav Hron|



### Useful `sacctmgr` commands:
list all accounts with coordinators for each account:
>`sacctmgr show accounts WithCoord`

list all users and their associations:
>`sacctmgr show association`

to restrict the output to specific user:
>`sacctmgr show association user=user_name`


### For coordinators 
associated with account (ACCOUNT in examples) and sub accounts:

#### Can do:
- add users to account they are administrating,
>`sacctmgr add user user_name account=ACCOUNT`

- remove users from accounts
>`sacctmgr  delete user user_name account=ACCOUNT`

- create sub-accounts
>`sacctmgr add account Name=some-sub-account Description="description of subaccount" Parent=ACCOUNT`

#### Cannot do:
- change default account:
>`sacctmgr modify user where user=user_name set defaultaccount=ACCOUNT`
- delete subaccounts:
>`sacctmgr delete account Name=some-sub-account`

please contact administrator if you need to do these operations.
