# Exercise 1


## Study how to submit jobs in SGE, how to check their state and how to cancel them.

### Submit jobs: 
1. Write SGE job script
2. Submit job ```qsub name_of_script```

### Check the state
```qstat```

### Cancel
```qdel job_id_list```


## Prepare a submission script that starts an arbitrary executable, e.g. /bin/hostname
```
#!/bin/bash

# Execute job in the queue "std.q" unless you have special requirements.
#$ -q std.q

# The batch system should use the current directory as working directory.
#$ -cwd

# Name your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N my_test_job

# Redirect output stream to this file.
#$ -o output.dat

# Join the error stream to the output stream.
#$ -j yes

#$ -pe openmpi-2perhost 8

/bin/hostname
```

## In your opionion, what are the 5 most important parameters available when submitting a job and why? 
What are possible settings of these parameters, and what effect do they have?

1. ```-q queuename```: Selects the queue where the job will be inserted. (Defines things like, max. execution hours, number of cpu slot, memory requieremnts, ...)
* ```std.q```: This is the default. General purpose queue. Default/maximum runtime: 240 hours.
* ```short.q```: or small test jobs. Limited number of CPU slots. Default/maximum runtime: 10 hours.
* ```bigmem.q```: Leo3e only. Jobs with high main memory requirements. Will run on the nodes equipped with 512GB of memory. Default/maximum runtime: 240 hours.

**Why**: Basic configurations needed for the SGE to select right queue necessery for the job.

2. ```-cwd```: execute job in current working directory. If you omit this option, you job will execute in $HOME, which is usually a bad idea. Input/output file names are relative to this directory.

**Why**: Basic configurations needed to define correct working directory.

3. ```-t 1-n```: Trivial parallelisation using a job array. Start n independent instances of your job (e.g. for extensive parameter studies). When the job is run, you use the environment variable $SGE_TASK_ID, which is set to a unique integer value from 1 .. n, to distinguish between the individual job instances (e.g. to initialize a random number generator, select an input file or compute parameter values).

**Why**: Simple way to call a programm with different settings e.g. see how parameters will affect the result or performance.


4. ```-pe parallel-environment number-of-slots```: If you run parallelized programs (MPI or shared memory), you need to specify a parallel environment and the number of processes/threads (= SGE slots) on which your parallel (MPI/OpenMP) application should run. By selecting a parallel environment you can also control how jobs are distributed across nodes.
parallel-environment:
* openmpi-1perhost
* openmpi-2perhost
* openmpi-4perhost
* openmpi-8perhost

number-of-slots: Number of hosts

**Why**: Basic configuration needed to define how much computational power we need or want to work with.

5. ```qsh -now n [...]```: The submission of interactive jobs is useful in situations where a job requires some sort of direct intervention. This is usually the case for X-Windows applications or in situations in which further processing depends on your interpretation of immediate results. A typical example for both of these cases is a graphical debugging session.

**Why**: Sometimes its necessary e.g. to debug.



## How do you run your program in parallel? What environment setup is required?
1. Write execution script.
2. Load within the script openmpi with ```module load openmpi/4.0.1``` (sets up the openmpi environment)
3. Execute programm with ```mpiexec -n 8 COMMAND```