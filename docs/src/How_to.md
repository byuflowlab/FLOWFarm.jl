# How-to Guide

## Multi-threading
Multi-threading is available for the calculation of annual energy production (AEP). It can be
enabled as follows in a bash terminal in Linux/OS prior to launching a julia session:

```
export JULIA_NUM_THREADS=<number of threads>
```
For enabling multi-threading on other shells/systems please see the julia parallel-computing
docs here: https://docs.julialang.org/en/v1/manual/parallel-computing/.

## Distributed Processing
Distributed parallel processing is available for the calculation of annual energy production (AEP). 

You may have to add `using Distributed` to your julia script and use the `@everywhere` macro 
in front of any functions you define that all processors will need access to. For an example, 
see `example_opt_6_38turb_round_distributed.jl`.

**Using Distributed Processing without an HPC Cluster Manager (e.g. on your local system)**

Distributed parallel processing can be enabled as follows when launching a julia session:

```
julia -p <number of processors>
```

**Using Distributed Processing with an HPC Cluster Manager (e.g. SLURM)**

The `-p` option to the julia call is unnecessary when running with a cluster manager. 
To work with cluster managers, add the following to your julia script (this example is for 
SLURM, but other managers are available as well):

```
using Distributed
using ClusterManagers

addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
@everywhere import FLOWFarm; const ff = FLOWFarm
```

Also include the `@everywhere` macro in front of any function definitions or include statements
in your julia script that all processors will need access to.

Your SLURM job script should look something like this:

```
#!/bin/bash -l
#SBATCH --ntasks=100
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH --time=01:00:00 # time=HH:MM:SS
#SBATCH -J "Your job name here"   # job name

module load julia

julia julia_script.jl
```

**References**
For more information on using julia in a distributed environment, please see https://docs.julialang.org/en/v1/manual/parallel-computing/.
