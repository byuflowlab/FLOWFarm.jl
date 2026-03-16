# How-to Guide

## Multi-threading
Multi-threading is available for the calculation of annual energy production (AEP). It can be
enabled as follows in a bash terminal in Linux/OS prior to launching a julia session:

```
export JULIA_NUM_THREADS=<number of threads>
```
For enabling multi-threading on other shells/systems please see the julia parallel-computing
docs here: https://docs.julialang.org/en/v1/manual/parallel-computing/.

**References**
For more information on using julia in a distributed environment, please see https://docs.julialang.org/en/v1/manual/parallel-computing/.
