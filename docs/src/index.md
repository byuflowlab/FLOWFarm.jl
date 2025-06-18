# FLOWFarm.jl

**Summary:** Wind farm simulation tool for gradient-based optimization.

**Authors:** Jared J. Thomas, Andrew P.J. Stanley

## Features
- Swap out models without changing anything else in the simulation setup.
- Smooth/continous model implementations.
- Runs on a single core, across multiple cores (threaded), or on multiple machines (distributed).
- Designed so that new model implementations can be included by adding a single method.
- Allows for Wake Expansion Continuation (WEC) as described [here](Explanation.md).

## Installation

### Install FLOWFarm

```
julia
(v1.x) pkg> add FLOWFarm
```

* Begin with the [quick start tutorial](Tutorial.md).
* An explanation of how to utilize sparsity can be found in [wind farm struct](WindFarmStruct.md).
* More advanced topics are covered in the [how-to guide](How_to.md).
* Theory, details, and sources can be found in the [theory](Explanation.md) page.
* Doc strings can be found in the [references](Reference.md) page.

**Citing:**
Thomas, McOmber, and Ning "Wake Expansion Continuation: Multi-Modality Reduction in the Wind Farm Layout Optimization Problem" *Wind Energy* (in review), -->

