[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://byuflowlab.github.io/FLOWFarm.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/FLOWFarm.jl/dev)
[![Build Status](https://github.com/byuflowlab/FLOWFarm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/FLOWFarm.jl/actions/workflows/CI.yml?query=branch%3Amain)

# FLOWFarm.jl

**Summary:** Wind farm simulation tool for gradient-based optimization.

**Authors:** Jared J. Thomas, Andrew P.J. Stanley

## Features
- Compatible with ForwardDiff for gradiant-based optimization
- Swap out models without changing anything else in the simulation setup
- Smooth/continous model implementations
- Runs on a single core, across multiple cores (threaded), or on multiple machines (distributed).
- Designed so that new model implementations can be included by adding a single method
- Allows for Wake Expansion Continuation (WEC) as described [here](http://flowlab.groups.et.byu.net/preprints/Thomas2021.pdf)

## Installation

### Install FLOWFarm

```
julia
(v1.x) pkg> add FLOWFarm
```

## Documentation

* Begin with the [quick start tutorial](https://flow.byu.edu/FLOWFarm.jl/Tutorial/).
* An explanation of how to utilize sparsity can be found in [wind farm struct](https://flow.byu.edu/FLOWFarm.jl/Wind_Farm_Struct/).
* More advanced topics are covered in the [how-to guide](https://flow.byu.edu/FLOWFarm.jl/How_to/).
* Theory details, and links, can be found in the [theory](https://flow.byu.edu/FLOWFarm.jl/Explanation) page.
* Doc strings can be found in the [references](https://flow.byu.edu/FLOWFarm.jl/Reference/) page.

**Citing:**
Thomas, McOmber, and Ning "Wake Expansion Continuation: Multi-Modality Reduction in the Wind Farm Layout Optimization Problem" *Wind Energy* (in review), -->

