![Tests](https://github.com/byuflowlab/FLOWFarm.jl/actions/workflows/test.yml/badge.svg)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://byuflowlab.github.io/FLOWFarm.jl/)


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
(v1.x) pkg> dev https://github.com/byuflowlab/FLOWFarm.jl.git
```

### Enable NaN Safe Mode in ForwardDiff
NaN Safe Mode must be enables in ForwardDiff for ForwardDiff to work properly with FLOWFarm.

```
julia
(v1.x) pkg> dev ForwardDiff
```
```
$ cd ~/.julia/dev/ForwardDiff/src/
```
In `prelude.jl`, on the first line, set `const NANSAFE_MODE_ENABLED = true` and save the file. 
For more information see the ForwardDiff documentation at 
http://www.juliadiff.org/ForwardDiff.jl/latest/user/advanced.html

## Testing

To test FLOWFarm, run the following from the top directory:

```
julia
julia
]
activate .
test
```

## Documentation

* Begin with the [quick start tutorial](https://flow.byu.edu/FLOWFarm.jl/Tutorial/).
* More advanced topics are covered in the [how-to guide](https://flow.byu.edu/FLOWFarm.jl/How_to/).
* Theory details, and links, can be found in the [theory](https://flow.byu.edu/FLOWFarm.jl/Explanation) page.
* Doc strings can be found in the [references](https://flow.byu.edu/FLOWFarm.jl/Reference/) page.

**Citing:**
Thomas, McOmber, and Ning "Wake Expansion Continuation: Multi-Modality Reduction in the Wind Farm Layout Optimization Problem" *Wind Energy* (in review), -->

