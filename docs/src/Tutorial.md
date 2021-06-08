# Tutorial

This section is designed to familiarize the user with the basic function of Flow Farm. For help with specifics problems or obtaining desired outputs refer to the [How-to guide](How_to.md).

Before attempting to use Flow Farm add the Flow Farm package as well as all of its dependencies. Refer to the [introduction](index.md) for installation procedure.

For basic functionality Flow Farm requires two inputs, the model set, and the problem set. The model set includes all of the different parameters which flow farm will use to analyze your wind farm. The problem set includes your farm definition. Theory behind each one of these selections can be found in the [theory](Explanation.md) section.

## Model Set
Model set will recquire a Wake Deficit Model, Wake Deflection Model, Wake Combination Model, and a Local Turbulence Intensity Model
* Deficit Model options: JensenTopHat, JensenCosine, MultiZone, GaussOriginal, GaussYaw, GaussYawVariableSpread, GaussSimple
* Deflection Model options: GaussYawDeflection, GaussYawVariableSpreadDeflection, JiminezYawDeflection, MultizoneDeflection
* Combination Model options: LinearFreestreamSuperposition, SumOfSquaresFreestreamSuperposition, SumOfSquaresLocalVelocitySuperposition, LinearLocalVelocitySuperposition
* Turbulence Model options: LocalTIModelNoLocalTI, LocalTIModelMaxTI
## Problem set 
Problem set will recquire number of turbines, turbine selection, ambient wind speeds, wind directions, and wind probabilities. Unlike the model set most of these inputs are user defined. Turbine selection is more limited. Flow Farm turbine choices include NREL 5 MW, Vestas v80, and Simens SWP 3.6. Refer to the [How-to guide](How_to.md) for help inputing files for other turbines.