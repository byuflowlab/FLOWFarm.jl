# Tutorial

This tutorial covers the basics of FlowFARM. For more specifics refer to the [How-to guide](How_to.md).

There are three main steps to setting up and running an analysis in FLOWFarm:
- (1) setting up the problem description
- (2) setting up the analysis model set
- (3) running the analysis

Details for setting up an optimization will depend heavily on the
optimization package you are using, your objective, and your design variables. Optimization
examples using various packages are provided in the example scripts located in the test 
directory.

### (1) Setting up the problem description

```julia
import FLOWFarm; const ff = FLOWFarm

# define the rotor diameter
diameter = 80.0

# set initial turbine x and y locations
turbine_x = [-3.0, 0.0, 3.0, 0.0, 0.0, -1.5, 0.0, 1.5, 0.0].*diameter
turbine_y = [0.0, 3.0, 0.0, -3.0, 0.0, 0.0, 1.5, 0.0, -1.5].*diameter

# calculate the number of turbines
nturbines = length(turbine_x)

# set turbine base heights
turbine_z = zeros(nturbines)

# set turbine yaw values
turbine_yaw = zeros(nturbines)

# set turbine design parameters
rotor_diameter = zeros(nturbines) .+ diameter   # m
hub_height = zeros(nturbines) .+ 70.0           # m
cut_in_speed = zeros(nturbines) .+4.            # m/s
cut_out_speed = zeros(nturbines) .+25.          # m/s
rated_speed = zeros(nturbines) .+16.            # m/s
rated_power = zeros(nturbines) .+2.0E6          # W
generator_efficiency = zeros(nturbines) .+ 0.944

# Rotor swept area sample points (normalized by rotor radius). These arrays define which
# which points on the rotor swept area should be used to estimate the effective inflow
# wind speed for each wind turbine. Values of 0.0 are at the rotor hub, 1.0 is at the blade
tip. z is vertical, and y is horizontal. These points track the rotor yaw.
rotor_points_y = [0.0]
rotor_points_z = [0.0]

# set flow parameters
wind_speed = 8.0        # m/2
air_density = 1.1716    # kg/m^3
ambient_ti = 0.077      # %
shearexponent = 0.15
winddirections = [275.0*pi/180.0, 0.0, pi]          # radians
windspeeds = [wind_speed, wind_speed, wind_speed]   # m/s
windprobabilities = [1.0/3.0,1.0/3.0,1.0/3.0]       # %
ambient_tis = [ambient_ti, ambient_ti, ambient_ti]  # %
measurementheight = [hub_height[1], hub_height[1], hub_height[1]]   # m

# initialize the wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, 
measurementheight, air_density, ambient_tis, wind_shear_model)

```

### (2) Setting up the analysis models

A model set requires a Wake Deficit Model, Wake Deflection Model, Wake Combination Model, and a Local Turbulence Intensity Model
* Deficit Models: JensenTopHat, JensenCosine, MultiZone, GaussOriginal, GaussYaw, GaussYawVariableSpread, GaussSimple
* Deflection Models: GaussYawDeflection, GaussYawVariableSpreadDeflection, JiminezYawDeflection, MultizoneDeflection
* Combination Models: LinearFreestreamSuperposition, SumOfSquaresFreestreamSuperposition, SumOfSquaresLocalVelocitySuperposition, LinearLocalVelocitySuperposition
* Turbulence Models: LocalTIModelNoLocalTI, LocalTIModelMaxTI

The model set can be set up as follows:

Initialize power model (this is a simple power model based only on turbine design and is not accurate. For examples on how to use more accurate power models, look at the example optimization scripts)
```julia
power_model = ff.PowerModelPowerCurveCubic()
```

The user can define different power models for different wind turbines, but here we use the same power model for every turbine. The initialization of the power_models vector is important for optmization using algorithmic differentiation via the ForwardDiff.jl package.
```julia
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end
```

Initialize thrust model. The user can provide a complete thrust curve. See the example scripts for details on initializing them. The initialization of ct_models vector is important for optmization using algorithmic differentiation via the ForwardDiff.jl package.
```julia
ct_model = ff.ThrustModelConstantCt(0.65)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end
```

Set up wake and related models. Here we will use the default values provided in FLOWFarm.
However, it is important to use the correct model parameters. More information and references
are provided in the doc strings attached to each model.

The wake deficit model predicts the impact of wind turbines wake on the wind speed.
```julia
wakedeficitmodel = ff.GaussYaw()
```

The wake deflection model predicts the cross-wind location of the center of a wind turbine wake.
```julia
wakedeflectionmodel = ff.GaussYawDeflection()
```

The wake combination model defines how the predicted deficits in each wake should be combined to predict the total deficit at a point
```julia
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
```

The local turbulence intensity models can be used to estimate the local turbulence intensity at each wind turbine or point to provide more accurate input information to the wake and deflection models if applicable.
```julia
localtimodel = ff.LocalTIModelMaxTI()
```

Initialize model set. This is just a convenience container for the analysis models.
```julia
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
```

### (3) Running the analysis

Calculate AEP
```julia

    AEP = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
        rotor_sample_points_y=rotor_points_y, rotor_sample_points_z=rotor_points_z)

```
