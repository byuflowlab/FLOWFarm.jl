# Tutorial

This tutorial covers the basics of FlowFARM. For more specifics refer to the [How-to guide](How_to.md).

This tutorial discusses how to do the following with FLOWFarm:
- (1) setting up a problem description
- (2) setting up an analysis model set
- (3) running analyses
- (4) setting up and running an optimization
- (5) calculating and visualizing a flow field

Details for setting up an optimization will depend heavily on the
optimization package you are using, your objective, and your design variables. Optimization
examples using various packages are provided in the example scripts located in the test 
directory.

## (1) Setting up the problem description

The problems description involves the physically description of the wind farm, the turbines, 
and the wind resource. While this tutorial uses the same design across all the wind turbines
and mostly equal properties across all wind flow states, all turbines and flow states can 
be unique.

For API demonstration purposes, we have directly assigned all values. However, values may 
be loaded from .csv and/or .yaml files.

### Set up the running environment
```@example 1
using FLOWFarm; const ff = FLOWFarm
using PyPlot; const plt = PyPlot
using VectorizedRoutines.Matlab: meshgrid
using SNOW
```

### Initialize the wind farm design
```@example 1

# set initial turbine x and y locations
turbinex = [-240.0, -240.0, -240.0, 0.0, 0.0, 0.0, 240.0, 240.0, 240.0]
turbiney = [-240.0, 0.0, 240.0, -240.0, 0.0, 240.0, -240.0, 0.0, 240.0]

# get the number of turbines
nturbines = length(turbinex)

# set turbine base heights
turbinez = zeros(nturbines)

# set turbine yaw values
turbineyaw = zeros(nturbines)

# set wind farm boundary parameters in meters (we won't really need this until we optimize)
boundarycenter = [0.0,0.0]
boundaryradius = hypot(300, 300) 
println("") # hide
```

### Initialize wind turbine design
```@example 1
# set turbine design parameters (these values correspond to the Vestas V80 turbine)
rotordiameter = zeros(nturbines) .+ 80.0   # m
hubheight = zeros(nturbines) .+ 70.0           # m
cutinspeed = zeros(nturbines) .+ 4.0           # m/s
cutoutspeed = zeros(nturbines) .+ 25.0         # m/s
ratedspeed = zeros(nturbines) .+ 16.0          # m/s
ratedpower = zeros(nturbines) .+ 2.0E6         # W
generatorefficiency = ones(nturbines)
println("") # hide
```

### Visualize the wind farm layout
```@example 1 
# initialize axis 
fig, ax = plt.subplots(1)

# plot layout using FLOWFarm
ff.plotlayout!(ax, turbinex, turbiney, rotordiameter)

# Label the axes
ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

# add farm boundary
circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
ax.add_patch(circle)

# set plot limits
ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01, aspect="equal")

plt.tight_layout() # hide
plt.savefig("initiallayout.png") # hide
```
![](initiallayout.png)

### Determine how to sample the flow field to determine effective inflow speeds
Rotor swept area sample points are normalized by the rotor radius. These arrays define which
which points on the rotor swept area should be used to estimate the effective inflow
wind speed for each wind turbine. Values of 0.0 are at the rotor hub, 1.0 is at the blade
tip, `z` is vertical, and `y` is horizontal. These points track the rotor when yawed. 
A single sample point will always be placed at the hub. More points can be arranged in 
either a grid pattern or a sunflower packing pattern with various options. 
See doc strings for more information.

```@example 1
# get the sample points
nsamplepoints = 50
rotorsamplepointsy, rotorsamplepointsz = ff.rotor_sample_points(nsamplepoints, method="sunflower")

# visualize the sample points
fig, ax = plt.subplots(1)
ff.plotrotorsamplepoints!(ax, rotorsamplepointsy, rotorsamplepointsz)
ax.set(xlabel="y/radius", ylabel="z/radius", aspect="equal")
plt.tight_layout() # hide
plt.savefig("rotorsamplepoints.png") # hide
println("") # hide
```
![](rotorsamplepoints.png)

### Setting up the wind resource
The wind resource determines the properties of the flowfield at all wind states. A wind 
state is any combination of wind speed, wind direction, turbulence intensity, etc...

```@example 1
# set flow parameters
windspeed = 8.0        # m/2
airdensity = 1.1716    # kg/m^3
ambientti = 0.1      # %
shearexponent = 0.15
ndirections = 5
winddirections = collect(range(0, 2*pi*(1-1/ndirections), length=ndirections))   # radians
windspeeds = ones(ndirections).*windspeed   # m/s
windprobabilities = ones(ndirections).*(1.0/ndirections)       # %
ambienttis = ones(ndirections).*ambientti  # %
measurementheight = ones(ndirections).*hubheight[1] # m

# initialize the wind shear model
windshearmodel = ff.PowerLawWindShear(shearexponent)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, 
measurementheight, airdensity, ambienttis, windshearmodel)

# visualize the wind resource
ff.plotwindresource!(windresource)
plt.tight_layout() # hide
plt.savefig("windresource.png") # hide
```
![](windresource.png)

## (2) Setting up the analysis models

A model set requires a Wake Deficit Model, Wake Deflection Model, Wake Combination Model, and a Local Turbulence Intensity Model. There are several options for each model type. To facilitate research studies, any of the models in each type can be used with any of the models in any other type. However, behavior is not guaranteed. It is recommended that common, validated, model combinations be used in most cases.

Model types and options are:
* Deficit Models: JensenTopHat, JensenCosine, MultiZone, GaussOriginal, GaussYaw, GaussYawVariableSpread, GaussSimple
* Deflection Models: GaussYawDeflection, GaussYawVariableSpreadDeflection, JiminezYawDeflection, MultizoneDeflection
* Combination Models: LinearFreestreamSuperposition, SumOfSquaresFreestreamSuperposition SumOfSquaresLocalVelocitySuperposition, LinearLocalVelocitySuperposition
* Turbulence Models: LocalTIModelNoLocalTI, LocalTIModelMaxTI

The model set can be set up as follows:

Initialize power model (this is a simple power model based only on turbine design and is not very accurate. For examples on how to use more accurate power models, look at the example optimization scripts)
```@example 1
powermodel = ff.PowerModelPowerCurveCubic()
```

The user can define different power models for different wind turbines, but here we use the same power model for every turbine. The initialization of the power_models vector is important for optmization using algorithmic differentiation via the ForwardDiff.jl package.
```@example 1
powermodels = Vector{typeof(powermodel)}(undef, nturbines)
for i = 1:nturbines
    powermodels[i] = powermodel
end
```

Initialize thrust model(s). The user can provide a complete thrust curve. See the example scripts for details on initializing them. The initialization of the ct models vector is important for optmization using algorithmic differentiation via the ForwardDiff.jl package.
```@example 1
ctmodel = ff.ThrustModelConstantCt(0.65)
ctmodels = Vector{typeof(ctmodel)}(undef, nturbines)
for i = 1:nturbines
    ctmodels[i] = ctmodel
end
```

Set up wake and related models. Here we will use the default values provided in FLOWFarm.
However, it is important to use the correct model parameters. More information and references
are provided in the doc strings attached to each model.

The wake deficit model predicts the impact of wind turbines wake on the wind speed.
```@example 1
wakedeficitmodel = ff.GaussYawVariableSpread()
```

The wake deflection model predicts the cross-wind location of the center of a wind turbine wake.
```@example 1
wakedeflectionmodel = ff.GaussYawDeflection()
```

The wake combination model defines how the predicted deficits in each wake should be combined to predict the total deficit at a point
```@example 1
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
```

The local turbulence intensity models can be used to estimate the local turbulence intensity at each wind turbine or point to provide more accurate input information to the wake and deflection models if applicable.
```@example 1
localtimodel = ff.LocalTIModelMaxTI()
```

Initialize model set. This is just a convenience container for the analysis models.
```@example 1
modelset = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
println("") # hide
```

## (3) Running the analysis

Now that the wind farm and analysis models have been defined, we can calculate AEP. The output
is in Watt-hours.

```@example 1
aep = ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz)
```

We can also get the AEP in each direction using the following.

```@example 1
state_aeps = ff.calculate_state_aeps(turbinex, turbiney, turbinez, rotordiameter,
        hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
        cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
        rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz, 
        hours_per_year=365.25*24.0, weighted=true)
```

If we instead set `weighted=false` then we would get the power in each direction in Watts.

If we want to get the individual turbine powers in each directions, we use the following.

```@example 1
turbine_powers_by_direction = ff.calculate_state_turbine_powers(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz)
```
The output shows each turbine power in an array that is ndirections by nturbines.

## (4) setting up and running an optimization
FLOWFarm is specifically designed for efficient optimization using gradient-based optimization
methods. Besides the steps outlined above, we need to define the following before we can run 
and optimization:

- (1) Optimization related variables
- (1) A container for non-differentiated parameters
- (2) Objective function 
- (3) Constraint function(s) 
- (4) Optimization tool specific items

In this tutorial we demonstrate optimizing using the IPOPT algorithms via SNOW.jl for simplicity.

First, set up optimization related variables. We will have two constraints, one to keep 
turbines from getting too close to each other (spacing), and the other to keep turbines 
inside the desired area (boundary). FLOWFarm provides several different ways of handling 
boundary constraints, including concave boundaries. However, for this tutorial we will use 
a simple circular boundary.

```@example 1
# scale objective derivatives to be between 0 and 1
objectivescale = 1E-6

# scale boundary constraint derivatives to be between 0 and 1
constraintscaleboundary = 1.0E-3

# scale spacing constraint derivatives to be between 0 and 1
constraintscalespacing = 1.0

# set the minimum spacing between turbines 
minimumspacing = 160.0
println("") # hide
```

Next, set up a container for non-differentiated parameters

```@example 1
# set up a struct for use in optimization functions
mutable struct params_struct{}
    modelset
    rotorsamplepointsy
    rotorsamplepointsz
    turbinez
    ambientti
    rotordiameter
    boundarycenter
    boundaryradius
    objectivescale
    constraintscaleboundary
    constraintscalespacing
    minimumspacing
    hubheight
    turbineyaw
    ctmodels
    generatorefficiency
    cutinspeed
    cutoutspeed
    ratedspeed
    ratedpower
    windresource
    powermodels
end

params = params_struct(modelset, rotorsamplepointsy, rotorsamplepointsz, turbinez, ambientti, 
    rotordiameter, boundarycenter, boundaryradius, objectivescale, constraintscaleboundary,
    constraintscalespacing, minimumspacing, hubheight, turbineyaw, 
    ctmodels, generatorefficiency, cutinspeed, cutoutspeed, ratedspeed, ratedpower, 
    windresource, powermodels)
println("") # hide
```

Now we are ready to set up wrapper functions for the objective and constraints.

```@example 1
# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant params
    boundarycenter = params.boundarycenter
    boundaryradius = params.boundaryradius
    constraintscaleboundary = params.constraintscaleboundary

    # find the number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(boundarycenter, boundaryradius, turbinex, turbiney).*constraintscaleboundary
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # include relevant params
    rotordiameter = params.rotordiameter
    constraintscalespacing = params.constraintscalespacing
    minimumspacing = params.minimumspacing

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # get and return spacing distances
    return constraintscalespacing.*(minimumspacing .- ff.turbine_spacing(turbinex,turbiney))
end

# set up aep wrapper function
function aep_wrapper(x, params)

    # include relevant params
    turbinez = params.turbinez
    rotordiameter = params.rotordiameter
    hubheight = params.hubheight
    turbineyaw =params.turbineyaw
    ctmodels = params.ctmodels
    generatorefficiency = params.generatorefficiency
    cutinspeed = params.cutinspeed
    cutoutspeed = params.cutoutspeed
    ratedspeed = params.ratedspeed
    ratedpower = params.ratedpower
    windresource = params.windresource
    powermodels = params.powermodels
    modelset = params.modelset
    rotorsamplepointsy = params.rotorsamplepointsy
    rotorsamplepointsz = params.rotorsamplepointsy
    objectivescale = params.objectivescale

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines] 
    turbiney = x[nturbines+1:end]

    # calculate AEP
    aep = objectivescale*ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
                hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
                cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
                rotor_sample_points_y=rotorsamplepointsy,rotor_sample_points_z=rotorsamplepointsz)
    
    # return the AEP
    return aep
end

# set up optimization problem wrapper function
function wind_farm_opt!(g, x, params)

    nturbines = Int(length(x)/2)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x, params)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x, params)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[1:(end-nturbines)] = spacing_con[:]
    g[end-nturbines+1:end] = boundary_con[:]
    
    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    obj = -aep_wrapper(x, params)[1]
    
    return obj
end
println("") # hide
```

Because the optimizer will need to call the objective function without knowing about the params, we
need to set up a method that will know the params values by default.

```@example 1
# generate objective function wrapper
obj_func!(g, x) = wind_farm_opt!(g, x, params)
```

Next we set up the optimizer.

```@example 1
# initialize design variable vector
x0 = [copy(turbinex);copy(turbiney)]

# set general lower and upper bounds for design variables
lx = zeros(length(x0)) .- boundaryradius
ux = zeros(length(x0)) .+ boundaryradius

# set general lower and upper bounds for constraints
ng = Int(nturbines + (nturbines)*(nturbines - 1)/2)
lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines)]
ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines)]

# IPOPT options
ip_options = Dict(
    "max_iter" => 50,
    "tol" => 1e-6
)
solver = IPOPT(ip_options)

# if using SNOPT, you can do the following instead:
# snopt_opt = Dict(
#    "Derivative option" => 1,
#    "Major optimality tolerance" => 1e-4,
# )
# solver = SNOPT(options=snopt_opt)

# initialize SNOW options
options = Options(solver=solver, derivatives=ForwardAD())  # choose AD derivatives
println("") # hide
```

Now that the optimizer is set up, we are ready to optimize and check the results.

```@example 1
# optimize
t1 = time() # start time
xopt, fopt, info, out = minimize(obj_func!, x0, ng, lx, ux, lg, ug, options)
t2 = time() # end time
clk = t2-t1 # approximate run time

# get final aep
aepfinal = -fopt/objectivescale

# print optimization results 
println("Finished in : ", clk, " (s)")
println("info: ", info)
println("Initial AEP: ", aep)
println("Final AEP: ", aepfinal)
println("AEP improvement (%) = ", 100*(aepfinal - aep)/aep) 

# extract final turbine locations
turbinexopt = copy(xopt[1:nturbines])
turbineyopt = copy(xopt[nturbines+1:end])

# initialize figure and axes object
fig, ax = plt.subplots(1)

# plot layout using FLOWFarm
ff.plotlayout!(ax, turbinexopt, turbineyopt, rotordiameter)

# label the axes
ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

# and the wind farm boundary
circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
ax.add_patch(circle)

# set limits on the plot region
ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01, aspect="equal")

plt.tight_layout() # hide
plt.savefig("optlayout.png") # hide
```
![](optlayout.png)

## (5) Calculating and visualizing a flow field

It is helpful to visualize the whole flow-field, not just the turbine powers. Here we will 
visualize the flow field for a 2D horizontal cross-section at the hub height. FLOWFarm is 
capable of generating flow fields in 1D, 2D, and 3D.

```@example 1
# define how many points should be in the flow field
xres = 1000
yres = 1000
zres = 1

# define flow field domain
maxy = boundaryradius*1.5
miny = -boundaryradius*1.5
maxx = boundaryradius*1.5
minx = -boundaryradius*1.5

# set up point grid for flow field
xrange = minx:(maxx-minx)/xres:maxx
yrange = miny:(maxy-miny)/yres:maxy
zrange = hubheight[1]

# run flowfarm 
ffvelocities = ff.calculate_flow_field(xrange, yrange, zrange,
    modelset, turbinexopt, turbineyopt, turbinez, turbineyaw,
    rotordiameter, hubheight, ctmodels, rotorsamplepointsy, rotorsamplepointsz,
    windresource, wind_farm_state_id=5)

# visualize 

# initialize figure
fig, ax = plt.subplots(1)

# generate meshgrid from ranges for passing to pyplots
xg, yg = meshgrid(collect(xrange), collect(yrange))

# plot as filled contours
cs = ax.contourf(xg, yg, ffvelocities[1,:,:], cmap="Blues_r")

# add colorbar 
cbar = ax.figure.colorbar(cs, ax=ax, label="Wind Speed (m/s)", orientation="vertical")

# label the axes
ax.set(xlabel="Easting (m)", ylabel="Northing (m)", aspect="equal")

# add the wind farm boundary
circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
ax.add_patch(circle)

plt.savefig("flowfield.png") # hide
```
![](flowfield.png)
