using FLOWFarm; const ff = FLOWFarm
using SNOW
using SparseArrays

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

# set turbine design parameters (these values correspond to the Vestas V80 turbine)
rotordiameter = zeros(nturbines) .+ 80.0   # m
hubheight = zeros(nturbines) .+ 70.0           # m
cutinspeed = zeros(nturbines) .+ 4.0           # m/s
cutoutspeed = zeros(nturbines) .+ 25.0         # m/s
ratedspeed = zeros(nturbines) .+ 16.0          # m/s
ratedpower = zeros(nturbines) .+ 2.0E6         # W
generatorefficiency = ones(nturbines)

# get the sample points
nsamplepoints = 50
rotorsamplepointsy, rotorsamplepointsz = ff.rotor_sample_points(nsamplepoints, method="sunflower")

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

powermodel = ff.PowerModelPowerCurveCubic()

powermodels = Vector{typeof(powermodel)}(undef, nturbines)
for i = 1:nturbines
    powermodels[i] = powermodel
end

ctmodel = ff.ThrustModelConstantCt(0.65)
ctmodels = Vector{typeof(ctmodel)}(undef, nturbines)
for i = 1:nturbines
    ctmodels[i] = ctmodel
end

wakedeficitmodel = ff.GaussYawVariableSpread()

wakedeflectionmodel = ff.GaussYawDeflection()

wakecombinationmodel = ff.LinearLocalVelocitySuperposition()

localtimodel = ff.LocalTIModelMaxTI()

modelset = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

aep = ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz)

state_aeps = ff.calculate_state_aeps(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz,
    hours_per_year=365.25*24.0, weighted=true)

turbine_powers_by_direction = ff.calculate_state_turbine_powers(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz)
nothing