import FlowFarm; const ff = FlowFarm

# set initial turbine x and y locations
diam = 126.4
data = readdlm("inputfiles/layout_amalia.txt",  ' ', skipstart=1)
turbine_x = data[:, 1]
nturbines = length(turbine_x)
turbine_y = data[:, 2]
turbine_z = zeros(nturbines)

turbine_x = turbine_x .- turbine_x[1]
turbine_y = turbine_y .- turbine_y[1]

# calculate the number of turbines
nturbines = length(turbine_x)

# set turbine base heights
turbine_z = zeros(nturbines) .+ 0.0

# set turbine yaw values
turbine_yaw = zeros(nturbines)

# set turbine design parameters
rotor_diameter = zeros(nturbines) .+ diam # m
hub_height = zeros(nturbines) .+ 90.0   # m
cut_in_speed = zeros(nturbines) .+3.  # m/s
cut_out_speed = zeros(nturbines) .+25.  # m/s
rated_speed = zeros(nturbines) .+11.4  # m/s
rated_power = zeros(nturbines) .+5.0E6  # W
generator_efficiency = zeros(nturbines) .+ 0.944

# rotor swept area sample points (normalized by rotor radius)
rotor_points_y = [0.0]
rotor_points_z = [0.0]

# set flow parameters
data = readdlm("inputfiles/windrose_amalia_72dirs.txt",  ' ', skipstart=1)
winddirections = data[:, 1].*pi/180.0
windspeeds = data[:,2]
windprobabilities = data[:, 3]
nstates = length(windspeeds)

air_density = 1.1716  # kg/m^3
ambient_ti = 0.077
shearexponent = 0.15
ambient_tis = zeros(nstates) .+ ambient_ti
measurementheight = zeros(nstates) .+ hub_height[1]

# load power curve
cpctdata = readdlm("inputfiles/NREL5MWCPCT.txt",  ' ', skipstart=1)
velpoints = cpctdata[:,1]
cppoints = cpctdata[:,2]

# initialize power model
power_model = ff.PowerModelCpPoints(velpoints, cppoints)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# load thrust curve
ctpoints = cpctdata[:,3]

# initialize thurst model
ct_model = ff.ThrustModelCtPoints(velpoints, ctpoints)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end

# initialize wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# get sorted indecies 
sorted_turbine_index = sortperm(turbine_x)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

# set up wake and related models
wakedeficitmodel = ff.JensenTopHat()

wakedeflectionmodel = ff.JiminezYawDeflection()
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

# initialize model set
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

