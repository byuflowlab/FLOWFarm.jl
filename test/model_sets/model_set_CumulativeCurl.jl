import FLOWFarm; const ff = FLOWFarm

# set initial turbine x and y locations
diam = 126.0
data = readdlm("inputfiles/layout_38turb_round.txt",  ' ', skipstart=1)
turbine_x = data[:, 1].*diam
nturbines = length(turbine_x)
turbine_y = data[:, 2].*diam
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
cut_in_speed = zeros(nturbines) .+ 3.  # m/s 4.0
cut_out_speed = zeros(nturbines) .+ 25.  # m/s
rated_speed = zeros(nturbines) .+ 11.4  # m/s
rated_power = zeros(nturbines) .+ 5.0E6  # W
generator_efficiency = zeros(nturbines) .+ 1.0

# rotor swept area sample points (normalized by rotor radius)
nrotorpoints = 3
rotor_points_y, rotor_points_z = ff.rotor_sample_points(nrotorpoints, method="grid", pradius=0.5)

# set flow parameters
winddirections = [270.0*pi/180]
windspeeds = [8.0]
windprobabilities = [1.0]
nstates = length(windspeeds)

air_density = 1.225  # kg/m^3
ambient_ti = 0.06
shearexponent = 0.0
ambient_tis = zeros(nstates) .+ ambient_ti
measurementheight = zeros(nstates) .+ hub_height[1]

# load power curve
powerdata = readdlm("inputfiles/NREL5MWCPCT_FLORIS_V3.csv",',',skipstart=1)
velpoints = powerdata[:,1]
powerpoints = powerdata[:,2]

# initialize power model
power_model = ff.PowerModelCpPoints(velpoints, powerpoints)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# load thrust curve
ctdata = readdlm("inputfiles/NREL5MWCPCT_FLORIS_V3.csv",',',skipstart=1)
velpoints = ctdata[:,1]
ctpoints = ctdata[:,3]

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
wakedeficitmodel = ff.CumulativeCurl()

wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

# initialize model set
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
