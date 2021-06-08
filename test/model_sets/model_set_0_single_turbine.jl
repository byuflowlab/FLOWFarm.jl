import FLOWFarm; const ff = FLOWFarm

turbine_x = [0.0]
turbine_y = [0.0]
turbine_z = [0.0]

# set turbine design parameters
rotor_diameter = zeros(nturbines) .+ diam # m
hub_height = zeros(nturbines) .+ 90.0   # m
cut_in_speed = zeros(nturbines) .+3.  # m/s
cut_out_speed = zeros(nturbines) .+25.  # m/s
rated_speed = zeros(nturbines) .+11.4  # m/s
rated_power = zeros(nturbines) .+5.0E6  # W
generator_efficiency = zeros(nturbines) .+ 0.944

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

sorted_turbine_index = [1]

wind_speed = 8.1
air_density = 1.1716  # kg/m^3
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheights = [hub_height[1]]
wtvelocities = [wind_speed]
ambient_tis = [0.1]
shearexponent = 0.0
wind_shear_model = ff.PowerLawWindShear(shearexponent)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

# set up wake and related models
wakedeficitmodel = ff.GaussYawVariableSpread()
wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)