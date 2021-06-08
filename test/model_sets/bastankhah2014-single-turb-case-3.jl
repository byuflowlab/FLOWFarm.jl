import FLOWFarm; const ff = FLOWFarm

turbine_x = [0.0]
turbine_y = [0.0]
turbine_z = [0.0]
turbine_yaw = [0.0]
diam = 80.0
nturbines = 1

# set turbine design parameters
rotor_diameter = zeros(nturbines) .+ diam # m
hub_height = zeros(nturbines) .+ 70.0   # m
cut_in_speed = zeros(nturbines) .+4.  # m/s
cut_out_speed = zeros(nturbines) .+25.  # m/s
rated_speed = zeros(nturbines) .+16.  # m/s
rated_power = zeros(nturbines) .+2.0E6  # W
generator_efficiency = zeros(nturbines) .+0.944

# load power curve
powerdata = readdlm("./inputfiles/niayifar_vestas_v80_power_curve_observed.txt",  ',', skipstart=1)
velpoints = powerdata[:,1]
powerpoints = powerdata[:,2]*1E6

# initialize power model
power_model = ff.PowerModelPowerPoints(velpoints, powerpoints)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# load thrust curve
ctdata = readdlm("./inputfiles/mfg-ct-vestas-v80-niayifar2016.txt",  ',', skipstart=1)
velpoints = ctdata[:,1]
ctpoints = ctdata[:,2]

# initialize thurst model
ct_model = ff.ThrustModelConstantCt(0.8)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end

sorted_turbine_index = [1]

wind_speed = 9.0
air_density = 1.1716  # kg/m^3
# air_density = 1.2
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheights = [hub_height[1]]
wtvelocities = [wind_speed]
ambient_tis = [0.094]
shearexponent = 0.0
wind_shear_model = ff.PowerLawWindShear(shearexponent)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

# set up wake and related models
wakedeficitmodel = ff.GaussOriginal(0.04)
# wakedeficitmodel = ff.GaussYawVariableSpread()
wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)