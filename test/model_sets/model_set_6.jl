import FlowFarm; const ff = FlowFarm


ai = 1.0/3.0
ct = 0.689
wind_speed = 8.0
air_density = 1.1716  # kg/m^3
ambient_ti = 0.077
turbine_x = [-3.0, 0.0, 3.0, 0.0, 0.0, -1.5, 0.0, 1.5, 0.0].*80.0
nturbines = length(turbine_x)
turbine_y = [0.0, 3.0, 0.0, -3.0, 0.0, 0.0, 1.5, 0.0, -1.5]*80.0
turbine_z = zeros(nturbines)
turbine_yaw = zeros(nturbines)
turbine_ct = zeros(nturbines) .+ ct
turbine_ai = zeros(nturbines) .+ ai
winddirections = [275.0*pi/180.0, 0.0, pi]
windspeeds = [wind_speed, wind_speed, wind_speed]
windprobabilities = [1.0/3.0,1.0/3.0,1.0/3.0]
ambient_tis = [ambient_ti, ambient_ti, ambient_ti]
shearexponent = 0.15
turbine_inflow_velcities = zeros(nturbines) .+ wind_speed

rotor_diameter = zeros(nturbines) .+ 80.0
hub_height = zeros(nturbines) .+ 70.0
turbine_ct = zeros(nturbines)
turbine_ai = zeros(nturbines)
turbine_inflow_velcities = zeros(nturbines)
turbine_local_ti = zeros(nturbines)
turbine_yaw = zeros(nturbines)

measurementheight = [hub_height[1], hub_height[1], hub_height[1]]

cut_in_speed = zeros(nturbines) .+4.  # m/s
cut_out_speed = zeros(nturbines) .+25.  # m/s
rated_speed = zeros(nturbines) .+16.  # m/s
rated_power = zeros(nturbines) .+2.0E6  # W
generator_efficiency = zeros(nturbines) .+0.944


# rotor sample points
rotor_points_y = [0.0]
rotor_points_z = [0.0]

# load power curve
powerdata = readdlm("inputfiles/niayifar_vestas_v80_power_curve_observed.txt",  ',', skipstart=1)
velpoints = powerdata[:,1]
powerpoints = powerdata[:,2]*1E6

# initialize power model
power_model = ff.PowerModelPowerPoints(velpoints, powerpoints)
# power_model = Vector{typeof(power_model1)}(undef, nturbines)
# for i = 1:nturbines
#     power_model[i] = power_model1
# end

# load thrust curve
ctdata = readdlm("inputfiles/predicted_ct_vestas_v80_niayifar2016.txt",  ',', skipstart=1)
velpoints = ctdata[:,1]
ctpoints = ctdata[:,2]

# initialize thurst model
ct_model1 = ff.ThrustModelCtPoints(velpoints, ctpoints)
ct_model = Vector{typeof(ct_model1)}(undef, nturbines)
for i = 1:nturbines
    ct_model[i] = ct_model1
end

# initialize wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# get sorted indecies 
sorted_turbine_index = sortperm(turbine_x)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

wakedeficitmodel = ff.GaussYaw()
wakedeflectionmodel = ff.GaussYawDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelMaxTI()

model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)