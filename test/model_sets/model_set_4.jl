import FlowFarm; const ff = FlowFarm

rotor_diameter = 80.0
hub_height = 70.0
yaw = 0.0

cut_in_speed = 4.  # m/s
cut_out_speed = 25.  # m/s
rated_speed = 16.  # m/s
rated_power = 2.0E6  # W
generator_efficiency = 0.944

ai = 1.0/3.0
ct = 0.689
wind_speed = 8.0
air_density = 1.1716  # kg/m^3
ambient_ti = 0.077

data = readdlm("inputfiles/horns_rev_locations.txt",  ',', skipstart=1)
turbine_x = data[:, 1].*diam
nturbines = length(turbine_x)
turbine_y = data[:, 2].*diam
turbine_z = zeros(nturbines)

rotor_diameter = zeros(nturbines) .+ 80.0
hub_height = zeros(nturbines) .+ 70.0
sorted_turbine_index = sortperm(turbine_x)
turbine_ct = zeros(nturbines)
turbine_ai = zeros(nturbines)
turbine_inflow_velcities = zeros(nturbines)
turbine_local_ti = zeros(nturbines)
turbine_yaw = zeros(nturbines)
ambient_ti = 0.077
ai = 1.0/3.0
turbine_ai = zeros(nturbines) .+ ai

localtimodel = ff.LocalTIModelMaxTI()

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

rotor_sample_points_y = [0.0]
rotor_sample_points_z = [0.0]

winddirections = [270.0*pi/180.0]
wind_speed = 8.0
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheight = [hub_height[1]]
air_density = 1.1716  # kg/m^3
ambient_tis = [ambient_ti]
shearexponent = 0.15
# initialize wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)


wakedeficitmodel = ff.GaussYaw(0.022, 0.022, 2.32, 0.154, 1.0)
wakedeflectionmodel = ff.GaussYawDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelMaxTI()
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)



# yaw = 0.0
#
# cut_in_speed = 4.  # m/s
# cut_out_speed = 25.  # m/s
# rated_speed = 16.  # m/s
# rated_power = 2.0E6  # W
# generator_efficiency = 0.944
#
#
#
#
# ambient_ti = 0.077
#
# turbine_z = zeros(nturbines)
# turbine_yaw = zeros(nturbines)
# turbine_ct = zeros(nturbines) .+ ct
#
#
#
# turbine_inflow_velcities = zeros(nturbines) .+ wind_speed
#
# # rotor sample points
# rotor_points_y = [0.0]
# rotor_points_z = [0.0]
#
# # load power curve
# powerdata = readdlm("inputfiles/niayifar_vestas_v80_power_curve_observed.txt",  ',', skipstart=1)
# velpoints = powerdata[:,1]
# powerpoints = powerdata[:,2]*1E6
#
# # initialize power model
# power_model = ff.PowerModelPowerPoints(velpoints, powerpoints)
#
#
#
#
#
# # initialize turbine definition
# turbine1 = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], ct_model, power_model)
#
# turbine_definitions = [turbine1 for i in 1:nturbines]
# sorted_turbine_index = sortperm(turbine_x)
# turbine_definition_ids = ones(Int, nturbines)
#
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities, zeros(nturbines), (zeros(nturbines).+ambient_ti),sorted_turbine_index)
#
#
#
#
#
# pd4 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
