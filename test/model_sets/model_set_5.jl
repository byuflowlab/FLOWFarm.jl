import FlowFarm; const ff = FlowFarm

rotor_diameter = 0.57
hub_height = 0.7
yaw = 0.0
ct = 0.81
constcp = 0.8
generator_efficiency = 0.944
ai = 1.0/3.0
wind_speed = 8.1
air_density = 1.1716  # kg/m^3
ambient_ti = 0.1
turbine_x = [0.0]
turbine_y = [0.0]
turbine_z = [0.0]
turbine_yaw = [0.0]
turbine_ct = [ct]
turbine_ai = [ai]
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheights = [hub_height]
shearexponent = 0.15
turbine_inflow_velcities = [wind_speed]
ambient_tis = [ambient_ti]
nturbines = 1

ct_model = ff.ThrustModelConstantCt(ct)
power_model = ff.PowerModelConstantCp(constcp)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end
wind_shear_model = [ff.PowerLawWindShear(shearexponent)]

cut_in_speed = 0.0
cut_out_speed = 25.0
rated_speed = 12.0
rated_power = 1.0176371581904552e6

turbine1 = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], ct_model, power_model)
turbine_definitions = [turbine1]
turbine_definition_ids = [1]
sorted_turbine_index = [1]

windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities, [0.0], [ambient_ti],sorted_turbine_index)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

loc = [7.0*rotor_diameter, 0.0, hub_height]
alpha = 0.1
wakedeficitmodel = ff.JensenTopHat(alpha)
horizontal_spread_rate = alpha
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

ms1 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
# pd1 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
