import FlowFarm; const ff = FlowFarm

diam = 80.0
data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)
turbine_x = data[:, 1].*7.0*diam
nturbines = length(turbine_x)

rotor_diameter = zeros(nturbines).+diam
generator_efficiency = zeros(nturbines).+0.944
cut_in_speed = zeros(nturbines).+4.  # m/s
cut_out_speed = zeros(nturbines).+25.  # m/s
rated_speed = zeros(nturbines).+16.  # m/s
rated_power = zeros(nturbines).+2.0E6  # W
air_density = 1.1716  # kg/m^3

cp = 0.8
power_model = ff.PowerModelConstantCp(cp)





#
#
#
# hub_height = 70.0
# yaw = 0.0
# ct = 0.689
#
#
#
#
#
# ai = 1.0/3.0
# wind_speed = 12.0
#
# ambient_ti = 0.1
#
# nturbines = length(turbine_x)
# turbine_y = zeros(nturbines)
# turbine_z = zeros(nturbines)
# turbine_yaw = zeros(nturbines)
# turbine_ct = zeros(nturbines) .+ ct
# turbine_ai = zeros(nturbines) .+ ai
# winddirections = [270.0*pi/180.0]
# windspeeds = [wind_speed]
# windprobabilities = [1.0]
# measurementheight = [hub_height]
# ambient_tis = [ambient_ti]
# shearexponent = 0.15
# turbine_inflow_velcities = zeros(nturbines) .+ wind_speed
#
# # rotor sample points
# rotor_points_y = [0.0]
# rotor_points_z = [0.0]
#
# ct_model = ff.ThrustModelConstantCt(ct)
#
# wind_shear_model = ff.PowerLawWindShear(shearexponent)
#
# turbine1 = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], ct_model, power_model)
#
# turbine_definitions = [turbine1 for i in 1:nturbines]
# sorted_turbine_index = [i for i  in 1:nturbines]
# turbine_definition_ids = ones(Int, nturbines)
#
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities, zeros(nturbines), (zeros(nturbines).+ambient_ti),sorted_turbine_index)
# windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, [wind_shear_model])
#
# alpha = 0.1
# wakedeficitmodel = ff.JensenTopHat(alpha)
# horizontal_spread_rate = alpha
# wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
# wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
# localtimodel = ff.LocalTIModelNoLocalTI()
#
# ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
