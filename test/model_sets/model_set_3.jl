import FlowFarm; const ff = FlowFarm

# based on https://backend.orbit.dtu.dk/ws/portalfiles/portal/146793996/Pena18_wes_3_191_2018.pdf

rotor_diameter = 120.0

locdata = readdlm("inputfiles/AnholtOffshoreWindFarmLocations.txt", ' ', skipstart=1)
turbine_x = locdata[:, 1]*rotor_diameter
turbine_y = locdata[:, 2]*rotor_diameter

hub_height = 81.6
yaw = 0.0

cut_in_speed = 3.5  # m/s
cut_out_speed = 25.  # m/s
rated_speed = 12.0  # m/s
rated_power = 3.6E6  # W

# load cp data
cpdata = readdlm("inputfiles/power_curve_siemens_swp3.6-120.txt",  ',', skipstart=1)
vel_points = cpdata[:, 1]
power_points = cpdata[:, 2]

ct = 0.689
cp = 0.48
generator_efficiency = 0.944
ai = 1.0/3.0

winddirections = deg2rad.(range(0,stop=350,length=36))
# winddirections = zeros(36).+deg2rad(270.0)
windspeeds = zeros(36).+(9.23)
windprobabilities = [1.3,1.1,1.3,1.1,1.0,1.3,1.7,2.2,2.1,3.3,3.5,4.1,4.5,5.0,
                            3.6,3.3,2.8,2.6,2.8,3.7,4.7,4.3,3.6,3.5,3.1,3.6,4.0,3.6,
                            4.2,3.8,1.8,1.6,1.6,1.4,1.2,1.3]
windprobabilities = windprobabilities./(sum(windprobabilities))

air_density = 1.1716  # kg/m^3
ambient_ti = 0.1
nturbines = length(turbine_x)
turbine_z = zeros(nturbines)
turbine_yaw = zeros(nturbines)
turbine_ct = zeros(nturbines) .+ ct
turbine_ai = zeros(nturbines) .+ ai

measurementheight = ones(length(winddirections)).*hub_height
ambient_tis = ones(length(winddirections)).*ambient_ti
shearexponent = 0.15
turbine_inflow_velcities = zeros(nturbines) .+ 9.21

# rotor sample points
rotor_points_y = [0.0]
rotor_points_z = [0.0]

ct_model = ff.ThrustModelConstantCt(ct)
# power_model = ff.PowerModelConstantCp(cp)
power_model = ff.PowerModelPowerPoints(vel_points, power_points)
wind_shear_model = ff.PowerLawWindShear(shearexponent)

turbine1 = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], ct_model, power_model)

turbine_definitions = [turbine1 for i in 1:nturbines]
sorted_turbine_index = [i for i  in 1:nturbines]
turbine_definition_ids = ones(Int, nturbines)

windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, [wind_shear_model])

windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities, zeros(nturbines), (zeros(nturbines).+ambient_ti),sorted_turbine_index)
windfarmstate_array = [windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
                            windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
                            windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
                            windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate]

for i = 1:length(winddirections)
              windfarmstate_array[i] = ff.SingleWindFarmState(i, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities, zeros(nturbines), (zeros(nturbines).+ambient_ti),sorted_turbine_index)
end
#
alpha = 0.04
wakedeficitmodel = ff.JensenTopHat(alpha)
horizontal_spread_rate = alpha
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

ms3 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
pd3 = ff.WindFarmProblemDescription(windfarm, windresource, windfarmstate_array)
