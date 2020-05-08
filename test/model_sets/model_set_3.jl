import FlowFarm; const ff = FlowFarm


rotor_diameter = 120.0

turbine_x = [6.2874587,6.6267154,8.1331868,9.2637805,11.427071,12.744123,16.296166,17.480180,
              18.136485,18.338263,19.320500,20.045544,21.229559,22.610909,23.597587,26.162952,
              30.109666,28.586395,31.688352,33.005404,34.056381,34.648388,35.240395,35.635066,
              36.424409,37.608424,38.002128,37.870057,39.714819,40.223450,41.752474,42.344481,
              43.790129,45.896524,46.685867,47.672545,48.002919,49.186933,49.448567,51.027253,
              51.224588,51.816596,54.973967,56.552653,57.736667,59.448391,59.710024,61.550344,
              62.078053,63.262067,65.827432,67.208782,69.971482,71.155496,72.734182,74.115532,
              78.392620,87.408726,23.015521,23.281596,23.547672,23.547672,23.813747,23.813747,
              23.680710,23.281596,23.015521,22.084257,21.419069,20.886918,15.299335,13.037694,
              11.707317,10.643016,7.8492239,6.2527716,4.7893570,3.1929047,1.5964523,0.0000000,
              9.8447894,4.3902439,19.689579,22.882483,24.611973,26.341463,28.203991,30.199557,
              31.929047,42.172949,83.015521,87.538803,92.461197,82.483370,77.560976,67.716186,
              62.660754,57.605322,52.815965,47.893570,42.838137,37.915743,32.860310,37.516630,
              37.250554,30.864745,26.208426,22.217295,20.221729,64.789357,75.565410] .* rotor_diameter

turbine_y = [175.20193, 155.99196, 170.58074, 127.23912, 179.91177, 157.78818, 97.594207, 92.225670,
             179.86703, 86.857133, 81.367953, 75.941108, 170.50689, 52.766924, 36.318468, 95.334947,
             80.862267, 114.17836, 73.525267, 100.99428, 138.17140, 57.171113, 93.724386, 50.065058,
             42.790393, 86.857133, 129.12097, 20.729884, 80.199999, 124.62390, 72.879054, 101.26271,
             65.651413, 57.911772, 111.28398, 50.311316, 87.930840, 106.98915, 42.532108, 81.153063,
             34.702792, 102.87327, 98.746205, 67.575138, 95.021783, 60.819730, 80.056986, 90.083626,
             54.109059, 74.062120, 85.604475, 68.246206, 81.443859, 62.564504, 45.049653, 77.511902,
             73.538837, 58.573441, 0.0000000, 5.0670241, 10.254692, 15.563003, 20.750670, 25.817694,
             31.005362, 41.501340, 46.689008, 58.873995, 65.026810, 70.455764, 102.90885, 112.80161,
             117.62735, 122.45308, 131.98391, 136.68901, 141.51475, 146.21984, 151.04558, 155.63003,
             165.76408, 180.00000, 175.17426, 165.76408, 161.17962, 156.59517, 151.89008, 147.30563,
             142.84182, 120.16086, 69.973190, 66.353887, 63.337802, 54.048257, 49.584450, 40.536193,
             36.072386, 31.487936, 27.024129, 22.560322, 18.096515, 13.391421, 9.1689008, 28.109920,
             35.589812, 107.37265, 120.88472, 131.62198, 137.17158, 47.412869, 57.305630] .* rotor_diameter


hub_height = 81.6
yaw = 0.0

cut_in_speed = 3.5  # m/s
cut_out_speed = 25.  # m/s
rated_speed = 11.4  # m/s
rated_power = 3.6E6  # W

ct = 0.689
cp = 0.48
generator_efficiency = 0.944
ai = 1.0/3.0

winddirections = deg2rad.(range(0,stop=350,length=36))
# winddirections = zeros(36).+deg2rad(270.0)
windspeeds = zeros(36).+(9.23*1.05)
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
power_model = ff.PowerModelPowerCurveCubic(rated_speed,rated_power)
wind_shear_model = ff.PowerLawWindShear(shearexponent)

turbine1 = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], ct_model, power_model)

turbine_definitions = [turbine1 for i in 1:nturbines]
sorted_turbine_index = [i for i  in 1:nturbines]
turbine_definition_ids = ones(Int, nturbines)

windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, [wind_shear_model])

windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines), (zeros(nturbines).+ambient_ti))
windfarmstate_array = [windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
                            windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
                            windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
                            windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate]

for i = 1:length(winddirections)
              windfarmstate_array[i] = ff.SingleWindFarmState(i, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines), (zeros(nturbines).+ambient_ti))
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
