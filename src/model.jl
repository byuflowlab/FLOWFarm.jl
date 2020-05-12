# using FlowFarm
# using CCBlade
# using PyPlot
# using FLOWMath
# using Statistics

const ff=FlowFarm

r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
              28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
              56.1667, 58.9000, 61.6333]
chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
                  3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
                  6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

Rhub = 1.5
Rtip = 63.0
B = 3
pitch = 0.0
precone = 2.5*pi/180

af_path = "/Users/ningrsrch/Dropbox/Projects/waked-loads/5MW_AFFiles_julia"
path1 = af_path*"/Cylinder1.dat"
path2 = af_path*"/Cylinder2.dat"
path3 = af_path*"/DU40_A17.dat"
path4 = af_path*"/DU35_A17.dat"
path5 = af_path*"/DU30_A17.dat"
path6 = af_path*"/DU25_A17.dat"
path7 = af_path*"/DU21_A17.dat"
path8 = af_path*"/NACA64_A17.dat"

af1 = af_from_files(path1)
af2 = af_from_files(path2)
af3 = af_from_files(path3)
af4 = af_from_files(path4)
af5 = af_from_files(path5)
af6 = af_from_files(path6)
af7 = af_from_files(path7)
af8 = af_from_files(path8)

af = [af1,af2,af3,af4,af5,af6,af7,af8]

af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]
airfoils = af[af_idx]

pitch = 0.
Rhub = 1.5
hubHt = 90.

speeds = range(3.,stop=25.,length=23)
omegas = [6.972,7.183,7.506,7.942,8.469,9.156,10.296,11.431,11.89,
                    12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,
                    12.1,12.1,12.1]
pitches = -1.0.*[0.,0.,0.,0.,0.,0.,0.,0.,0.,3.823,6.602,8.668,10.45,12.055,
                        13.536,14.92,16.226,17.473,18.699,19.941,21.177,22.347,
                        23.469]


rotor_diameter = 126.4
hub_height = 90.0
yaw = 0.0
ct = 8.0/9.0
cp = 0.42

cut_in_speed = 3.0
cut_out_speed = 25.0
rated_speed = 11.4
rated_power = 5e6

generator_efficiency = 0.944
ai = 1.0/3.0
wind_speed = 11.
air_density = 1.225  # kg/m^3
nturbines = 1
# turbine_y = zeros(nturbines)

turbine_x = [0.0]
turbine_y = [0.0]
turbine_z = zeros(nturbines)
turbine_yaw = zeros(nturbines)
turbine_ct = zeros(nturbines) .+ ct
turbine_ai = zeros(nturbines) .+ ai
winddirections = [270.0*pi/180.0]
ambient_ti = ones(length(winddirections)) .* 0.046
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheight = [hub_height]
shearexponent = 0.15
turbine_inflow_velcities = zeros(nturbines) .+ wind_speed

# rotor sample points
rotor_points_y = [0.0]
rotor_points_z = [0.0]

ct_model = ff.ThrustModelConstantCt(ct)
power_model = ff.PowerModelConstantCp(cp)
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# turbine = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [ct_model], [power_model])
turbine = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], ct_model, power_model)
# turbine_definitions = [turbine for i in 1:nturbines]
turbine_definitions = [turbine]
sorted_turbine_index = [i for i  in 1:nturbines]
turbine_definition_ids = ones(Int, nturbines)


windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_ti, [wind_shear_model])

wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()

local_ti_model = ff.LocalTIModelNoLocalTI()
