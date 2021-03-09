import FlowFarm; const ff = FlowFarm
using DelimitedFiles

# Defining Turbine Parameters and locations
D = 126.4
locx = 7D
locy = 0
locz = 90
turbine_x = [7D,0]
turbine_y = [0,0]
turbine_z = [0,0]
nturbines = length(turbine_x)

upstream_turbine_id = 2
downstream_turbine_id = 1
hub_height = [90.0,90.0]
rotor_diameter = [D,D]
turbine_ai = [1.0/3.0,1.0/3.0]
turbine_yaw = [0, 0]
dt = rotor_diameter[upstream_turbine_id]
cut_in_speed = [0.0, 0.0]
cut_out_speed = [25.0, 25.0]
rated_speed = [12.0, 12.0]
rated_power = [1.0176371581904552e6 , 1.0176371581904552e6]
#            turbine_yaw = [-40,-20,0,20,40]

# define power parameters
u = 8.0
air_density = 1.225
generator_efficiency = [.944, .944]
rotor_area = .25*pi*D^2
Fd = 0
ai = 1.0/3.0
ambient_ti = .06


turbine_local_ti = [ambient_ti,ambient_ti]

# values defining wake interaction from Gebraad et. al
me = [-.5,.22,1]
ke = .065
MU = [.5,1,5.5]
aU = 5
bU = 1.66

# identify wake models
wakedeficitmodel = ff.Multizone()
wakedeflectionmodel = ff.MultizoneDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelMaxTI()
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

# define wind characteristics
winddirections = [0.0]
windspeeds = [u]
windprobabilities = [1.0]
measurementheights = [hub_height[1]]
ambient_tis = [ambient_ti]
shearexponent = 0.15
wind_shear_model = ff.PowerLawWindShear(shearexponent)
wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

# Sample wind point at the hub
rotor_sample_points_y = [0]
rotor_sample_points_z = [0]
sorted_turbine_index = sortperm(turbine_x)

eta = .768

# import cp and ct corresponding to wind velocity
<<<<<<< HEAD
powerdata = readdlm(".inputfiles/NREL5MWCPCT.txt", skipstart=1)
=======
powerdata = readdlm("./inputfiles/NREL5MWCPCT.txt", skipstart=1)
>>>>>>> b7d53ceaffd9ae50fb772b96cc0fab2dc2541096
vel_points = powerdata[:,1]
cp_points = powerdata[:,2]
ct_points = powerdata[:,3]

# initialize power model
power_model = ff.PowerModelCpPoints(vel_points, cp_points)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# initialize thurst model
ct_model = ff.ThrustModelCtPoints(vel_points, ct_points)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end