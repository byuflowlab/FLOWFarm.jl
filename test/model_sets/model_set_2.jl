import FLOWFarm; const ff = FLOWFarm



diam = 80.0
data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)
turbine_x = data[:, 1].*7.0*diam
nturbines = length(turbine_x)

turbine_y = zeros(nturbines)
turbine_z = zeros(nturbines)

ambient_ti = 0.1
ct = 0.689
ai = 1.0/3.0
turbine_ct = zeros(nturbines) .+ ct
turbine_ai = zeros(nturbines) .+ ai
turbine_yaw = zeros(nturbines)
sorted_turbine_index = [i for i  in 1:nturbines]

rotor_sample_points_y = [0.0]
rotor_sample_points_z = [0.0]

rotor_diameter = zeros(nturbines).+diam
hub_height = zeros(nturbines).+70.0
generator_efficiency = zeros(nturbines).+0.944
cut_in_speed = zeros(nturbines).+4.  # m/s
cut_out_speed = zeros(nturbines).+25.  # m/s
rated_speed = zeros(nturbines).+16.  # m/s
rated_power = zeros(nturbines).+2.0E6  # W
air_density = 1.1716  # kg/m^3

constcp = 0.8
power_model = ff.PowerModelConstantCp(constcp)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

ct_model1 = ff.ThrustModelConstantCt(ct)
ct_model = Vector{typeof(ct_model1)}(undef, nturbines)
for i = 1:nturbines
    ct_model[i] = ct_model1
end

wind_speed = 12.0
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheight = [hub_height[1]]
ambient_tis = [ambient_ti]
shearexponent = 0.15
wind_shear_model = ff.PowerLawWindShear(shearexponent)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

alpha = 0.1
wakedeficitmodel = ff.JensenTopHat(alpha)
horizontal_spread_rate = alpha
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
