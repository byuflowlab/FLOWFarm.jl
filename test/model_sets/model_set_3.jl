import FlowFarm; const ff = FlowFarm

# based on https://backend.orbit.dtu.dk/ws/portalfiles/portal/146793996/Pena18_wes_3_191_2018.pdf

diam = 120.0

locdata = readdlm("inputfiles/AnholtOffshoreWindFarmLocations.txt", ' ', skipstart=1)
turbine_x = locdata[:, 1]*diam
turbine_y = locdata[:, 2]*diam
nturbines = length(turbine_x)

turbine_z = zeros(nturbines)
rotor_diameter = zeros(nturbines).+diam
hub_height = zeros(nturbines).+81.6
turbine_yaw = zeros(nturbines).+0.0
turbine_ai = zeros(nturbines).+(1.0/3.0)

generator_efficiency = zeros(nturbines).+0.944
cut_in_speed = zeros(nturbines).+3.5  # m/s
cut_out_speed = zeros(nturbines).+25.  # m/s
rated_speed = zeros(nturbines).+12.0  # m/s
rated_power = zeros(nturbines).+3.6E6  # W


# load cp data
cpdata = readdlm("inputfiles/power_curve_siemens_swp3.6-120.txt",  ',', skipstart=1)
vel_points = cpdata[:, 1]
power_points = cpdata[:, 2]

power_model = ff.PowerModelPowerPoints(vel_points, power_points)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# rotor sample points
rotor_points_y = [0.0]
rotor_points_z = [0.0]

winddirections = deg2rad.(range(0,stop=350,length=36))
windspeeds = zeros(36).+(9.23)
windprobabilities = [1.3,1.1,1.3,1.1,1.0,1.3,1.7,2.2,2.1,3.3,3.5,4.1,4.5,5.0,
                            3.6,3.3,2.8,2.6,2.8,3.7,4.7,4.3,3.6,3.5,3.1,3.6,4.0,3.6,
                            4.2,3.8,1.8,1.6,1.6,1.4,1.2,1.3]
windprobabilities = windprobabilities./(sum(windprobabilities))

air_density = 1.1716  # kg/m^3
ambient_ti = 0.1
measurementheight = copy(hub_height)
ambient_tis = ones(length(winddirections)).*ambient_ti
shearexponent = 0.15
wind_shear_model = ff.PowerLawWindShear(shearexponent)
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

alpha = 0.04
wakedeficitmodel = ff.JensenTopHat(alpha)
wakedeflectionmodel = ff.JiminezYawDeflection(alpha)
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

ct = 0.689
ct_model = ff.ThrustModelConstantCt(ct)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end
