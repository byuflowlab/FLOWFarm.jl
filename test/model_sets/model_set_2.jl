import FlowFarm; const ff = FlowFarm

rotor_diameter = 80.0
hub_height = 70.0
yaw = 0.0
ct = 0.689 
cp = 0.8

generator_efficiency = 0.944
ai = 1.0/3.0
wind_speed = 12.0
air_density = 1.1716  # kg/m^3
data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)
turbine_x = data[:, 1].*7.0*rotor_diameter #[0.0, 7.0*rotor_diameter, 6.0*rotor_diameter]
nturbines = length(turbine_x)
turbine_y = zeros(nturbines)
turbine_z = zeros(nturbines)
turbine_yaw = zeros(nturbines)
turbine_ct = zeros(nturbines) .+ ct
turbine_ai = zeros(nturbines) .+ ai
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheight = [hub_height]
shearexponent = 0.15
turbine_inflow_velcities = zeros(nturbines) .+ wind_speed

# rotor sample points 
rotor_points_y = [0.0]
rotor_points_z = [0.0]

ct_model = ff.ConstantCt(ct)
power_model = ff.ConstantCp([cp], [generator_efficiency])
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# turbine_definitions = Array{ff.AbstractTurbine, length(turbine_x)}
turbine_definitions = [] # or `Vector{Coords}(undef, x)`
turbine = ff.Turbine(1, [rotor_diameter], [hub_height], [ct_model], [power_model])
turbine_definitions = [turbine for i in 1:nturbines]
sorted_turbine_index = zeros(Int64, nturbines)
for i in 1:nturbines
    turbine = ff.Turbine(i, [rotor_diameter], [hub_height], [ct_model], [power_model])
    turbine_definitions[i] = turbine
    sorted_turbine_index[i] = Int64(i)
end 

windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definitions)
windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])

alpha = 0.1
wakedeficitmodel = ff.JensenTopHat(alpha)
horizontal_spread_rate = alpha
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()

ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)
pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])