@everywhere import FlowFarm; const ff = FlowFarm

# based on IEA case study 4

# set initial turbine x and y locations
layout_file_name = "./inputfiles/iea37-ex-opt4.yaml"
turbine_x, turbine_y, fname_turb, fname_wr = ff.get_turb_loc_YAML(layout_file_name)

# calculate the number of turbines
nturbines = length(turbine_x)

# set turbine base heights
turbine_z = zeros(nturbines)

# set turbine yaw values
turbine_yaw = zeros(nturbines)

# set turbine design parameters
turbine_file_name = string("./inputfiles/",fname_turb)
turb_ci, turb_co, rated_ws, rated_pwr, turb_diam, turb_hub_height = ff.get_turb_atrbt_YAML(turbine_file_name)

rotor_diameter = zeros(nturbines) .+ turb_diam # m
hub_height = zeros(nturbines) .+ turb_hub_height   # m
cut_in_speed = zeros(nturbines) .+ turb_ci  # m/s
cut_out_speed = zeros(nturbines) .+ turb_co  # m/s
rated_speed = zeros(nturbines) .+ rated_ws # m/s
rated_power = zeros(nturbines) .+ rated_pwr # W
generator_efficiency = zeros(nturbines) .+ 1.0

# rotor swept area sample points (normalized by rotor radius)
rotor_points_y = [0.0]
rotor_points_z = [0.0]

# set flow parameters
windrose_file_name = string("./inputfiles/",fname_wr)
winddirections, windspeeds, windprobabilities, ambient_ti = ff.get_reduced_wind_rose_YAML(windrose_file_name)

# winddirections = [winddirections[1]]
# windspeeds = [windspeeds[1]]
# windprobabilities = [windprobabilities[1]]

nstates = length(winddirections)
winddirections *= pi/180.0

air_density = 1.1716  # kg/m^3
shearexponent = 0.15
ambient_tis = zeros(nstates) .+ ambient_ti
measurementheight = zeros(nstates) .+ turb_hub_height

# initialize power model
power_model = ff.PowerModelPowerCurveCubic()
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# load thrust curve
ct = 4.0*(1.0/3.0)*(1.0 - 1.0/3.0)

# initialize thurst model
ct_model = ff.ThrustModelConstantCt(ct)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end

# initialize wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# get sorted indecies 
sorted_turbine_index = sortperm(turbine_x)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

# set up wake and related models
k = 0.0324555
wakedeficitmodel = ff.GaussSimple(k)

wakedeflectionmodel = ff.JiminezYawDeflection()
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

# initialize model set
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)