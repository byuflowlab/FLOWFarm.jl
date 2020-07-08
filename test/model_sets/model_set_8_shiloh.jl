import FlowFarm; const ff = FlowFarm
using DelimitedFiles
using CSV 
using PyPlot

# based on IEA case study 3

# set initial turbine x and y locations
layout_file_name = "./inputfiles/shiloh_layout.csv"
layout_data = CSV.read(layout_file_name)

# println(layout_data)
# println(layout_data.elevation[10])

nturbines = length(layout_data.elevation)
lat = zeros(nturbines)
long = zeros(nturbines)
hub_height = zeros(nturbines)
model = layout_data.model
global m 
m = 0
    
for i in 1:nturbines
    lat[i] = layout_data.latitude[i]
    long[i] = layout_data.longitude[i]
    # println(layout_data.model[i])
    if ismissing(layout_data.elevation[i])
        hub_height[i] = (layout_data.elevation[i-1] + layout_data.elevation[i+1])/2.0 + layout_data.hub_height[i]
        global m = m + 1
    else
        hub_height[i] = layout_data.elevation[i] + layout_data.hub_height[i]
    end
end
utm_zone = 10
turbine_x, turbine_y = ff.latlong_to_xy(lat, long, utm_zone)
# turbine_x = [0.0]
# turbine_y = [0.0]

# calculate the number of turbines
nturbines = length(turbine_x)

# set turbine base heights
turbine_z = zeros(nturbines)

# set turbine yaw values
turbine_yaw = zeros(nturbines)

# load ge turb data
turbge1p5file = string("./inputfiles/ge15-77_thrust_power.csv")

gedata = CSV.read(turbge1p5file)

ngepoints = length(gedata.id)
gespeed = zeros(ngepoints)
gecp= zeros(ngepoints)
gethrust = zeros(ngepoints)
gerd = gedata.rd[1]
gemodel = gedata.model[1]
for i in 1:ngepoints
    gespeed[i] = gedata.wind_speed[i]
    gecp[i] = gedata.power[i]
    gethrust[i] = gedata.thrust[i]
end

# https://www.thewindpower.net/turbine_en_56_ge-energy_1.5sl.php
gecutin = 3.0
gecutout = 25.0
geratedspeed =  15.0
geratedpower =  1.5E6

# load mm turb data
turbmm92file = string("./inputfiles/mm92_thrust_power.csv")

mmdata = CSV.read(turbmm92file)

nmmpoints = length(mmdata.id)
mmspeed = zeros(nmmpoints)
mmcp = zeros(nmmpoints)
mmthrust = zeros(nmmpoints)
mmrd = mmdata.rd[1]
mmmodel = mmdata.model[1]
for i in 1:nmmpoints
    mmspeed[i] = mmdata.wind_speed[i]
    mmcp[i] = mmdata.power[i]
    mmthrust[i] = mmdata.thrust[i]
end

# https://www.thewindpower.net/turbine_en_15_repower_mm92.php
mmcutin = 3.0
mmcutout = 25.0
mmratedspeed =  10.5
mmratedpower =  2.0E6

rotor_diameter = zeros(nturbines) # m
cut_in_speed = zeros(nturbines)   # m/s
cut_out_speed = zeros(nturbines)   # m/s
rated_speed = zeros(nturbines)  # m/s
rated_power = zeros(nturbines)  # W
generator_efficiency = zeros(nturbines) .+ 1.0

for i in 1:nturbines

    if layout_data.model[i] == gemodel
        rotor_diameter[i] = gerd    # m
        cut_in_speed[i] = gecutin   # m/s
        cut_out_speed[i] = gecutout   # m/s
        rated_speed[i] = geratedspeed  # m/s
        rated_power[i] = geratedpower  # W
    elseif layout_data.model[i] == mmmodel
        rotor_diameter[i] = mmrd    # m
        cut_in_speed[i] = mmcutin   # m/s
        cut_out_speed[i] = mmcutout   # m/s
        rated_speed[i] = mmratedspeed  # m/s
        rated_power[i] = mmratedpower  # W
    end
end

# rotor swept area sample points (normalized by rotor radius)
rotor_points_y = [0.0]
rotor_points_z = [0.0]

fig, ax = plt.subplots()
plot([minimum(turbine_x),minimum(turbine_x),maximum(turbine_x), maximum(turbine_x), 
minimum(turbine_x)],[minimum(turbine_y),maximum(turbine_y),maximum(turbine_y),
minimum(turbine_y),minimum(turbine_y)])

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[i]/2.0, fill=false,color="C1", linestyle="--")) 
end

# println(nturbines)
# println(m)
# println(minimum(layout_data.hub_height))
# # set up and show plot
# axis("square")
plt.show()

# set flow parameters
wind_rose_file_name = string("./inputfiles/wind_rose_jepson_prarie.txt")
winddata = readdlm(wind_rose_file_name, ',', skipstart=9)
# println(winddata)
nstates = length(winddata[:,1])*6

speeds = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0]
winddirections = zeros(nstates)
windspeeds = zeros(nstates)
windprobabilities = zeros(nstates)
for i in 1:length(winddata[:,1])
    dirstr = split(winddata[i,1] , r"[-]")
    for j in 1:6
        winddirections[(i-1)*6 + j] = (0.5*(parse(Float64, dirstr[2])+parse(Float64,dirstr[1])))
        windspeeds[(i-1)*6 + j] = speeds[j]
        windprobabilities[(i-1)*6 + j] = winddata[i,j+2] / 100.0
    end
end

winddirections *= pi/180.0

air_density = 1.1716  # kg/m^3
shearexponent = 0.15
ambient_tis = zeros(nstates) .+ 0.1
measurementheight = zeros(nstates) .+ sum(hub_height)/nturbines

# initialize power model
gepower_model = ff.PowerModelCpPoints(gespeed, gecp)
mmpower_model = ff.PowerModelCpPoints(mmspeed, mmcp)
power_models = Vector{typeof(gepower_model)}(undef, nturbines)
for i = 1:nturbines
    if model[i] == gemodel
        power_models[i] = gepower_model
    elseif model[i] == mmmodel
        power_models[i] = mmpower_model
    else
        throw(UndefVarError(:model))
    end
end

# initialize thurst model
gect_model = ff.ThrustModelCtPoints(gespeed, gethrust)
mmct_model = ff.ThrustModelCtPoints(mmspeed, mmthrust)
ct_models = Vector{typeof(gect_model)}(undef, nturbines)
for i = 1:nturbines
    if model[i] == gemodel
        ct_models[i] = gect_model
    elseif model[i] == mmmodel
        ct_models[i] = mmct_model
    else
        throw(UndefVarError(:model))
    end
end

# initialize wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# get sorted indecies 
sorted_turbine_index = sortperm(turbine_x)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource([winddirections[1]], [windspeeds[1]], [windprobabilities[1]], [measurementheight[1]], air_density, [ambient_tis[1]], wind_shear_model)

# set up wake and related models
wakedeficitmodel = ff.GaussYawVariableSpread()

wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelMaxTI()

# initialize model set
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)