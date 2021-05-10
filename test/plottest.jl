"""
Author: Jared J. Thomas
Date:   May 29, 2020

Test based on:
[1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu
"""

using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using PyPlot

# load data
data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)

# load problem set up
include("./model_sets/model_set_2.jl")

# calculate wind turbine velocities and corresponding aerodynamic operational states
turbine_inflow_velcities, turbine_ct, turbine_ai, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, windresource,
model_set, wind_farm_state_id=1)

# calculate the power production of each wind turbine
ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, air_density, power_model)

# set up the point locations to be used in generating the plot
stepsize = 5
xrange = 1:stepsize:10*rotor_diameter[1]*nturbines
yrange = -1*rotor_diameter[1]*nturbines:stepsize:1*rotor_diameter[1]*nturbines
zrange = hub_height[1]:stepsize:hub_height[1]

# calculate the wind speed at each point in the flow field 
velh = ff.calculate_flow_field(xrange, yrange, zrange,
model_set, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_inflow_velcities,
windresource; wind_farm_state_id=1)

# visualize the resulting data
flowfieldplot = contourf(xrange, yrange, velh, cmap="Blues_r")