
using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using PyPlot; const plt=PyPlot
using ForwardDiff


"""
Author: Jared J. Thomas
Date:   May 29, 2020

Test based on:
[1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu
"""
function plot_velocity_row_deficit()

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
end

function test_point_in_polygon()
    vertices = [0.0 0.0; 0.0 10.0; 10.0 10.0; 10.0 0.0] # square 
    # vertices = [0.0 0.0; 0.0 10.0; 10.0 10.0; 10.0 6.0; 5.0 6.0; 5.0 4.0; 10.0 4.0; 10.0 0.0] # square with cut out
    # vertices = [0.0 0.0; 0.0 10.0; 4.0 4.0; 10.0 2.0] # almost a triangle 
    # vertices = []

    normals = ff.boundary_normals_calculator(vertices)
    step = 5E-2
    pointx = -5:step:15.0
    pointy = -5:step:15.0
    
    # set up function for getting AD derivatives
    pointinpolygon_diff(x) = ff.pointinpolygon(x, vertices, normals, return_distance=true)

    distance = zeros((length(pointx), length(pointy)))
    xderiv = zeros((length(pointx), length(pointy)))
    yderiv = zeros((length(pointx), length(pointy)))

    for i = 1:length(pointx)
        for j = 1:length(pointy)
            point = [pointx[i], pointy[j]]
            distance[j,i] = ff.pointinpolygon(point, vertices, normals, return_distance=true)
            deriv = ForwardDiff.gradient(pointinpolygon_diff, point)
            xderiv[j,i] = deriv[1]
            yderiv[j,i] = deriv[2]
        end
    end

    fig, ax = plt.subplots(2,2)
    ff.plotboundary!(ax[1,1], vertices)
    maxmin = maximum(abs.(distance))
    cplot1 = ax[1,1].contourf(pointx, pointy, distance, cmap="PuOr", vmin=-maxmin, vmax=maxmin, levels=200)
    cbar1 = plt.colorbar(cplot1)

    ff.plotboundary!(ax[1,2], vertices)
    maxmin = maximum(abs.(xderiv))
    cplot2 = ax[1,2].contourf(pointx, pointy, xderiv, cmap="PuOr", vmin=-maxmin, vmax=maxmin, levels=200)
    cbar2 = plt.colorbar(cplot2)

    ff.plotboundary!(ax[2,1], vertices)
    maxmin = maximum(abs.(yderiv))
    cplot3 = ax[2,1].contourf(pointx, pointy, yderiv, cmap="PuOr", vmin=-maxmin, vmax=maxmin, levels=200)
    cbar3 = plt.colorbar(cplot3)
    
end