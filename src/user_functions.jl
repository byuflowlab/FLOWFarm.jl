"""file with functions used in wind farm optimization to generalize the setup up function calls
created January 29, 2024
author: Benjamin Varela
"""

function build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
            ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
            power_models,model_set;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],AEP_scale=0.0,
            update_function=dv->dv,opt_x=false,opt_y=false,opt_height=false,opt_yaw=false,opt_diameter=false)

    n_turbines = length(turbine_x)

    ideal_AEP = calculate_ideal_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set;
                rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    if AEP_scale == 0.0
        AEP_scale = 1.0/ideal_AEP
    end

    n_threads = Threads.nthreads()

    input_type = eltype(turbine_x)

    preallocations = preallocations_struct(zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads))

    cfg = ForwardDiff.GradientConfig(nothing,x)
    cfg_jac = ForwardDiff.JacobianConfig(nothing,x)
    t = eltype(cfg)

    preallocations_dual = preallocations_struct(zeros(t,n_turbines,n_threads),zeros(t,n_turbines,n_threads),
                zeros(t,n_turbines,n_threads),zeros(t,n_turbines,n_threads),zeros(t,n_turbines,n_turbines,n_threads),
                zeros(t,n_turbines,n_turbines,n_threads),zeros(t,n_turbines,n_turbines,n_threads),
                zeros(t,n_turbines,n_turbines,n_threads))

    AEP_gradient = zeros(input_type,length(x))

    wind_farm_constants = wind_farm_constants_struct(turbine_z,ct_models,generator_efficiency,cut_in_speed,cut_out_speed,
                rated_speed,rated_power,wind_resource,power_models,model_set,rotor_sample_points_y,rotor_sample_points_z)

    turbine_x_dual = nothing
    turbine_y_dual = nothing
    hub_height_dual = nothing
    turbine_yaw_dual = nothing
    rotor_diameter_dual = nothing

    if opt_x
        turbine_x_dual = Vector{t}(turbine_x)
    end
    if opt_y
        turbine_y_dual = Vector{t}(turbine_y)
    end
    if opt_height
        hub_height_dual = Vector{t}(hub_height)
    end
    if opt_yaw
        turbine_yaw_dual = Vector{t}(turbine_yaw)
    end
    if opt_diameter
        rotor_diameter_dual = Vector{t}(rotor_diameter)
    end

    design_variable_duals = wind_farm_derivative_struct(turbine_x_dual,turbine_y_dual,hub_height_dual,turbine_yaw_dual,rotor_diameter_dual,cfg)

    return wind_farm_struct(turbine_x, turbine_y, hub_height, turbine_yaw, rotor_diameter, AEP_gradient,
                wind_farm_constants, AEP_scale, ideal_AEP, preallocations, preallocations_dual, update_function, design_variable_duals)
end

function build_spacing_struct(x,n_turbines,space,scale)
    n_constraints = n_turbines * (n_turbines - 1) ÷ 2
    spacing_vec = zeros(eltype(x),n_constraints)
    spacing_jacobian = zeros(eltype(x),n_constraints,length(x))
    cfg = ForwardDiff.JacobianConfig(nothing,spacing_vec,x)
    return spacing_struct(space,scale,spacing_vec,spacing_jacobian,false,cfg)
end

function build_boundary_struct(x,n_turbines,scaling,constraint_function)
    boundary_vec = zeros(eltype(x),n_turbines)
    boundary_jacobian = zeros(eltype(x),n_turbines,length(x))
    cfg = ForwardDiff.JacobianConfig(nothing,boundary_vec,x)
    return boundary_struct(scaling,constraint_function,boundary_vec,boundary_jacobian,false,cfg)
end

function update_turbine_x!(farm,x)
    farm.turbine_x .= x
end

function update_turbine_x!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.duals.turbine_x_dual .= x
end

function update_turbine_y!(farm,x)
    farm.turbine_y .= x
end

function update_turbine_y!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.duals.turbine_y_dual .= x
end

function update_hub_height!(farm,x)
    farm.hub_height .= x
end

function update_hub_height!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.duals.hub_height_dual .= x
end

function update_turbine_yaw!(farm,x)
    farm.turbine_yaw .= x
end

function update_turbine_yaw!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.duals.turbine_yaw_dual .= x
end

function update_rotor_diameter!(farm,x)
    farm.rotor_diameter .= x
end

function update_rotor_diameter!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.duals.rotor_diameter_dual .= x
end

function get_duals(farm)
    if isnothing(farm.duals.turbine_x_dual)
        x = farm.turbine_x
    else
        x = farm.duals.turbine_x_dual
    end

    if isnothing(farm.duals.turbine_y_dual)
        y = farm.turbine_y
    else
        y = farm.duals.turbine_y_dual
    end

    if isnothing(farm.duals.hub_height_dual)
        hub_height = farm.hub_height
    else
        hub_height = farm.duals.hub_height_dual
    end

    if isnothing(farm.duals.turbine_yaw_dual)
        yaw = farm.turbine_yaw
    else
        yaw = farm.duals.turbine_yaw_dual
    end

    if isnothing(farm.duals.rotor_diameter_dual)
        rotor_diameter = farm.rotor_diameter
    else
        rotor_diameter = farm.duals.rotor_diameter_dual
    end

    return x, y, hub_height, yaw, rotor_diameter
end

function calculate_aep(x,farm)
    farm.update_function(farm,x)

    AEP = calculate_aep(farm.turbine_x, farm.turbine_y, farm.constants.turbine_z, farm.rotor_diameter,
                farm.hub_height, farm.turbine_yaw, farm.constants.ct_models, farm.constants.generator_efficiency,
                farm.constants.cut_in_speed, farm.constants.cut_out_speed, farm.constants.rated_speed,
                farm.constants.rated_power, farm.constants.wind_resource, farm.constants.power_models,
                farm.constants.model_set,rotor_sample_points_y=farm.constants.rotor_sample_points_y,
                rotor_sample_points_z=farm.constants.rotor_sample_points_z,
                prealloc_turbine_velocities=farm.preallocations.prealloc_turbine_velocities,
                prealloc_turbine_ct=farm.preallocations.prealloc_turbine_ct,
                prealloc_turbine_ai=farm.preallocations.prealloc_turbine_ai,
                prealloc_turbine_local_ti=farm.preallocations.prealloc_turbine_local_ti,
                prealloc_wake_deficits=farm.preallocations.prealloc_wake_deficits,
                prealloc_contribution_matrix=farm.preallocations.prealloc_contribution_matrix,
                prealloc_deflections=farm.preallocations.prealloc_deflections,
                prealloc_sigma_squared=farm.preallocations.prealloc_sigma_squared
                ) .* farm.AEP_scale

    return AEP
end

function calculate_aep(x::Vector{T},farm) where T <: ForwardDiff.Dual
    farm.update_function(farm,x)

    dx, dy, dh, dyaw, d_diameter = get_duals(farm)

    AEP = calculate_aep(dx, dy, farm.constants.turbine_z, d_diameter, dh, dyaw, farm.constants.ct_models,
                farm.constants.generator_efficiency, farm.constants.cut_in_speed, farm.constants.cut_out_speed,
                farm.constants.rated_speed, farm.constants.rated_power, farm.constants.wind_resource,
                farm.constants.power_models, farm.constants.model_set,
                rotor_sample_points_y=farm.constants.rotor_sample_points_y,
                rotor_sample_points_z=farm.constants.rotor_sample_points_z,
                prealloc_turbine_velocities=farm.preallocations_dual.prealloc_turbine_velocities,
                prealloc_turbine_ct=farm.preallocations_dual.prealloc_turbine_ct,
                prealloc_turbine_ai=farm.preallocations_dual.prealloc_turbine_ai,
                prealloc_turbine_local_ti=farm.preallocations_dual.prealloc_turbine_local_ti,
                prealloc_wake_deficits=farm.preallocations_dual.prealloc_wake_deficits,
                prealloc_contribution_matrix=farm.preallocations_dual.prealloc_contribution_matrix,
                prealloc_deflections=farm.preallocations_dual.prealloc_deflections,
                prealloc_sigma_squared=farm.preallocations_dual.prealloc_sigma_squared
                ) .* farm.AEP_scale

    return AEP
end

function calculate_spacing!(spacing_vec,x,farm,spacing_struct)
    farm.update_function(farm,x)

    turbine_spacing!(spacing_vec,farm.turbine_x,farm.turbine_y)

    spacing_vec .= (spacing_struct.constraint_spacing .- spacing_vec) .* spacing_struct.constraint_scaling

    return spacing_vec
end

function calculate_spacing!(spacing_vec,x::Vector{T},farm,spacing_struct) where T <: ForwardDiff.Dual
    farm.update_function(farm,x)

    dx, dy, _, _, _ = get_duals(farm)

    turbine_spacing!(spacing_vec,dx,dy)

    spacing_vec .= (spacing_struct.constraint_spacing .- spacing_vec) .* spacing_struct.constraint_scaling

    return spacing_vec
end

function calculate_boundary!(boundary_vec,x,farm,boundary_struct)
    farm.update_function(farm,x)

    boundary_struct.boundary_function(boundary_vec,farm.turbine_x,farm.turbine_y)

    boundary_vec .*= boundary_struct.boundary_scaling_factor

    return boundary_vec
end

function calculate_boundary!(boundary_vec,x::Vector{T},farm,boundary_struct) where T <: ForwardDiff.Dual
    farm.update_function(farm,x)

    dx, dy, _, _, _ = get_duals(farm)

    boundary_struct.boundary_function(boundary_vec,dx,dy)

    boundary_vec .*= boundary_struct.boundary_scaling_factor

    return boundary_vec
end
