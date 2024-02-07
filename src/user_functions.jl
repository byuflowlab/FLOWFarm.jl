"""file with functions used in wind farm optimization to generalize the setup up function calls
created January 29, 2024
author: Benjamin Varela
"""

function build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
            ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
            power_models,model_set;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
            hours_per_year=365.25*24.0,objective_scale=0.0,boundary=nothing,
            spacing=nothing,update_function=x->x)

    n_turbines = length(turbine_x)

    ideal_AEP = calculate_ideal_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set;
                rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    if objective_scale == 0.0
        objective_scale = 1.0/ideal_AEP
    end

    n_threads = Threads.nthreads()

    input_type = eltype(turbine_x)

    preallocations = preallocations_struct(zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads))

    t = eltype(ForwardDiff.GradientConfig(nothing,x))

    preallocations_dual = preallocations_struct(zeros(t,n_turbines,n_threads),zeros(t,n_turbines,n_threads),
                zeros(t,n_turbines,n_threads),zeros(t,n_turbines,n_threads),zeros(t,n_turbines,n_turbines,n_threads),
                zeros(t,n_turbines,n_turbines,n_threads),zeros(t,n_turbines,n_turbines,n_threads),
                zeros(t,n_turbines,n_turbines,n_threads))

    turbine_x_dual = Vector{t}(turbine_x)
    turbine_y_dual = Vector{t}(turbine_y)
    hub_height_dual = Vector{t}(hub_height)
    turbine_yaw_dual = Vector{t}(turbine_yaw)

    aep_gradient = zeros(length(x))

    return wind_farm_struct(turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
                ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,
                wind_resource,power_models,model_set,rotor_sample_points_y,rotor_sample_points_z,
                hours_per_year,objective_scale,ideal_AEP,boundary,spacing,preallocations,n_turbines,
                turbine_x_dual,turbine_y_dual,hub_height_dual,turbine_yaw_dual,
                preallocations_dual,update_function,aep_gradient)
end

function update_turbine_x!(farm,x::Vector{T}) where T
    farm.turbine_x .= x
end

function update_turbine_x!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.turbine_x_dual .= x
end

function update_turbine_y!(farm,x::Vector{T}) where T
    farm.turbine_y .= x
end

function update_turbine_y!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.turbine_y_dual .= x
end

function update_hub_height!(farm,x::Vector{T}) where T
    farm.hub_height .= x
end

function update_hub_height!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.hub_height_dual .= x
end

function update_turbine_yaw!(farm,x::Vector{T}) where T
    farm.turbine_yaw .= x
end

function update_turbine_yaw!(farm,x::Vector{T}) where T <: ForwardDiff.Dual
    farm.turbine_yaw_dual .= x
end

function calculate_aep(x::Vector{T},farm) where T
    farm.update_function(farm,x)

    AEP = calculate_aep(farm.turbine_x, farm.turbine_y, farm.turbine_z, farm.rotor_diameter,
                farm.hub_height, farm.turbine_yaw, farm.ct_models, farm.generator_efficiency, farm.cut_in_speed,
                farm.cut_out_speed, farm.rated_speed, farm.rated_power, farm.wind_resource, farm.power_models, farm.model_set,
                rotor_sample_points_y=farm.rotor_sample_points_y,rotor_sample_points_z=farm.rotor_sample_points_z,
                prealloc_turbine_velocities=farm.preallocations.prealloc_turbine_velocities,
                prealloc_turbine_ct=farm.preallocations.prealloc_turbine_ct,
                prealloc_turbine_ai=farm.preallocations.prealloc_turbine_ai,
                prealloc_turbine_local_ti=farm.preallocations.prealloc_turbine_local_ti,
                prealloc_wake_deficits=farm.preallocations.prealloc_wake_deficits,
                prealloc_contribution_matrix=farm.preallocations.prealloc_contribution_matrix,
                prealloc_deflections=farm.preallocations.prealloc_deflections,
                prealloc_sigma_squared=farm.preallocations.prealloc_sigma_squared
                ) .* farm.objective_scale

    return AEP
end

function calculate_aep(x::Vector{T},farm) where T <: ForwardDiff.Dual
    farm.update_function(farm,x)

    AEP = calculate_aep(farm.turbine_x_dual, farm.turbine_y_dual, farm.turbine_z, farm.rotor_diameter,
                farm.hub_height_dual, farm.turbine_yaw_dual, farm.ct_models, farm.generator_efficiency, farm.cut_in_speed,
                farm.cut_out_speed, farm.rated_speed, farm.rated_power, farm.wind_resource, farm.power_models, farm.model_set,
                rotor_sample_points_y=farm.rotor_sample_points_y,rotor_sample_points_z=farm.rotor_sample_points_z,
                prealloc_turbine_velocities=farm.preallocations_dual.prealloc_turbine_velocities,
                prealloc_turbine_ct=farm.preallocations_dual.prealloc_turbine_ct,
                prealloc_turbine_ai=farm.preallocations_dual.prealloc_turbine_ai,
                prealloc_turbine_local_ti=farm.preallocations_dual.prealloc_turbine_local_ti,
                prealloc_wake_deficits=farm.preallocations_dual.prealloc_wake_deficits,
                prealloc_contribution_matrix=farm.preallocations_dual.prealloc_contribution_matrix,
                prealloc_deflections=farm.preallocations_dual.prealloc_deflections,
                prealloc_sigma_squared=farm.preallocations_dual.prealloc_sigma_squared
                ) .* farm.objective_scale

    return AEP
end
