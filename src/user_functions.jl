"""file with functions used in wind farm optimization to generalize the setup up function calls
created January 29, 2024
author: Benjamin Varela
"""

function build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
            ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
            power_models,model_set;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],AEP_scale=0.0,
            update_function=dv->dv,input_type=nothing,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false)

    n_turbines = length(turbine_x)

    ideal_AEP = calculate_ideal_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set;
                rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    if AEP_scale == 0.0
        AEP_scale = 1.0/ideal_AEP
    end

    n_threads = Threads.nthreads()
    cfg = ForwardDiff.GradientConfig(nothing,x)

    if input_type === nothing
        input_type = eltype(x)
    elseif input_type == "ForwardDiff"
        input_type = eltype(cfg)
    end

    preallocations = preallocations_struct(zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads),
                zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads))

    results = DiffResults.GradientResult(x)

    AEP_gradient = spzeros(Float64,length(x))
    AEP = Array{Float64,0}(undef)

    wind_farm_constants = wind_farm_constants_struct(turbine_z,ct_models,generator_efficiency,cut_in_speed,cut_out_speed,
                rated_speed,rated_power,wind_resource,power_models,model_set,rotor_sample_points_y,rotor_sample_points_z)

    if opt_x
        turbine_x = Vector{input_type}(turbine_x)
    end
    if opt_y
        turbine_y = Vector{input_type}(turbine_y)
    end
    if opt_hub
        hub_height = Vector{input_type}(hub_height)
    end
    if opt_yaw
        turbine_yaw = Vector{input_type}(turbine_yaw)
    end
    if opt_diam
        rotor_diameter = Vector{input_type}(rotor_diameter)
    end

    return wind_farm_struct(turbine_x, turbine_y, hub_height, turbine_yaw, rotor_diameter, results,
                wind_farm_constants, AEP_scale, ideal_AEP, preallocations, update_function, AEP_gradient, AEP, cfg)
end

function build_spacing_struct(x,n_turbines,space,scale,update_function)
    n_constraints = n_turbines * (n_turbines - 1) ÷ 2
    spacing_vec = zeros(eltype(x),n_constraints)
    spacing_jacobian = spzeros(eltype(x),n_constraints,length(x))
    cfg = ForwardDiff.JacobianConfig(nothing,spacing_vec,x)
    turbine_x = zeros(eltype(cfg),n_turbines)
    turbine_y = zeros(eltype(cfg),n_turbines)
    return spacing_struct(turbine_x,turbine_y,space,scale,spacing_vec,spacing_jacobian,cfg,update_function)
end

function build_boundary_struct(x,n_turbines,n_constraints,scaling,constraint_function,update_function;using_sparsity=true)
    boundary_vec = zeros(Float64,n_constraints)
    boundary_jacobian = spzeros(Float64,n_constraints,length(x))

    cfg = ForwardDiff.JacobianConfig(nothing,boundary_vec,x)
    T = eltype(cfg)
    turbine_x = zeros(T,n_turbines)
    turbine_y = zeros(T,n_turbines)
    b_struct = boundary_struct(turbine_x,turbine_y,scaling,constraint_function,boundary_vec,
    boundary_jacobian,cfg,update_function)
    if !using_sparsity
        return b_struct
    end

    calculate_boundary(a,b) = calculate_boundary!(a,b,b_struct)
    x_temp = copy(x)
    for i = 1:3
        x_temp .+= rand(length(x_temp))
        ForwardDiff.jacobian!(b_struct.jacobian,calculate_boundary,b_struct.boundary_vec,x_temp,b_struct.config)
        boundary_jacobian .+= b_struct.jacobian
    end
    boundary_jacobian .= dropzeros(boundary_jacobian)
    ad = AutoSparseForwardDiff()
    sd = JacPrototypeSparsityDetection(; jac_prototype=boundary_jacobian)
    cache = sparse_jacobian_cache(ad, sd, nothing, boundary_vec, x)
    T = eltype(cache.cache.t)
    turbine_x = zeros(T,n_turbines)
    turbine_y = zeros(T,n_turbines)
    return sparse_boundary_struct(turbine_x,turbine_y,boundary_jacobian,ad,cache,boundary_vec,constraint_function,update_function,scaling)
end

function calculate_aep!(farm,x)
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

function calculate_aep_gradient!(farm,x)
    calculate_aep(x) = calculate_aep!(farm,x)
    ForwardDiff.gradient!(farm.results,calculate_aep,x,farm.config)
    farm.AEP .= DiffResults.value(farm.results)
    farm.AEP_gradient .= DiffResults.gradient(farm.results)
    return farm.AEP[1], farm.AEP_gradient
end

function calculate_spacing!(spacing_vec,x,spacing_struct)
    spacing_struct.update_function(spacing_struct,x)

    turbine_spacing!(spacing_vec,spacing_struct.turbine_x,spacing_struct.turbine_y)

    spacing_vec .= (spacing_struct.constraint_spacing .- spacing_vec) .* spacing_struct.constraint_scaling

    return spacing_vec
end

function calculate_spacing_jacobian!(spacing_struct,x)
    calculate_spacing(a,b) = calculate_spacing!(a,b,spacing_struct)
    ForwardDiff.jacobian!(spacing_struct.jacobian,calculate_spacing,spacing_struct.spacing,x,spacing_struct.config)
    spacing_struct.jacobian .= dropzeros(spacing_struct.jacobian)
    return spacing_struct.spacing, spacing_struct.jacobian
end

function calculate_boundary!(boundary_vec,x,boundary_struct)
    boundary_struct.update_function(boundary_struct,x)

    boundary_struct.boundary_function(boundary_vec,boundary_struct.turbine_x,boundary_struct.turbine_y)

    boundary_vec .*= boundary_struct.boundary_scaling_factor

    return boundary_vec
end

function calculate_boundary_jacobian!(boundary_struct,x)
    calculate_boundary(a,b) = calculate_boundary!(a,b,boundary_struct)
    ForwardDiff.jacobian!(boundary_struct.jacobian,calculate_boundary,boundary_struct.boundary_vec,x,boundary_struct.config)
    boundary_struct.jacobian .= dropzeros(boundary_struct.jacobian)
    return boundary_struct.boundary_vec, boundary_struct.jacobian
end

function calculate_boundary_jacobian!(boundary_struct::T,x) where T <: AbstractSparseMethod
    calculate_boundary(a,b) = calculate_boundary!(a,b,boundary_struct)
    sparse_jacobian!(boundary_struct.jacobian, boundary_struct.ad, boundary_struct.cache, calculate_boundary, boundary_struct.boundary_vec, x)
    return boundary_struct.boundary_vec, boundary_struct.jacobian
end
