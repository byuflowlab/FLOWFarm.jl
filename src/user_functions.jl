export build_wind_farm_struct, build_spacing_struct, build_boundary_struct, calculate_aep!, calculate_aep_gradient!,
    calculate_spacing!, calculate_spacing_jacobian!, calculate_boundary!, calculate_boundary_jacobian!
    
"""file with functions used in wind farm optimization to generalize the setup up of function calls
created January 29, 2024
author: Benjamin Varela
"""

"""
build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
            ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
            power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
            AEP_scale=0.0,input_type=nothing,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,
            force_single_thread=false)

function to build a wind_farm_struct

# Arguments
- `x`: Vector containing the desired design variables for optimization (if any)
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `turbine_z`: Vector containing z positions of turbines (this is where the base meets the ground)
- `hub_height`: Vector containing hub heights of turbines
- `turbine_yaw`: Vector containing yaw angles of turbines
- `rotor_diameter`: Vector containing rotor diameters of turbines
- `ct_models`: Vector containing ct_models for each turbine
- `generator_efficiency`: Vector containing generator efficiencies for each turbine
- `cut_in_speed`: Vector containing cut in speeds for each turbine
- `cut_out_speed`: Vector containing cut out speeds for each turbine
- `rated_speed`: Vector containing rated speeds for each turbine
- `rated_power`: Vector containing rated powers for each turbine
- `wind_resource`: The DiscretizedWindResource struct
- `power_models`: Vector containing power models for each turbine
- `model_set`: The WindFarmModelSet for the wind farm
- `update_function`: Function that updates the wind farm struct with the new design variables
- `rotor_sample_points_y`: Vector containing the y positions of the rotor sample points
- `rotor_sample_points_z`: Vector containing the z positions of the rotor sample points
- `AEP_scale`: Single float that scales the AEP, if 0.0 will be set to 1.0/ideal_AEP
- `input_type`: default is nothing and will be set to the type of x, if "ForwardDiff" then the input type will be set to ForwardDiff.dual
- `opt_x`: Boolean to optimize x positions of turbines
- `opt_y`: Boolean to optimize y positions of turbines
- `opt_hub`: Boolean to optimize hub heights of turbines
- `opt_yaw`: Boolean to optimize yaw angles of turbines
- `opt_diam`: Boolean to optimize rotor diameters of turbines
- `force_single_thread`: Boolean to force single thread calculation
"""
function build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
            ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
            power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
            AEP_scale=0.0,input_type=nothing,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,
            force_single_thread=false)

    n_turbines = length(turbine_x)
    n_threads = Threads.nthreads()
    force_single_thread && (n_threads = 1)
    results = DiffResults.GradientResult(x)
    AEP_gradient = zeros(eltype(x),length(x))
    AEP = Array{eltype(x),0}(undef)

    wind_farm_constants = wind_farm_constants_struct(turbine_z,ct_models,generator_efficiency,cut_in_speed,cut_out_speed,
                rated_speed,rated_power,wind_resource,power_models,model_set,rotor_sample_points_y,rotor_sample_points_z)

    ideal_AEP = calculate_ideal_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set;
                rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    (AEP_scale == 0.0) && (AEP_scale = 1.0/ideal_AEP)

    cfg = nothing
    if input_type === nothing
        cfg = nothing
        input_type = eltype(x)
    elseif input_type == "ForwardDiff"
        cfg = ForwardDiff.GradientConfig(nothing,x)
        input_type = eltype(cfg)
    elseif input_type == "ForwardDiffJacobian"
        y = eltype(x).(collect(1:n_turbines))
        cfg = ForwardDiff.JacobianConfig(nothing,y,x)
        input_type = eltype(cfg)
    end

    opt_x && (turbine_x = Vector{input_type}(turbine_x))
    opt_y && (turbine_y = Vector{input_type}(turbine_y))
    opt_hub && (hub_height = Vector{input_type}(hub_height))
    opt_yaw && (turbine_yaw = Vector{input_type}(turbine_yaw))
    opt_diam && (rotor_diameter = Vector{input_type}(rotor_diameter))

    preallocations = preallocations_struct(zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),
                    zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,n_threads),zeros(input_type,n_turbines,
                    n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads),
                    zeros(input_type,n_turbines,n_turbines,n_threads),zeros(input_type,n_turbines,n_turbines,n_threads))

    return wind_farm_struct(turbine_x, turbine_y, hub_height, turbine_yaw, rotor_diameter, results,
                wind_farm_constants, AEP_scale, ideal_AEP, preallocations, update_function, AEP_gradient, AEP, cfg, force_single_thread)
end

"""
build_spacing_struct(x,n_turbines,space,scale,update_function)

function to build a spacing_struct

# Arguments
- `x`: Vector containing the desired design variables for optimization
- `n_turbines`: Number of turbines in the wind farm
- `space`: The minimum spacing between turbines
- `scale`: The scaling factor for the spacing constraint
- `update_function`: Function that updates the spacing struct with the new design variables
"""
function build_spacing_struct(x,n_turbines,space,scale,update_function)
    n_constraints = n_turbines * (n_turbines - 1) ÷ 2
    spacing_vec = zeros(eltype(x),n_constraints)
    spacing_jacobian = zeros(eltype(x),n_constraints,length(x))
    cfg = ForwardDiff.JacobianConfig(nothing,spacing_vec,x)
    turbine_x = zeros(eltype(cfg),n_turbines)
    turbine_y = zeros(eltype(cfg),n_turbines)
    return spacing_struct(turbine_x,turbine_y,space,scale,spacing_vec,spacing_jacobian,cfg,update_function)
end

"""
build_boundary_struct

build_boundary_struct(x,n_turbines,n_constraints,scaling,constraint_function,update_function;using_sparsity=true)

# Arguments
- `x`: Vector containing the desired design variables for optimization
- `n_turbines`: Number of turbines in the wind farm
- `n_constraints`: Number of boundary constraints (n_turbines * number of sides of the boundary)
- `scaling`: The scaling factor for the boundary constraint
- `constraint_function`: Function that calculates the boundary constraints
- `update_function`: Function that updates the boundary struct with the new design variables
- `using_sparsity`: Boolean to use sparsity in the jacobian calculation (default is true)
"""
function build_boundary_struct(x,n_turbines,n_constraints,scaling,constraint_function,update_function;using_sparsity=true)
    boundary_vec = zeros(eltype(x),n_constraints)
    boundary_jacobian = zeros(eltype(x),n_constraints,length(x))

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
    x_temp = similar(x)
    for i = 1:3
        x_temp .= x + rand(length(x_temp)) * 1E-5
        ForwardDiff.jacobian!(b_struct.jacobian,calculate_boundary,b_struct.boundary_vec,x_temp,b_struct.config)
        boundary_jacobian .+= b_struct.jacobian
    end
    boundary_jacobian = dropzeros(sparse(boundary_jacobian))
    ad = AutoSparseForwardDiff()
    sd = JacPrototypeSparsityDetection(; jac_prototype=boundary_jacobian)
    cache = sparse_jacobian_cache(ad, sd, nothing, boundary_vec, x)
    T = eltype(cache.cache.t)
    turbine_x = zeros(T,n_turbines)
    turbine_y = zeros(T,n_turbines)
    return sparse_boundary_struct(turbine_x,turbine_y,boundary_jacobian,ad,cache,boundary_vec,constraint_function,update_function,scaling)
end

"""
calculate_aep!

function calculate_aep!(farm,x)

# Arguments
- `farm`: The wind_farm_struct
- `x`: Vector containing the design variables
"""
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
                prealloc_sigma_squared=farm.preallocations.prealloc_sigma_squared,
                force_single_thread=farm.force_single_thread
                ) .* farm.AEP_scale

    return AEP
end

"""
calculate_aep_gradient!(farm,x)

function to calculate the AEP and its gradient for the wind farm, results are stored within the wind_farm_struct and take into account the scaling factor

# Arguments
- `farm`: The wind_farm_struct
- `x`: Vector containing the design variables
"""
function calculate_aep_gradient!(farm,x)
    calculate_aep(x) = calculate_aep!(farm,x)
    ForwardDiff.gradient!(farm.results,calculate_aep,x,farm.config)
    farm.AEP .= DiffResults.value(farm.results)
    farm.AEP_gradient .= DiffResults.gradient(farm.results)
    return farm.AEP[1], farm.AEP_gradient
end

"""
calculate_aep_gradient!(farm,x,sparse_struct::T)

function to calculate the AEP and its gradient for the wind farm using sparse methods

# Arguments
- `farm`: The wind_farm_struct
- `x`: Vector containing the design variables
"""
function calculate_aep_gradient!(farm,x,sparse_struct::T) where T <: AbstractSparseMethod
    calculate_aep_gradient!(farm,x,sparse_struct)
    return farm.AEP[1], farm.AEP_gradient
end

"""
calculate_spacing!(spacing_vec,x,spacing_struct)

function to calculate the spacing constraints for the wind farm

# Arguments
- `spacing_vec`: Vector containing the spacing constraints (just us spacing_struct.spacing_vec)
- `x`: Vector containing the design variables
- `spacing_struct`: The spacing_struct
"""
function calculate_spacing!(spacing_vec,x,spacing_struct)
    spacing_struct.update_function(spacing_struct,x)

    turbine_spacing!(spacing_vec,spacing_struct.turbine_x,spacing_struct.turbine_y)

    spacing_vec .= (spacing_struct.constraint_spacing .- spacing_vec) .* spacing_struct.constraint_scaling

    return spacing_vec
end

"""
calculate_spacing_jacobian!(spacing_struct,x)

function to calculate the spacing constraints and the jacobian for the wind farm, results stored in spacing_struct

# Arguments
- `spacing_struct`: The spacing_struct
- `x`: Vector containing the design variables
"""
function calculate_spacing_jacobian!(spacing_struct,x)
    calculate_spacing(a,b) = calculate_spacing!(a,b,spacing_struct)
    ForwardDiff.jacobian!(spacing_struct.jacobian,calculate_spacing,spacing_struct.spacing_vec,x,spacing_struct.config)
    return spacing_struct.spacing_vec, spacing_struct.jacobian
end

"""
calculate_boundary!(boundary_vec,x,boundary_struct)

function to calculate the boundary constraints for the wind farm

# Arguments
- `boundary_vec`: Vector containing the boundary constraints (just us boundary_struct.boundary_vec)
- `x`: Vector containing the design variables
- `boundary_struct`: The boundary_struct
"""
function calculate_boundary!(boundary_vec,x,boundary_struct)
    boundary_struct.update_function(boundary_struct,x)

    boundary_struct.boundary_function(boundary_vec,boundary_struct.turbine_x,boundary_struct.turbine_y)

    boundary_vec .*= boundary_struct.boundary_scaling_factor

    return boundary_vec
end

"""
calculate_boundary_jacobian!(boundary_struct,x)

function to calculate the boundary constraints and the jacobian for the wind farm, results stored in boundary_struct

# Arguments
- `boundary_struct`: The boundary_struct
- `x`: Vector containing the design variables
"""
function calculate_boundary_jacobian!(boundary_struct,x)
    calculate_boundary(a,b) = calculate_boundary!(a,b,boundary_struct)
    ForwardDiff.jacobian!(boundary_struct.jacobian,calculate_boundary,boundary_struct.boundary_vec,x,boundary_struct.config)
    return boundary_struct.boundary_vec, boundary_struct.jacobian
end
