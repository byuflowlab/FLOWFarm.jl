"""file with functions used in wind farm optimization employing sparse methods
created January 26, 2024
author: Benjamin Varela
"""

abstract type AbstractSparseMethod end
abstract type StableSparseMethod <: AbstractSparseMethod end
abstract type UnstableSparseMethod <: AbstractSparseMethod end

#TODO: documentation and examples

# ∇AEP Optimization (stable pattern) #######################################################
struct sparse_AEP_struct_stable_pattern{T1,T2,T3,T4,T5} <: StableSparseMethod
    caches::T1 # vector of caches
    jacobians::T2 # vector of sparse jacobians
    state_gradients::T3 # 2d array, each row is a state gradient (used for threads)
    turbine_powers::T4 # 2d array that holds the powers or each turbine (used for threads)
    adtype::T5
end

function build_stable_sparse_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
    ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
    power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
    AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16)

    farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,input_type="ForwardDiff")

    sparse_struct = build_stable_sparse_struct(x,farm;tolerance=tolerance)

    farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type=eltype(sparse_struct.caches[1].cache.t))

    return farm,sparse_struct
end

function build_stable_sparse_struct(x,farm;tolerance=1E-16)
    n_states = length(farm.constants.wind_resource.wind_probabilities)
    n_turbines = length(farm.turbine_x)
    pow = zeros(n_turbines,n_states)
    jacobians = Array{SparseMatrixCSC{Float64, Int64},1}(undef,n_states)
    state_gradients = zeros(Float64,n_states,length(x))
    caches = nothing
    adtype = AutoSparseForwardDiff()
    for i = 1:n_states
        jacobians[i] = define_stable_jacobian_pattern(x,farm,tolerance,pow[:,i],i)
        sd = JacPrototypeSparsityDetection(; jac_prototype=jacobians[i])
        cache = sparse_jacobian_cache(adtype, sd, nothing, pow[:,i], x)
        if isnothing(caches)
            T = typeof(cache)
            adtype = AutoSparseForwardDiff(chunksize=cache.cache.chunksize)
            caches = Vector{T}(undef,n_states)
        end
        caches[i] = cache
    end

    return sparse_AEP_struct_stable_pattern(caches,jacobians,state_gradients,pow,adtype)
end

function define_stable_jacobian_pattern(x,farm,tolerance,pow,state_id)
    p(a,x) = calculate_wind_state_power!(a,x,farm,state_id)
    x_temp = similar(x)
    cfg = ForwardDiff.JacobianConfig(nothing,pow,x)
    jacobian = ForwardDiff.jacobian(p,pow,x,cfg)
    for i = 1:3
        x_temp .= x .+ (rand(size(x_temp)) .- 0.5) * 1E-4
        jacobian .+= ForwardDiff.jacobian(p,pow,x,cfg)
        jacobian[abs.(jacobian) .<= tolerance] .= 0.0
    end
    return dropzeros(sparse(jacobian))
end

function calculate_wind_state_power!(pow,x,farm,state_id;prealloc_id=1,hours_per_year=365.25*24.0)
    farm.update_function(farm,x)

    rot_x, rot_y = rotate_to_wind_direction(farm.turbine_x, farm.turbine_y,
                    farm.constants.wind_resource.wind_directions[state_id])
    sorted_turbine_index = sortperm(rot_x)
    turbine_velocities = ff.turbine_velocities_one_direction(rot_x, rot_y, farm.constants.turbine_z,
                    farm.rotor_diameter, farm.hub_height, farm.turbine_yaw, sorted_turbine_index,
                    farm.constants.ct_models, farm.constants.rotor_sample_points_y, farm.constants.rotor_sample_points_z,
                    farm.constants.wind_resource, farm.constants.model_set, wind_farm_state_id=state_id,
                    velocity_only=true,
                    turbine_velocities=view(farm.preallocations.prealloc_turbine_velocities,:,prealloc_id),
                    turbine_ct=view(farm.preallocations.prealloc_turbine_ct,:,prealloc_id),
                    turbine_ai=view(farm.preallocations.prealloc_turbine_ai,:,prealloc_id),
                    turbine_local_ti=view(farm.preallocations.prealloc_turbine_local_ti,:,prealloc_id),
                    wake_deficits=view(farm.preallocations.prealloc_wake_deficits,:,:,prealloc_id),
                    contribution_matrix=view(farm.preallocations.prealloc_contribution_matrix,:,:,prealloc_id),
                    deflections=view(farm.preallocations.prealloc_deflections,:,:,prealloc_id),
                    sigma_squared=view(farm.preallocations.prealloc_sigma_squared,:,:,prealloc_id))

    pow .= turbine_powers_one_direction(farm.constants.generator_efficiency, farm.constants.cut_in_speed,
                    farm.constants.cut_out_speed, farm.constants.rated_speed, farm.constants.rated_power,
                    farm.rotor_diameter, turbine_velocities, farm.turbine_yaw, farm.constants.wind_resource.air_density,
                    farm.constants.power_models)

    pow .= pow .* hours_per_year .* farm.constants.wind_resource.wind_probabilities[state_id] .* farm.AEP_scale
end

function calculate_aep_gradient!(farm,x,sparse_struct::T) where T <: StableSparseMethod
    n_threads = Threads.nthreads()
    n_states = length(farm.constants.wind_resource.wind_probabilities)

    if n_threads > 1
        n_per_thread, rem = divrem(n_states,n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:n_states

        Threads.@threads for i_assignment in eachindex(assignments)
            i_start = assignments[i_assignment]
            i_stop = min(i_start+n_per_thread-1, n_states)
            for i = i_start:i_stop
                p(a,x) = calculate_wind_state_power!(a,x,farm,i;prealloc_id=i_assignment)
                SparseDiffTools.sparse_jacobian!(sparse_struct.jacobians[i],sparse_struct.adtype,
                                sparse_struct.caches[i],p,view(sparse_struct.turbine_powers,:,i),x)
                sparse_struct.state_gradients[i,:] .= sum(sparse_struct.jacobians[i],dims=1)[:]
            end
        end
    else
        for i = 1:n_states
            p(a,x) = calculate_wind_state_power!(a,x,farm,i;prealloc_id=1)
            SparseDiffTools.sparse_jacobian!(sparse_struct.jacobians[i],sparse_struct.adtype,
                            sparse_struct.caches[i],p,view(sparse_struct.turbine_powers,:,i),x)
            sparse_struct.state_gradients[i,:] .= sum(sparse_struct.jacobians[i],dims=1)[:]
            update_turbine_powers!(sparse_struct,i)
        end
    end

    farm.AEP .= sum(sparse_struct.turbine_powers)
    farm.AEP_gradient .= sum(sparse_struct.state_gradients,dims=1)[:]

    return farm.AEP[1], farm.AEP_gradient
end

function update_turbine_powers!(sparse_struct,i)
    n = length(sparse_struct.caches[i].cache.fx)
    for j = 1:n
        sparse_struct.turbine_powers[j,i] = sparse_struct.caches[i].cache.fx[j].value
    end
end

# ∇AEP Optimization (update pattern) #######################################################

# Sparse spacing constraint methods #######################################################

# Sparse boundary constraint methods #######################################################
struct sparse_boundary_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9} <: AbstractSparseMethod
    turbine_x::T1
    turbine_y::T2
    jacobian::T3
    ad::T4
    cache::T5
    boundary_vec::T6
    boundary_function::T7
    update_function::T8
    boundary_scaling_factor::T9
end
