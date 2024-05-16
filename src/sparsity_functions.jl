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
                AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16,
                force_single_thread=false)

    farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type="ForwardDiff",force_single_thread=force_single_thread)

    sparse_struct = build_stable_sparse_struct(x,farm;tolerance=tolerance)

    farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type=eltype(sparse_struct.caches[1].cache.t),
                force_single_thread=force_single_thread)

    return farm,sparse_struct
end

function build_stable_sparse_struct(x,farm;tolerance=1E-16)
    n_states = length(farm.constants.wind_resource.wind_probabilities)
    n_turbines = length(farm.turbine_x)
    pow = zeros(eltype(x),n_turbines,n_states)
    jacobians = Array{SparseMatrixCSC{eltype(x), Int64},1}(undef,n_states)
    state_gradients = zeros(eltype(x),n_states,length(x))
    caches = nothing
    adtype = AutoSparseForwardDiff()

    define_patterns!(jacobians,x,farm,tolerance,pow,n_states)

    for i = 1:n_states
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

function define_patterns!(jacobians,x,farm,tolerance,pow,n_states)
    n_threads = Threads.nthreads()
    if n_threads > 1 && !farm.force_single_thread
        n_per_thread, rem = divrem(n_states,n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:n_states
        Threads.@threads for i_assignment in eachindex(assignments)
            i_start = assignments[i_assignment]
            i_stop = min(i_start+n_per_thread-1, n_states)
            for i = i_start:i_stop
                jacobians[i] = define_stable_jacobian_pattern(x,farm,tolerance,pow[:,i],i;prealloc_id=i_assignment)
            end
        end
    else
        for i = 1:n_states
            jacobians[i] = define_stable_jacobian_pattern(x,farm,tolerance,pow[:,i],i)
        end
    end
    return nothing
end

function define_stable_jacobian_pattern(x,farm,tolerance,pow,state_id;prealloc_id=1)
    p(a,x) = calculate_wind_state_power!(a,x,farm,state_id;prealloc_id=prealloc_id)
    x_temp = similar(x)
    cfg = ForwardDiff.JacobianConfig(nothing,pow,x)
    jacobian = ForwardDiff.jacobian(p,pow,x,cfg)
    jacobian[abs.(jacobian) .<= tolerance] .= 0.0
    jac_temp = similar(jacobian)
    for i = 1:3
        x_temp .= x .+ (rand(size(x_temp)) .- 0.5) * 1E-4
        jac_temp .= ForwardDiff.jacobian(p,pow,x,cfg)
        jac_temp[abs.(jac_temp) .<= tolerance] .= 0.0
        jacobian .+= jac_temp
    end
    return dropzeros(sparse(jacobian))
end

function calculate_wind_state_power!(pow,x,farm,state_id;prealloc_id=1,hours_per_year=365.25*24.0,lock=nothing)
    if !isnothing(lock)
        Threads.lock(lock)
    end

    farm.update_function(farm,x)
    turbine_x = copy(farm.turbine_x)
    turbine_y = copy(farm.turbine_y)
    rotor_diameter = copy(farm.rotor_diameter)
    turbine_yaw = copy(farm.turbine_yaw)
    hub_height = copy(farm.hub_height)

    if !isnothing(lock)
        Threads.unlock(lock)
    end

    rot_x, rot_y = rotate_to_wind_direction(turbine_x, turbine_y,
                    farm.constants.wind_resource.wind_directions[state_id])

    sorted_turbine_index = sortperm(rot_x)
    turbine_velocities = ff.turbine_velocities_one_direction(rot_x, rot_y, farm.constants.turbine_z,
                    rotor_diameter, hub_height, turbine_yaw, sorted_turbine_index,
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
                    rotor_diameter, turbine_velocities, farm.turbine_yaw, farm.constants.wind_resource.air_density,
                    farm.constants.power_models)

    pow .*= hours_per_year .* farm.constants.wind_resource.wind_probabilities[state_id] .* farm.AEP_scale
end

function calculate_aep_gradient!(farm,x,sparse_struct::T) where T <: StableSparseMethod
    n_threads = Threads.nthreads()
    n_states = length(farm.constants.wind_resource.wind_probabilities)

    if n_threads > 1 && !farm.force_single_thread && n_states > 1
        n_per_thread, rem = divrem(n_states,n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:n_states
        l = Threads.SpinLock()

        Threads.@threads for i_assignment in eachindex(assignments)
            i_start = assignments[i_assignment]
            i_stop = min(i_start+n_per_thread-1, n_states)
            for i = i_start:i_stop
                p(a,x) = calculate_wind_state_power!(a,x,farm,i;prealloc_id=i_assignment,lock=l)
                SparseDiffTools.sparse_jacobian!(sparse_struct.jacobians[i],sparse_struct.adtype,
                                sparse_struct.caches[i],p,sparse_struct.turbine_powers[:,i],x)
                sparse_struct.state_gradients[i,:] .= sum(sparse_struct.jacobians[i],dims=1)[:]
                update_turbine_powers!(sparse_struct,i)
            end
        end
    else
        for i = 1:n_states
            p(a,x) = calculate_wind_state_power!(a,x,farm,i;prealloc_id=1)
            SparseDiffTools.sparse_jacobian!(sparse_struct.jacobians[i],sparse_struct.adtype,
                            sparse_struct.caches[i],p,sparse_struct.turbine_powers[:,i],x)
            sparse_struct.state_gradients[i,:] .= sum(sparse_struct.jacobians[i],dims=1)[:]
            update_turbine_powers!(sparse_struct,i)
        end
    end

    farm.AEP .= sum(sparse_struct.turbine_powers)
    farm.AEP_gradient .= sum(sparse_struct.state_gradients,dims=1)[:]

    return farm.AEP[1], farm.AEP_gradient
end

function update_turbine_powers!(sparse_struct::T,i) where T <: StableSparseMethod
    n = length(sparse_struct.caches[i].cache.fx)
    for j = 1:n
        sparse_struct.turbine_powers[j,i] = sparse_struct.caches[i].cache.fx[j].value
    end
end

# ∇AEP Optimization (unstable pattern) #######################################################

struct sparse_AEP_struct_unstable_pattern{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10} <: UnstableSparseMethod
    deficit_thresholds::T1
    patterns::T2
    state_gradients::T3
    jacobians::T4
    turbine_powers::T5
    farm::T6 #farm of floats
    old_patterns::T7
    colors::T8
    state_powers::T9
    chunksize::T10
end

function build_unstable_sparse_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
                ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
                power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
                AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16,
                force_single_thread=false)

    farm_floats = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type=nothing,force_single_thread=force_single_thread)

    farm_forwarddiff = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type="ForwardDiffJacobian",force_single_thread=force_single_thread)

    sparse_struct, cache = build_unstable_sparse_struct(x,farm_floats,farm_forwarddiff;tolerance=tolerance)

    farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type=eltype(cache.t),force_single_thread=force_single_thread)

    return farm, sparse_struct

end

function build_unstable_sparse_struct(x,farm,farm_forwarddiff;tolerance=1E-16)
    n_states = length(farm.constants.wind_resource.wind_probabilities)
    n_turbines = length(farm.turbine_x)
    pow = zeros(eltype(x),n_turbines,n_states)
    jacobians = Array{SparseMatrixCSC{eltype(x), Int64},1}(undef,n_states)
    state_gradients = zeros(eltype(x),n_states,length(x))
    thresholds = zeros(eltype(x),n_states)
    patterns = zeros(eltype(x),n_turbines,length(x),n_states)
    old_patterns = zeros(eltype(x),n_turbines,length(x),n_states)
    colors = zeros(Int64,length(x),n_states)
    state_powers = zeros(eltype(x),n_states)

    calculate_thresholds!(jacobians,thresholds,x,farm_forwarddiff,farm,tolerance,pow,n_states)

    for i = 1:n_states
        colors[:,i] .= matrix_colors(jacobians[i])
    end

    cache = ForwardColorJacCache(nothing,x;
                    dx = pow[:,1],
                    colorvec = colors[:,1],
                    sparsity = jacobians[1])

    return sparse_AEP_struct_unstable_pattern(thresholds,patterns,state_gradients,jacobians,pow,farm,old_patterns,colors,state_powers,cache.chunksize), cache
end

function calculate_thresholds!(jacobians,thresholds,x,farm_forwarddiff,farm,tolerance,pow,n_states)
    n_threads = Threads.nthreads()

    if n_threads > 1 && !farm.force_single_thread
        n_per_thread, rem = divrem(n_states,n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:n_states
        l = Threads.SpinLock()
        Threads.@threads for i_assignment in eachindex(assignments)
            i_start = assignments[i_assignment]
            i_stop = min(i_start+n_per_thread-1, n_states)
            for i = i_start:i_stop
                jacobians[i],thresholds[i] = calculate_threshold(x,farm_forwarddiff,farm,tolerance,pow[:,i],i;prealloc_id=i_assignment,lock=l)
            end
        end
    else
        for i = 1:n_states
            jacobians[i],thresholds[i] = calculate_threshold(x,farm_forwarddiff,farm,tolerance,pow[:,i],i)
        end
    end
end

function calculate_threshold(x,farm_forwarddiff,farm,tolerance,pow,state_id;prealloc_id=1,lock=nothing)
    n_turbines = length(farm.turbine_x)
    jacobian = zeros(eltype(x),n_turbines,length(x))
    if farm.constants.wind_resource.wind_speeds[state_id] == 0.0 || farm.constants.wind_resource.wind_probabilities[state_id] == 0.0
        return sparse(jacobian), 0.0
    end
    x_temp = x
    n_variables = length(x)÷n_turbines
    p(a,x) = calculate_wind_state_power!(a,x,farm_forwarddiff,state_id;prealloc_id=prealloc_id,lock=lock)
    cfg = deepcopy(farm_forwarddiff.config)
    ForwardDiff.jacobian!(jacobian,p,pow,x_temp,cfg)
    jacobian[abs.(jacobian) .< tolerance] .= 0.0
    calculate_wind_state_power!(pow,x_temp,farm,state_id;prealloc_id=prealloc_id,lock=lock)
    deficits = view(farm.preallocations.prealloc_wake_deficits,:,:,prealloc_id)
    pattern = zeros(eltype(x),size(deficits))
    jac = deepcopy(reshape(jacobian,n_turbines,n_turbines,n_variables))
    for j = 1:n_turbines, i = 1:n_turbines
        for k = 1:n_variables
            if !iszero(jac[i,j,k])
                pattern[i,j] = 1.0
                break
            end
        end
    end
    deficits[pattern .== 0.0] .= 0.0
    deficits[deficits .== 0.0] .= Inf
    threshold = minimum(deficits) .* 1E-1
    return dropzeros(sparse(jacobian)), threshold
end

function calculate_aep_gradient!(farm,x,sparse_struct::T) where T <: UnstableSparseMethod
    n_threads = Threads.nthreads()
    n_states = length(farm.constants.wind_resource.wind_probabilities)

    if n_threads > 1 && !farm.force_single_thread && n_states > 1
        n_per_thread, rem = divrem(n_states,n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:n_states
        l = Threads.SpinLock()
        Threads.@threads for i_assignment in eachindex(assignments)
            i_start = assignments[i_assignment]
            i_stop = min(i_start+n_per_thread-1, n_states)
            for i = i_start:i_stop
                unstable_sparse_aep_gradient!(sparse_struct,x,farm,i;prealloc_id=i_assignment,lock=l)
            end
        end
    else
        for i = 1:n_states
            unstable_sparse_aep_gradient!(sparse_struct,x,farm,i)
        end
    end

    farm.AEP[1] = sum(sparse_struct.state_powers)
    farm.AEP_gradient .= sum(sparse_struct.state_gradients,dims=1)[:]

    return farm.AEP[1], farm.AEP_gradient
end

function unstable_sparse_aep_gradient!(sparse_struct::T,x,farm,wind_state_id;prealloc_id=1,lock=nothing) where T <: UnstableSparseMethod
    if farm.constants.wind_resource.wind_speeds[wind_state_id] == 0.0 || farm.constants.wind_resource.wind_probabilities[wind_state_id] == 0.0
        return
    end

    # calculate deficits and powers
    pow = view(sparse_struct.turbine_powers,:,wind_state_id)
    calculate_wind_state_power!(pow,x,sparse_struct.farm,wind_state_id;prealloc_id=prealloc_id,lock=lock)

    # calculate sparsity pattern
    calculate_unstable_sparsity_pattern!(sparse_struct,x,wind_state_id,prealloc_id)

    # calculate sparse jacobian
    calculate_unstable_sparse_jacobian!(sparse_struct,x,farm,wind_state_id,prealloc_id,lock)

    # sum turbine powers into state power
    sparse_struct.state_powers[wind_state_id] = sum(pow)

    # sum jacobian
    sparse_struct.state_gradients[wind_state_id,:] .= sum(sparse_struct.jacobians[wind_state_id],dims=1)[:]

    return nothing
end

function calculate_unstable_sparsity_pattern!(sparse_struct::T,x,wind_state_id,prealloc_id) where T <: UnstableSparseMethod
    n_turbines = length(sparse_struct.farm.turbine_x)
    n_variables = length(x)÷n_turbines
    pattern = view(sparse_struct.patterns,:,:,wind_state_id)
    pattern .= 0
    deficits = view(sparse_struct.farm.preallocations.prealloc_wake_deficits,:,:,prealloc_id)
    deficits[deficits .< sparse_struct.deficit_thresholds[wind_state_id]] .= 0.0

    for j in 1:n_turbines, i in 1:n_turbines
        if i == j
            for k = 1:n_variables
                pattern[i,j+(k-1)*n_turbines] = 1.0
            end
        else
            if !iszero(deficits[i,j])
                for k = 1:n_variables
                    pattern[i,j+(k-1)*n_turbines] = 1.0
                end
            end
        end
    end

    sparse_struct.jacobians[wind_state_id] .= dropzeros(sparse(pattern))

    # recolor if necessary
    if sparse_struct.patterns[:,:,wind_state_id] != sparse_struct.old_patterns[:,:,wind_state_id]
        recolor_jacobian!(sparse_struct,wind_state_id,n_variables,n_turbines)
        sparse_struct.old_patterns[:,:,wind_state_id] .= sparse_struct.patterns[:,:,wind_state_id]
    end
    return nothing
end

function recolor_jacobian!(sparse_struct::T,wind_state_id,n_variables,n_turbines) where T <: UnstableSparseMethod
    start_i(x) = n_turbines*(x-1)+1
    stop_i(x) = start_i(x) + n_turbines-1
    max_color = 0
    for i = 1:n_variables
        if i == 1
            sparse_struct.colors[start_i(i):stop_i(i),wind_state_id] .= matrix_colors(sparse_struct.jacobians[wind_state_id][:,start_i(i):stop_i(i)]) .+ max_color
        else
            sparse_struct.colors[start_i(i):stop_i(i),wind_state_id] .= sparse_struct.colors[start_i(i-1):stop_i(i-1),wind_state_id] .+ max_color
        end
        max_color = maximum(sparse_struct.colors[start_i(i):stop_i(i),wind_state_id])
    end
end

function calculate_unstable_sparse_jacobian!(sparse_struct::T,x,farm,wind_state_id,prealloc_id,lock) where T <: UnstableSparseMethod
    p(a,x) = calculate_wind_state_power!(a,x,farm,wind_state_id;prealloc_id=prealloc_id,lock=lock)

    cache = ForwardColorJacCache(nothing,x,sparse_struct.chunksize;
                              dx = sparse_struct.turbine_powers[:,wind_state_id],
                              colorvec=sparse_struct.colors[:,wind_state_id],
                              sparsity = sparse_struct.jacobians[wind_state_id])

    forwarddiff_color_jacobian!(sparse_struct.jacobians[wind_state_id],p,x,cache)
end

# Sparse spacing constraint methods ########################################################
struct sparse_spacing_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10} <: AbstractSparseMethod
    turbine_x::T1
    turbine_y::T2
    constraint_spacing::T3 # Single float that defines the minimum spacing between turbines in meters
    constraint_scaling::T4 # Single float that scales the constraint
    spacing_vec::T5 # In place vector
    jacobian::T6
    cache::T7
    update_function::T8
    relevant_list::T9
    ad::T10
end

function build_sparse_spacing_struct(x,turbine_x,turbine_y,space,scale,update_function;first_opt=true,relevant_spacing_factor=2)
    if first_opt
        n_constraints = 0
        relevant_list = nothing
        spacing_vec = zeros(eltype(x),n_constraints)
        spacing_jacobian = zeros(eltype(x),n_constraints,length(x))
        cache = nothing
        ad = nothing
        return sparse_spacing_struct(turbine_x,turbine_y,space,scale,spacing_vec,spacing_jacobian,cache,update_function,relevant_list,ad)
    end

    relevant_list,idx = build_relevant_list(turbine_x,turbine_y,space,relevant_spacing_factor)
    n_constraints = size(relevant_list,1)
    spacing_vec = zeros(eltype(x),n_constraints)

    # calculate jacobian pattern
    s_struct = build_spacing_struct(x,length(turbine_x),space,scale,update_function)
    calculate_spacing_jacobian!(s_struct,x)
    spacing_jacobian = dropzeros(sparse(s_struct.jacobian[idx,:]))

    ad = AutoSparseForwardDiff()
    sd = JacPrototypeSparsityDetection(; jac_prototype=spacing_jacobian)
    cache = sparse_jacobian_cache(ad, sd, nothing, spacing_vec, x)
    T = eltype(cache.cache.t)

    turbine_x = Vector{T}(turbine_x)
    turbine_y = Vector{T}(turbine_y)

    return sparse_spacing_struct(turbine_x,turbine_y,space,scale,spacing_vec,spacing_jacobian,cache,update_function,relevant_list,ad)
end

function build_relevant_list(turbine_x,turbine_y,space,factor)
    n_turbines = length(turbine_x)
    spacing = turbine_spacing(turbine_x,turbine_y)
    n_constraints = length(spacing)
    relevant = zeros(Int64,n_constraints,2)
    k = 1
    for i = 1:n_turbines
        for j = i+1:n_turbines
            relevant[k,1] = j
            relevant[k,2] = i
            k += 1
        end
    end

    for i = axes(relevant,1)
        if spacing[i] >= space*factor
            spacing[i] = Inf
        end
    end
    idx = spacing .!= Inf
    relevant = relevant[idx,:]
    return relevant,idx
end

function calculate_spacing!(spacing_vec,x,spacing_struct::T) where T<: AbstractSparseMethod
    spacing_struct.update_function(spacing_struct,x)

    sparse_spacing!(spacing_vec,spacing_struct.turbine_x,spacing_struct.turbine_y,spacing_struct.relevant_list)

    spacing_vec .= (spacing_struct.constraint_spacing .- spacing_vec) .* spacing_struct.constraint_scaling
end

function sparse_spacing!(spacing_vec,turbine_x,turbine_y,relevant)
    for i in axes(relevant,1)
        j = relevant[i,1]
        k = relevant[i,2]
        spacing_vec[i] = sqrt((turbine_x[j] - turbine_x[k])^2+(turbine_y[j] - turbine_y[k])^2)
    end
    return spacing_vec
end

function calculate_spacing_jacobian!(spacing_struct::T,x) where T <: AbstractSparseMethod
    if isnothing(spacing_struct.cache)
        return spacing_struct.spacing_vec, spacing_struct.jacobian
    end
    calculate_spacing(a,b) = calculate_spacing!(a,b,spacing_struct)
    sparse_jacobian!(spacing_struct.jacobian,spacing_struct.ad,spacing_struct.cache,calculate_spacing,spacing_struct.spacing_vec,x)

    for i = eachindex(spacing_struct.cache.cache.fx)
        spacing_struct.spacing_vec[i] = spacing_struct.cache.cache.fx[i].value
    end

    return spacing_struct.spacing_vec, spacing_struct.jacobian
end

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

function calculate_boundary_jacobian!(boundary_struct::T,x) where T <: AbstractSparseMethod
    calculate_boundary(a,b) = calculate_boundary!(a,b,boundary_struct)
    sparse_jacobian!(boundary_struct.jacobian, boundary_struct.ad, boundary_struct.cache, calculate_boundary, boundary_struct.boundary_vec, x)

    for i = eachindex(boundary_struct.cache.cache.fx)
        boundary_struct.boundary_vec[i] = boundary_struct.cache.cache.fx[i].value
    end

    return boundary_struct.boundary_vec, boundary_struct.jacobian
end
