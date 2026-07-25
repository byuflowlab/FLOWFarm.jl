export build_stable_sparse_struct, build_unstable_sparse_struct, build_sparse_spacing_struct,
calculate_aep_gradient!, calculate_spacing_jacobian!, calculate_boundary_jacobian!, calculate_spacing!

"""file with functions used in wind farm optimization employing sparse methods
created January 26, 2024
author: Benjamin Varela
"""

abstract type AbstractSparseMethod end
abstract type StableSparseMethod <: AbstractSparseMethod end
abstract type UnstableSparseMethod <: AbstractSparseMethod end

# ∇AEP Optimization (stable pattern) #######################################################
"""
sparse_AEP_struct_stable_pattern

Struct that holds all the necessary variables to calculate the AEP gradient using a stable sparsity pattern

# Arguments
- `preps`: vector of Ref cells, one per wind state, holding the DifferentiationInterface Jacobian
  preparation for that state (built lazily on first use, `nothing` beforehand)
- `jacobians`: vector of sparse matracies containing jacobians for each wind state
- `state_gradients`: 2d array, each row is a state gradient (used for threads)
- `turbine_powers`: 2d array that holds the powers for each turbine (used for threads)
- `adtypes`: vector of AutoSparse(AutoForwardDiff()) objects, one per wind state, all sharing a
  common chunk width so that the wind farm's preallocated dual buffers stay valid across states
- `value_farm`: a plain (non-dual) WindFarm struct used to compute the actual turbine power
  values each call; kept separate from the dual-typed farm used for the jacobian sweep because
  `jacobian!` never has to fall back to a plain-Float64 evaluation internally (unlike
  `value_and_jacobian!`, which does, and which is therefore incompatible with a farm whose
  preallocated buffers are permanently dual-typed)
"""
struct sparse_AEP_struct_stable_pattern{T1,T2,T3,T4,T5,T6} <: StableSparseMethod
    preps::T1 # vector of Ref cells holding DI jacobian preparations
    jacobians::T2 # vector of sparse jacobians
    state_gradients::T3 # 2d array, each row is a state gradient (used for threads)
    turbine_powers::T4 # 2d array that holds the powers or each turbine (used for threads)
    adtypes::T5
    value_farm::T6
end

"""
build_stable_sparse_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
                ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
                power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
                AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16)

Function that builds a wind_farm_struct and a sparse_AEP_struct_stable_pattern struct that goes with it

# Arguments
- `x`: Vector containing the design variables
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
- `tolerance`: Single float that defines the tolerance for the jacobian pattern (default is 1E-16), set to 0.0 to use traditional sparsity
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for the jacobian coloring (default is `GreedyColoringAlgorithm()`)
"""
function build_stable_sparse_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
                ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
                power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
                AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16,
                coloring_algorithm=GreedyColoringAlgorithm())

    probe_farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type="ForwardDiff")

    n_states = length(probe_farm.constants.wind_resource.wind_probabilities)
    n_turbines = length(probe_farm.turbine_x)
    pow = zeros(eltype(x),n_turbines,n_states)
    jacobians = Array{SparseMatrixCSC{eltype(x), Int64},1}(undef,n_states)
    define_patterns!(jacobians,x,probe_farm,tolerance,pow,n_states)

    # chunk width shared by every wind state, so the farm's preallocated dual buffers
    # (built once, below) stay valid no matter which state is being differentiated;
    # derived directly from each state's jacobian sparsity pattern, no differentiated
    # call needed
    n_colors = shared_chunk_width(jacobians,n_states,coloring_algorithm)
    T_dual = dual_type(_wind_state_power_ad!, eltype(x), n_colors)

    farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type=T_dual)

    # plain (non-dual) farm used only to evaluate the actual turbine power values each call
    value_farm = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam)

    sparse_struct = build_stable_sparse_struct(jacobians,pow,x,n_states,value_farm;coloring_algorithm=coloring_algorithm)

    return farm, sparse_struct
end

"""
build_stable_sparse_struct(x,farm;tolerance=1E-16)

Helper function that builds a sparse_AEP_struct_stable_pattern struct

Note: `farm` here is used both to detect the jacobian sparsity pattern and as the struct's
`value_farm` (the plain, non-dual evaluator used to get actual turbine power values each call),
so it must be a genuinely plain (non-dual) WindFarm struct. The dual-typed farm actually used
for differentiation is a separate object the caller must build (with a `ForwardDiff.Dual` chunk
width matching `shared_chunk_width` of the resulting jacobian patterns and tagged to
`_wind_state_power_ad!`) and pass to `calculate_aep_gradient!` directly - it is not this `farm`.
The full-argument `build_stable_sparse_struct` method builds and wires up both farms correctly
and should be preferred; use this method only if you already know what you're doing.

# Arguments
- `x`: Vector containing the  design variables
- `farm`: plain (non-dual) WindFarm struct
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for the jacobian coloring (default is `GreedyColoringAlgorithm()`)
"""
function build_stable_sparse_struct(x,farm;tolerance=1E-16,coloring_algorithm=GreedyColoringAlgorithm())
    n_states = length(farm.constants.wind_resource.wind_probabilities)
    n_turbines = length(farm.turbine_x)
    pow = zeros(eltype(x),n_turbines,n_states)
    jacobians = Array{SparseMatrixCSC{eltype(x), Int64},1}(undef,n_states)

    define_patterns!(jacobians,x,farm,tolerance,pow,n_states)

    return build_stable_sparse_struct(jacobians,pow,x,n_states,farm;coloring_algorithm=coloring_algorithm)
end

"""
build_stable_sparse_struct(jacobians,pow,x,n_states,value_farm;coloring_algorithm=GreedyColoringAlgorithm())

Helper function that builds a sparse_AEP_struct_stable_pattern struct from already-computed
jacobian sparsity patterns

# Arguments
- `jacobians`: Vector of sparse arrays holding the jacobian pattern for each wind state
- `pow`: 2d array that holds the powers or each turbine (used for threads)
- `x`: Vector containing the design variables
- `n_states`: Number of wind states
- `value_farm`: a plain (non-dual) WindFarm struct used to compute actual turbine power values
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for the jacobian coloring (default is `GreedyColoringAlgorithm()`)
"""
function build_stable_sparse_struct(jacobians,pow,x,n_states,value_farm;coloring_algorithm=GreedyColoringAlgorithm())
    state_gradients = zeros(eltype(x),n_states,length(x))
    preps = [Ref{Any}(nothing) for _ = 1:n_states]

    # every state must share one chunk width, since they all read/write through the same
    # (thread-indexed, not state-indexed) preallocated dual buffers on the wind farm struct
    n_colors = shared_chunk_width(jacobians,n_states,coloring_algorithm)
    adtypes = [AutoSparse(AutoForwardDiff(chunksize=n_colors); sparsity_detector=KnownJacobianSparsityDetector(jacobians[i]),
                    coloring_algorithm=coloring_algorithm) for i = 1:n_states]

    return sparse_AEP_struct_stable_pattern(preps,jacobians,state_gradients,pow,adtypes,value_farm)
end

"""
shared_chunk_width(jacobians,n_states,coloring_algorithm)

Helper function that returns the largest ForwardDiff chunk width that remains valid
(<= number of colors) for every wind state's jacobian sparsity pattern, computed directly
from the patterns themselves (no differentiated function call required)
"""
function shared_chunk_width(jacobians,n_states,coloring_algorithm)
    return minimum(maximum(column_colors(coloring(jacobians[i], ColoringProblem(), coloring_algorithm))) for i = 1:n_states)
end

"""
dual_type(f,V,N)

Helper function that builds the exact `ForwardDiff.Dual` type DifferentiationInterface will use
to differentiate the top-level function `f` (with primal element type `V`) at chunk width `N`.

Note: this deliberately does NOT go through `ForwardDiff.JacobianConfig(f,y,x,ForwardDiff.Chunk(N))`
- `ForwardDiff.Chunk(N::Integer)` silently clamps `N` down to `ForwardDiff.DEFAULT_CHUNK_THRESHOLD`
  (12), which would silently produce farm preallocations narrower than the chunk width DI actually
  uses whenever a jacobian's coloring needs more than 12 simultaneous colors - DI's own chunk
  selection is not subject to that cap. Constructing the tag/type directly avoids it.
"""
function dual_type(f,V,N)
    return ForwardDiff.Dual{typeof(ForwardDiff.Tag(f,V)),V,N}
end

"""
define_patterns!(jacobians,x,farm,tolerance,pow,n_states)

Helper function that defines the jacobian patterns for each wind state

# Arguments
- `jacobians`: Vector of sparse arrays holding jacobians
- `x`: Vector containing the  design variables
- `farm`: WindFarm struct
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
- `pow`: 2d array that holds the powers or each turbine (used for threads)
- `n_states`: Number of wind states
"""
function define_patterns!(jacobians,x,farm,tolerance,pow,n_states)
    n_threads = Threads.nthreads()
    if n_threads > 1
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

"""
define_stable_jacobian_pattern(x,farm,tolerance,pow,state_id;prealloc_id=1)

Helper function that defines the jacobian pattern for a single wind state

# Arguments
- `x`: Vector containing the  design variables
- `farm`: WindFarm struct
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
- `pow`: 1d array that holds the powers or each turbine
- `state_id`: Wind state id
- `prealloc_id`: Preallocation id (to select the correct preallocated memory inside the wind farm struct)
"""
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

"""
calculate_wind_state_power!(pow,x,farm,state_id;prealloc_id=1,hours_per_year=365.25*24.0,lock=nothing)

Helper function that calculates the power fora a single wind state

# Arguments
- `pow`: 1d array that holds the powers or each turbine
- `x`: Vector containing the  design variables
- `farm`: WindFarm struct
- `state_id`: Wind state id
- `prealloc_id`: Preallocation id (to select the correct preallocated memory inside the wind farm struct)
- `hours_per_year`: Single float that defines the hours per year
- `lock`: SpinLock object to lock the farm struct for multithreadeding
"""
function calculate_wind_state_power!(pow,x,farm,state_id;prealloc_id=1,hours_per_year=365.25*24.0,lock=nothing)
    if below_cutin_speed(farm.constants.wind_resource, state_id, farm.constants.turbine_z,
                         farm.hub_height, farm.constants.cut_in_speed)
        pow .= 0
        return pow
    end

    if !isnothing(lock)
        Threads.lock(lock)
    end

    farm.update_function(farm,x)
    prealloc_rot_x = view(farm.preallocations.prealloc_rot_x, :, prealloc_id)
    prealloc_rot_y = view(farm.preallocations.prealloc_rot_y, :, prealloc_id)
    rot_x, rot_y = rotate_to_wind_direction!(prealloc_rot_x, prealloc_rot_y, farm.turbine_x, farm.turbine_y,
                    farm.constants.wind_resource.wind_directions[state_id])
    
    for i in 1:length(rot_x)
        farm.preallocations.prealloc_diam[i, prealloc_id] = farm.rotor_diameter[i]
        farm.preallocations.prealloc_yaw[i, prealloc_id] = farm.turbine_yaw[i]
        farm.preallocations.prealloc_hub[i, prealloc_id] = farm.hub_height[i]
    end
    rotor_diameter = view(farm.preallocations.prealloc_diam, :, prealloc_id)
    turbine_yaw = view(farm.preallocations.prealloc_yaw, :, prealloc_id)
    hub_height = view(farm.preallocations.prealloc_hub, :, prealloc_id)

    if !isnothing(lock)
        Threads.unlock(lock)
    end

    sorted_turbine_index = view(farm.preallocations.prealloc_sort_index, :, prealloc_id)
    sortperm!(sorted_turbine_index, rot_x)
    turbine_velocities = turbine_velocities_one_direction_vel(rot_x, rot_y, farm.constants.turbine_z,
                    rotor_diameter, hub_height, turbine_yaw, sorted_turbine_index,
                    farm.constants.ct_models, farm.constants.rotor_sample_points_y, farm.constants.rotor_sample_points_z,
                    farm.constants.wind_resource, farm.constants.model_set, wind_farm_state_id=state_id,
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
                    farm.constants.power_models; wt_power=pow)

    pow .*= hours_per_year .* farm.constants.wind_resource.wind_probabilities[state_id] .* farm.AEP_scale
end

"""
_wind_state_power_ad!(pow,x,farm,state_id,prealloc_id,lock)

Top-level (non-closure) wrapper around `calculate_wind_state_power!` used as the function
differentiated by DifferentiationInterface. `farm`, `state_id`, `prealloc_id`, and `lock` are
passed as DI `Constant` contexts rather than captured in a closure: DI ties a Jacobian
preparation to the exact type of the function being differentiated, and a plain closure's type
would depend on `typeof(farm)`, which itself depends on the dual type baked into `farm`'s
preallocated arrays - a circular requirement. Using this fixed top-level function instead makes
the differentiation tag independent of `farm`'s type entirely.
"""
function _wind_state_power_ad!(pow,x,farm,state_id,prealloc_id,lock)
    return calculate_wind_state_power!(pow,x,farm,state_id;prealloc_id=prealloc_id,lock=lock)
end

"""
calculate_aep_gradient!(farm,x,sparse_struct::T)

Function that calculates the AEP gradient using a stable sparsity pattern

# Arguments
- `farm`: WindFarm struct
- `x`: Vector containing the scaled design variables
- `sparse_struct`: sparse_AEP_struct_stable_pattern struct
"""
function calculate_aep_gradient!(farm,x,sparse_struct::T) where T <: StableSparseMethod
    n_threads = Threads.nthreads()
    n_states = length(farm.constants.wind_resource.wind_probabilities)

    if n_threads > 1 && n_states > 1
        calculate_aep_gradient_multithreads!(farm, x, sparse_struct, n_threads, n_states)
    else
        @views @inbounds for i = 1:n_states
            # actual power value: cheap plain (non-dual) evaluation against the plain value_farm
            calculate_wind_state_power!(sparse_struct.turbine_powers[:,i],x,sparse_struct.value_farm,i;prealloc_id=1)

            # jacobian: dual evaluation against the caller-supplied dual-typed farm
            if isnothing(sparse_struct.preps[i][])
                sparse_struct.preps[i][] = prepare_jacobian(_wind_state_power_ad!, sparse_struct.turbine_powers[:,i],
                                sparse_struct.adtypes[i], x, Constant(farm), Constant(i), Constant(1), Constant(nothing))
            end
            jacobian!(_wind_state_power_ad!, sparse_struct.turbine_powers[:,i], sparse_struct.jacobians[i],
                            sparse_struct.preps[i][], sparse_struct.adtypes[i], x, Constant(farm), Constant(i), Constant(1), Constant(nothing))
            sum_jacobians!(sparse_struct,i)
        end
    end

    # Avoid vector operations for assignment
    total_aep = zero(eltype(x))
    for i = 1:size(sparse_struct.turbine_powers, 1), j = 1:size(sparse_struct.turbine_powers, 2)
        total_aep += sparse_struct.turbine_powers[i, j]
    end
    farm.AEP[1] = total_aep

    for j = 1:size(sparse_struct.state_gradients, 2)
        s = zero(eltype(x))
        for i = 1:size(sparse_struct.state_gradients, 1)
            s += sparse_struct.state_gradients[i, j]
        end
        farm.AEP_gradient[j] = s
    end

    return farm.AEP[1], farm.AEP_gradient
end

function sum_jacobians!(sparse_struct::T, state_idx) where T <: StableSparseMethod
    for j = 1:size(sparse_struct.jacobians[state_idx], 2)
        s = zero(eltype(sparse_struct.jacobians[state_idx]))
        for k = 1:size(sparse_struct.jacobians[state_idx], 1)
            s += sparse_struct.jacobians[state_idx][k, j]
        end
        sparse_struct.state_gradients[state_idx, j] = s
    end
end

function calculate_aep_gradient_multithreads!(farm, x, sparse_struct::T, n_threads, n_states) where T <: StableSparseMethod
    n_per_thread, rem = divrem(n_states,n_threads)
    n = n_per_thread + (rem > 0)
    assignments = 1:n:n_states
    l = Threads.SpinLock()

    Threads.@threads for i_assignment in eachindex(assignments)
        i_start = assignments[i_assignment]
        i_stop = min(i_start+n-1, n_states)
        for i = i_start:i_stop
            calculate_wind_state_power!(view(sparse_struct.turbine_powers,:,i),x,sparse_struct.value_farm,i;
                            prealloc_id=i_assignment,lock=l)
            # note: value_farm and the dual-typed farm each have their own preallocation
            # buffers (separate structs), so this lock only needs to guard concurrent access
            # within a single one of them, not between the two calls above

            if isnothing(sparse_struct.preps[i][])
                sparse_struct.preps[i][] = prepare_jacobian(_wind_state_power_ad!, view(sparse_struct.turbine_powers,:,i),
                                sparse_struct.adtypes[i], x, Constant(farm), Constant(i), Constant(i_assignment), Constant(l))
            end
            jacobian!(_wind_state_power_ad!, view(sparse_struct.turbine_powers, :, i), sparse_struct.jacobians[i],
                            sparse_struct.preps[i][], sparse_struct.adtypes[i], x, Constant(farm), Constant(i), Constant(i_assignment), Constant(l))
            sum_jacobians!(sparse_struct,i)
        end
    end
end

# ∇AEP Optimization (unstable pattern) #######################################################

"""
sparse_AEP_struct_unstable_pattern

Struct that holds all the necessary variables to calculate the AEP gradient using an unstable sparsity pattern

# Arguments
- `deficit_thresholds`: Vector of floats that define the deficit thresholds for each wind state
- `patterns`: 3d array that holds the sparsity patterns for each wind state
- `state_gradients`: 2d array, each row is a state gradient (used for threads)
- `jacobians`: Vector of sparse arrays containing jacobians for each wind state
- `turbine_powers`: 2d array that holds the powers or each turbine (used for threads)
- `farm`: WindFarm struct
- `old_patterns`: 3d array that holds the old sparsity patterns for each wind state
- `preps`: vector of Ref cells, one per wind state, holding the DifferentiationInterface Jacobian
  preparation for that state (rebuilt whenever that state's sparsity pattern changes)
- `state_powers`: 1d array that holds the state powers
- `adtypes`: vector of Ref cells, one per wind state, holding the current AutoSparse(AutoForwardDiff())
  object for that state (rebuilt alongside `preps` whenever the pattern changes)
- `tier_widths`: Vector of Int, descending chunk widths (e.g. [8,4,2,1]) pre-built as candidate
  dual chunk widths, derived once from the initial shared color count across states
- `tier_farms`: Vector of WindFarm structs, one per entry in `tier_widths`, each dual-typed with
  that tier's chunk width baked into its preallocated buffers
- `current_tier`: Vector of Int, one per wind state, index into `tier_widths`/`tier_farms`
  selecting which pre-built tier is currently valid for that state's sparsity pattern
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for jacobian coloring, stored
  so that `recolor_jacobian!` reuses the same algorithm the struct was built with

Note: since the sparsity pattern can change at runtime and DI requires chunk width <= current
color count, no single chunk width is guaranteed valid for the whole run. Rather than falling
back to width 1 (fully sequential dual sweeps) for the entire run, a small fixed set of
pre-built farm/adtype tiers lets each state use the widest pre-built tier still valid for its
*current* pattern, recomputed whenever that pattern changes (see `recolor_jacobian!`). This adds
some upfront build cost and a little extra memory (one farm's preallocated buffers per tier) but
avoids ever rebuilding a farm at runtime.
"""
struct sparse_AEP_struct_unstable_pattern{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14} <: UnstableSparseMethod
    deficit_thresholds::T1
    patterns::T2
    state_gradients::T3
    jacobians::T4
    turbine_powers::T5
    farm::T6 #farm of floats
    old_patterns::T7
    preps::T8
    state_powers::T9
    adtypes::T10
    tier_widths::T11
    tier_farms::T12
    current_tier::T13
    coloring_algorithm::T14
end

"""
build_unstable_sparse_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
                ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
                power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
                AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16)

Function that builds a wind_farm_struct and a sparse_AEP_struct_unstable_pattern struct

# Arguments
- `x`: Vector containing the  design variables
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
- `tolerance`: Single float that defines the tolerance for the jacobian pattern (default is 1E-16), set to 0.0 to use traditional sparsity
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for the jacobian coloring (default is `GreedyColoringAlgorithm()`)
"""
function build_unstable_sparse_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,rotor_diameter,
                ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,rated_power,wind_resource,
                power_models,model_set,update_function;rotor_sample_points_y=[0.0],rotor_sample_points_z=[0.0],
                AEP_scale=0.0,opt_x=false,opt_y=false,opt_hub=false,opt_yaw=false,opt_diam=false,tolerance=1E-16,
                coloring_algorithm=GreedyColoringAlgorithm())

    farm_floats = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type=nothing)

    farm_forwarddiff = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                input_type="ForwardDiffJacobian")

    jacobians, thresholds, patterns, old_patterns, state_gradients, state_powers, pow, n_states, n_turbines =
                    build_unstable_sparse_arrays(x,farm_floats,farm_forwarddiff,tolerance)

    # small fixed set of candidate chunk widths, descending from the initial shared color count
    # down to 1 (see sparse_AEP_struct_unstable_pattern docstring); one dual-typed farm is built
    # per tier, each tagged to the same top-level function that will actually be differentiated
    # (_wind_state_power_ad!)
    tier_widths = build_tier_widths(jacobians,n_states,coloring_algorithm)
    tier_farms = Vector{Any}(undef,length(tier_widths))
    for (ti,w) = enumerate(tier_widths)
        T_dual = dual_type(_wind_state_power_ad!, eltype(x), w)
        tier_farms[ti] = build_wind_farm_struct(x,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
                    rotor_diameter,ct_models,generator_efficiency,cut_in_speed,
                    cut_out_speed,rated_speed,rated_power,wind_resource,power_models,
                    model_set,update_function;rotor_sample_points_y=rotor_sample_points_y,
                    rotor_sample_points_z=rotor_sample_points_z,AEP_scale=AEP_scale,
                    opt_x=opt_x,opt_y=opt_y,opt_hub=opt_hub,opt_yaw=opt_yaw,opt_diam=opt_diam,
                    input_type=T_dual)
    end

    preps = [Ref{Any}(nothing) for _ = 1:n_states]
    adtypes = [Ref{Any}(nothing) for _ = 1:n_states]
    current_tier = ones(Int,n_states)

    sparse_struct = sparse_AEP_struct_unstable_pattern(thresholds,patterns,state_gradients,jacobians,pow,farm_floats,
                    old_patterns,preps,state_powers,adtypes,tier_widths,tier_farms,current_tier,coloring_algorithm)

    return farm_floats, sparse_struct

end

"""
build_tier_widths(jacobians,n_states,coloring_algorithm)

Helper function that builds the descending list of candidate chunk widths for the
unstable-pattern struct, halving from the initial shared color count (the smallest color count
across all states, so every state has at least one valid tier to start from) down to 1.
"""
function build_tier_widths(jacobians,n_states,coloring_algorithm)
    w = minimum(maximum(column_colors(coloring(jacobians[i], ColoringProblem(), coloring_algorithm))) for i = 1:n_states)
    widths = Int[]
    while w > 1
        push!(widths,w)
        w = w ÷ 2
    end
    push!(widths,1)
    return unique(widths)
end

"""
build_unstable_sparse_struct(x,farm,farm_forwarddiff;tolerance=1E-16)

Helper function that builds a sparse_AEP_struct_unstable_pattern struct with a single chunk-width
tier (chunksize=1), using `farm` itself as that tier's farm.

Note: `farm`'s preallocated arrays must already be typed with a `ForwardDiff.Dual` chunk width of
1 tagged to `_wind_state_power_ad!`, or differentiation will error with a chunk-size mismatch.
This single-tier form exists for callers that only have one farm on hand; the full-argument
`build_unstable_sparse_struct` method builds multiple width tiers (see
`sparse_AEP_struct_unstable_pattern` docstring) and should be preferred for performance.

# Arguments
- `x`: Vector containing the  design variables
- `farm`: WindFarm struct
- `farm_forwarddiff`: WindFarm struct with ForwardDiff input type for deficit tolerance calculation
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for the jacobian coloring (default is `GreedyColoringAlgorithm()`)
"""
function build_unstable_sparse_struct(x,farm,farm_forwarddiff;tolerance=1E-16,coloring_algorithm=GreedyColoringAlgorithm())
    jacobians, thresholds, patterns, old_patterns, state_gradients, state_powers, pow, n_states, n_turbines =
                    build_unstable_sparse_arrays(x,farm,farm_forwarddiff,tolerance)

    preps = [Ref{Any}(nothing) for _ = 1:n_states]
    adtypes = [Ref{Any}(nothing) for _ = 1:n_states]
    current_tier = ones(Int,n_states)

    return sparse_AEP_struct_unstable_pattern(thresholds,patterns,state_gradients,jacobians,pow,farm,old_patterns,
                    preps,state_powers,adtypes,[1],[farm],current_tier,coloring_algorithm)
end

"""
build_unstable_sparse_arrays(x,farm,farm_forwarddiff,tolerance)

Helper function that computes the jacobian sparsity patterns/thresholds and allocates the plain
(non-tier-specific) arrays shared by `sparse_AEP_struct_unstable_pattern`, without constructing
the struct itself (since building the chunk-width tiers requires the raw farm-construction
arguments, not just `farm`/`farm_forwarddiff`).

# Arguments
- `x`: Vector containing the  design variables
- `farm`: WindFarm struct
- `farm_forwarddiff`: WindFarm struct with ForwardDiff input type for deficit tolerance calculation
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
"""
function build_unstable_sparse_arrays(x,farm,farm_forwarddiff,tolerance)
    n_states = length(farm.constants.wind_resource.wind_probabilities)
    n_turbines = length(farm.turbine_x)
    pow = zeros(eltype(x),n_turbines,n_states)
    jacobians = Array{SparseMatrixCSC{eltype(x), Int64},1}(undef,n_states)
    state_gradients = zeros(eltype(x),n_states,length(x))
    thresholds = zeros(eltype(x),n_states)
    patterns = zeros(eltype(x),n_turbines,length(x),n_states)
    old_patterns = zeros(eltype(x),n_turbines,length(x),n_states)
    state_powers = zeros(eltype(x),n_states)

    calculate_thresholds!(jacobians,thresholds,x,farm_forwarddiff,farm,tolerance,pow,n_states)

    return jacobians, thresholds, patterns, old_patterns, state_gradients, state_powers, pow, n_states, n_turbines
end

"""
calculate_thresholds!(jacobians,thresholds,x,farm_forwarddiff,farm,tolerance,pow,n_states)

Helper function that calculates the thresholds for each wind state

# Arguments
- `jacobians`: Vector of sparse arrays holding jacobians
- `thresholds`: Vector of floats that define the deficit thresholds for each wind state
- `x`: Vector containing the  design variables
- `farm_forwarddiff`: WindFarm struct with ForwardDiff input type for deficit tolerance calculation
- `farm`: WindFarm struct
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
- `pow`: 2d array that holds the powers or each turbine (used for threads)
- `n_states`: Number of wind states
"""
function calculate_thresholds!(jacobians,thresholds,x,farm_forwarddiff,farm,tolerance,pow,n_states)
    n_threads = Threads.nthreads()

    if n_threads > 1
        n_per_thread, rem = divrem(n_states,n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:n_states
        l = Threads.SpinLock()
        Threads.@threads for i_assignment = eachindex(assignments)
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

"""
calculate_threshold(x,farm_forwarddiff,farm,tolerance,pow,state_id;prealloc_id=1,lock=nothing)

Helper function that calculates the deficit threshold for a single wind state

# Arguments
- `x`: Vector containing the design variables
- `farm_forwarddiff`: WindFarm struct with ForwardDiff input type for deficit tolerance calculation
- `farm`: WindFarm struct
- `tolerance`: Single float that defines the tolerance for the jacobian pattern
- `pow`: 1d array that holds the powers or each turbine
- `state_id`: Wind state id
- `prealloc_id`: Preallocation id (to select the correct preallocated memory inside the wind farm struct)
- `lock`: SpinLock object to lock the farm struct for multithreadeding
"""
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

"""
calculate_aep_gradient!(farm,x,sparse_struct::T)

Function that calculates the AEP gradient using an unstable sparsity pattern

# Arguments
- `farm`: WindFarm struct
- `x`: Vector containing the scaled design variables
- `sparse_struct`: sparse_AEP_struct_unstable_pattern struct
"""
function calculate_aep_gradient!(farm,x,sparse_struct::T) where T <: UnstableSparseMethod
    n_threads = Threads.nthreads()
    n_states = length(farm.constants.wind_resource.wind_probabilities)

    if n_threads > 1 && n_states > 1
        calculate_aep_gradient_multithreads!(farm, x, sparse_struct, n_threads, n_states)
    else
        for i = 1:n_states
            unstable_sparse_aep_gradient!(sparse_struct,x,i)
        end
    end

    # Avoid vector operations for assignment
    total_aep = zero(eltype(x))
    for i = 1:size(sparse_struct.turbine_powers, 1), j = 1:size(sparse_struct.turbine_powers, 2)
        total_aep += sparse_struct.turbine_powers[i, j]
    end
    farm.AEP[1] = total_aep

    for j = 1:size(sparse_struct.state_gradients, 2)
        s = zero(eltype(x))
        for i = 1:size(sparse_struct.state_gradients, 1)
            s += sparse_struct.state_gradients[i, j]
        end
        farm.AEP_gradient[j] = s
    end

    ##HERE

    return farm.AEP[1], farm.AEP_gradient
end

function calculate_aep_gradient_multithreads!(farm, x, sparse_struct::T, n_threads, n_states) where T <: UnstableSparseMethod
    n_per_thread, rem = divrem(n_states,n_threads)
    n = n_per_thread + (rem > 0)
    assignments = 1:n:n_states
    l = Threads.SpinLock()

    Threads.@threads for i_assignment in eachindex(assignments)
        i_start = assignments[i_assignment]
        i_stop = min(i_start+n-1, n_states)
        for i = i_start:i_stop
            unstable_sparse_aep_gradient!(sparse_struct,x,i;prealloc_id=i_assignment,lock=l)
        end
    end
end

"""
unstable_sparse_aep_gradient!(sparse_struct::T,x,wind_state_id;prealloc_id=1,lock=nothing)

Function that calculates the AEP gradient for a single wind state using an unstable sparsity pattern

# Arguments
- `sparse_struct`: sparse_AEP_struct_unstable_pattern struct
- `x`: Vector containing the design variables
- `wind_state_id`: Wind state id
- `prealloc_id`: Preallocation id (to select the correct preallocated memory inside the wind farm struct)
- `lock`: SpinLock object to lock the farm struct for multithreadeding
"""
function unstable_sparse_aep_gradient!(sparse_struct::T,x,wind_state_id;prealloc_id=1,lock=nothing) where T <: UnstableSparseMethod
    if sparse_struct.farm.constants.wind_resource.wind_speeds[wind_state_id] == 0.0 || sparse_struct.farm.constants.wind_resource.wind_probabilities[wind_state_id] == 0.0
        return
    end

    # calculate deficits and powers
    pow = view(sparse_struct.turbine_powers,:,wind_state_id)
    calculate_wind_state_power!(pow,x,sparse_struct.farm,wind_state_id;prealloc_id=prealloc_id,lock=lock)

    # calculate sparsity pattern
    calculate_unstable_sparsity_pattern!(sparse_struct,x,wind_state_id,prealloc_id)

    # calculate sparse jacobian
    calculate_unstable_sparse_jacobian!(sparse_struct,x,wind_state_id,prealloc_id,lock)

    # sum turbine powers into state power
    sparse_struct.state_powers[wind_state_id] = sum(pow)

    # sum jacobian
    for j = 1:size(sparse_struct.jacobians[wind_state_id], 2)
        s = zero(eltype(x))
        for i = 1:size(sparse_struct.jacobians[wind_state_id], 1)
            s += sparse_struct.jacobians[wind_state_id][i, j]
        end
        sparse_struct.state_gradients[wind_state_id, j] = s
    end

    return nothing
end

"""
calculate_unstable_sparsity_pattern!(sparse_struct::T,x,wind_state_id,prealloc_id)

Helper function that calculates the sparsity pattern for a single wind state using an unstable sparsity pattern

# Arguments
- `sparse_struct`: sparse_AEP_struct_unstable_pattern struct
- `x`: Vector containing the design variables
- `wind_state_id`: Wind state id
- `prealloc_id`: Preallocation id (to select the correct preallocated memory inside the wind farm struct)
"""
function calculate_unstable_sparsity_pattern!(sparse_struct::T,x,wind_state_id,prealloc_id) where T <: UnstableSparseMethod
    n_turbines = length(sparse_struct.farm.turbine_x)
    n_variables = length(x)÷n_turbines
    pattern = view(sparse_struct.patterns,:,:,wind_state_id)
    pattern .= 0
    deficits = view(sparse_struct.farm.preallocations.prealloc_wake_deficits,:,:,prealloc_id)
    for i in eachindex(deficits)
        if deficits[i] < sparse_struct.deficit_thresholds[wind_state_id]
            deficits[i] = 0.0
        end
    end

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

    for row in axes(pattern, 1), col in axes(pattern, 2)
        sparse_struct.jacobians[wind_state_id][row, col] = pattern[row, col]
    end
    dropzeros!(sparse_struct.jacobians[wind_state_id])

    # recolor if necessary
    # Avoid vector operations for pattern comparison and assignment
    patterns_changed = false
    for i in 1:n_turbines, j in 1:(n_turbines * n_variables)
        if sparse_struct.patterns[i, j, wind_state_id] != sparse_struct.old_patterns[i, j, wind_state_id]
            patterns_changed = true
            break
        end
    end
    if patterns_changed
        recolor_jacobian!(sparse_struct, wind_state_id)
        for i in 1:n_turbines, j in 1:(n_turbines * n_variables)
            sparse_struct.old_patterns[i, j, wind_state_id] = sparse_struct.patterns[i, j, wind_state_id]
        end
    end
    return nothing
end

"""
recolor_jacobian!(sparse_struct::T,wind_state_id)

Helper function that rebuilds the DifferentiationInterface Jacobian preparation for a single
wind state after its sparsity pattern has changed. DI ties a preparation to the exact sparsity
pattern (and coloring) it was built against, so a changed pattern invalidates the old
preparation; unlike SparseDiffTools' `matrix_colors`, there is no separate low-level recoloring
call - rebuilding the `AutoSparse` detector and re-preparing does the recoloring internally.

# Arguments
- `sparse_struct`: sparse_AEP_struct_unstable_pattern struct
- `wind_state_id`: Wind state id
"""
function recolor_jacobian!(sparse_struct::T,wind_state_id) where T <: UnstableSparseMethod
    n_colors = maximum(column_colors(coloring(sparse_struct.jacobians[wind_state_id], ColoringProblem(), sparse_struct.coloring_algorithm)))

    # largest pre-built tier still valid (<=) for the current color count; tier_widths is
    # descending and always ends in 1, so this always finds a match
    tier_idx = findfirst(w -> w <= n_colors, sparse_struct.tier_widths)
    sparse_struct.current_tier[wind_state_id] = tier_idx

    sparse_struct.adtypes[wind_state_id][] = AutoSparse(AutoForwardDiff(chunksize=sparse_struct.tier_widths[tier_idx]);
                    sparsity_detector=KnownJacobianSparsityDetector(sparse_struct.jacobians[wind_state_id]),
                    coloring_algorithm=sparse_struct.coloring_algorithm)
    sparse_struct.preps[wind_state_id][] = nothing
    return nothing
end

"""
calculate_unstable_sparse_jacobian!(sparse_struct::T,x,wind_state_id,prealloc_id,lock)

Helper function that calculates the sparse jacobian for a single wind state using an unstable sparsity pattern

# Arguments
- `sparse_struct`: sparse_AEP_struct_unstable_pattern struct
- `x`: Vector containing the design variables
- `wind_state_id`: Wind state id
- `prealloc_id`: Preallocation id (to select the correct preallocated memory inside the wind farm struct)
- `lock`: SpinLock object to lock the farm struct for multithreadeding
"""
function calculate_unstable_sparse_jacobian!(sparse_struct::T,x,wind_state_id,prealloc_id,lock) where T <: UnstableSparseMethod
    tier_farm = sparse_struct.tier_farms[sparse_struct.current_tier[wind_state_id]]
    adtype = sparse_struct.adtypes[wind_state_id][]
    if isnothing(sparse_struct.preps[wind_state_id][])
        sparse_struct.preps[wind_state_id][] = prepare_jacobian(_wind_state_power_ad!, view(sparse_struct.turbine_powers,:,wind_state_id),
                        adtype, x, Constant(tier_farm), Constant(wind_state_id), Constant(prealloc_id), Constant(lock))
    end

    jacobian!(_wind_state_power_ad!, view(sparse_struct.turbine_powers,:,wind_state_id), sparse_struct.jacobians[wind_state_id],
                    sparse_struct.preps[wind_state_id][], adtype, x, Constant(tier_farm), Constant(wind_state_id), Constant(prealloc_id), Constant(lock))
end

# Sparse spacing constraint methods ########################################################
"""
sparse_spacing_struct

Struct that holds all the necessary variables to calculate the spacing constraints using sparse methods

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `constraint_spacing`: Single float that defines the minimum spacing between turbines
- `constraint_scaling`: Single float that scales the constraint
- `spacing_vec`: Vector containing the spacing constraints
- `jacobian`: Sparse matrix containing the jacobian of the spacing constraints
- `prep`: Ref cell holding the DifferentiationInterface Jacobian preparation (built lazily on first use, `nothing` beforehand)
- `update_function`: Function that updates the spacing struct with the new design variables
- `relevant_list`: 2d array that holds the relevant turbine pairs for the spacing constraint (column 1 holds the first turbine and column 2 holds the second turbine in the pair)
- `ad`: AutoSparse(AutoForwardDiff()) object
- `safe_design_variables`: Vector containing the last set of design variables that satisfy the constraints
- `full_spacing_vec`: Vector containing the full spacing constraints of the farm for final evaluation
"""
struct sparse_spacing_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13} <: AbstractSparseMethod
    turbine_x::T1
    turbine_y::T2
    constraint_spacing::T3 # Single float that defines the minimum spacing between turbines
    constraint_scaling::T4 # Single float that scales the constraint
    spacing_vec::T5 # In place vector
    jacobian::T6
    prep::T7
    update_function::T8
    relevant_list::T9
    ad::T10
    safe_design_variables::T11 # used to hold the last set of design vaiables that satisfy the constraints
    full_spacing_vec::T12
    first_opt::Bool
    value_context::T13 # plain (non-dual) stand-in used to compute the actual spacing_vec value each call
end

"""
spacing_value_context

Lightweight, plain (non-dual) stand-in for `sparse_spacing_struct` holding just the fields
`calculate_spacing!` needs, used to compute the actual (non-differentiated) spacing_vec value.
Kept separate from `sparse_spacing_struct` because that struct's `turbine_x`/`turbine_y` are
permanently dual-typed for differentiation, and mixing a plain-Float64 evaluation into them
would error.
"""
struct spacing_value_context{T1,T2,T3,T4,T5,T6} <: AbstractSparseMethod
    turbine_x::T1
    turbine_y::T2
    update_function::T3
    relevant_list::T4
    constraint_spacing::T5
    constraint_scaling::T6
end

"""
build_sparse_spacing_struct(x,turbine_x,turbine_y,space,scale,update_function;first_opt=true,relevant_spacing_factor=2)

Function that builds a sparse_spacing_struct

# Arguments
- `x`: Vector containing the design variables
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `space`: Single float that defines the minimum spacing between turbines
- `scale`: Single float that scales the constraint
- `update_function`: Function that updates the spacing struct with the new design variables
- `first_opt`: Boolean to determine if this is the first optimization (if true uses no spacing constraints)
- `relevant_spacing_factor`: Single float that defines the factor of space to be considered relevant
- `coloring_algorithm`: SparseMatrixColorings coloring algorithm used for the jacobian coloring (default is `GreedyColoringAlgorithm()`)
"""
function build_sparse_spacing_struct(x,turbine_x,turbine_y,space,scale,update_function;first_opt=true,relevant_spacing_factor=2,
                coloring_algorithm=GreedyColoringAlgorithm())
    if first_opt
        n_constraints = 0
        relevant_list = nothing
        spacing_vec = zeros(eltype(x),n_constraints)
        spacing_jacobian = zeros(eltype(x),n_constraints,length(x))
        prep = Ref{Any}(nothing)
        ad = nothing
        value_context = spacing_value_context(Float64.(turbine_x),Float64.(turbine_y),update_function,relevant_list,space,scale)
        return sparse_spacing_struct(turbine_x,turbine_y,space,scale,spacing_vec,spacing_jacobian,prep,update_function,relevant_list,ad,deepcopy(x),turbine_spacing(turbine_x,turbine_y),first_opt,value_context)
    end

    relevant_list,idx = build_relevant_list(turbine_x,turbine_y,space,relevant_spacing_factor)
    n_constraints = size(relevant_list,1)
    spacing_vec = zeros(eltype(x),n_constraints)

    # calculate jacobian pattern
    s_struct = build_spacing_struct(x,length(turbine_x),space,scale,update_function)
    calculate_spacing_jacobian!(s_struct,x)
    spacing_jacobian = dropzeros(sparse(s_struct.jacobian[idx,:]))

    ad = AutoSparse(AutoForwardDiff(); sparsity_detector=KnownJacobianSparsityDetector(spacing_jacobian),
                    coloring_algorithm=coloring_algorithm)
    prep = Ref{Any}(nothing)
    value_context = spacing_value_context(Float64.(turbine_x),Float64.(turbine_y),update_function,relevant_list,space,scale)

    # pre-type turbine_x/turbine_y to the Dual type that will be produced when x is
    # differentiated through, using the chunk width implied by the jacobian's own coloring
    # (computed directly from the sparsity pattern, no differentiated call needed), tagged to
    # the same top-level function that will actually be differentiated (calculate_spacing!)
    n_colors = maximum(column_colors(coloring(spacing_jacobian, ColoringProblem(), coloring_algorithm)))
    T = dual_type(calculate_spacing!, eltype(x), n_colors)

    turbine_x = Vector{T}(turbine_x)
    turbine_y = Vector{T}(turbine_y)

    return sparse_spacing_struct(turbine_x,turbine_y,space,scale,spacing_vec,spacing_jacobian,prep,update_function,relevant_list,ad,deepcopy(x),turbine_spacing(turbine_x,turbine_y),first_opt,value_context)
end

"""
build_relevant_list(turbine_x,turbine_y,space,factor)

Helper function that builds the relevant list for the sparse spacing constraints

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `space`: Single float that defines the minimum spacing between turbines
- `factor`: Single float that defines the factor of space to be considered relevant
"""
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

"""
build_spacing_struct(x,n_turbines,space,scale,update_function)

Calculates the spacing constraints using sparse methods
"""
function calculate_spacing!(spacing_vec,x,spacing_struct::T) where T<: AbstractSparseMethod
    spacing_struct.update_function(spacing_struct,x)

    sparse_spacing!(spacing_vec,spacing_struct.turbine_x,spacing_struct.turbine_y,spacing_struct.relevant_list)

    spacing_vec .= (spacing_struct.constraint_spacing .- spacing_vec) .* spacing_struct.constraint_scaling
end

"""
sparse_spacing!(spacing_vec,turbine_x,turbine_y,relevant)

Helper function that calculates the relevant spacing constraints

# Arguments
- `spacing_vec`: Vector containing the spacing constraints
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `relevant`: 2d array that holds the relevant turbine pairs for the spacing constraint
"""
function sparse_spacing!(spacing_vec,turbine_x,turbine_y,relevant)
    for i in axes(relevant,1)
        j = relevant[i,1]
        k = relevant[i,2]
        spacing_vec[i] = sqrt((turbine_x[j] - turbine_x[k])^2+(turbine_y[j] - turbine_y[k])^2)
    end
    return spacing_vec
end

"""
calculate_spacing_jacobian!(spacing_struct,x)

Function that calculates the spacing constraints jacobian using sparse methods

# Arguments
- `spacing_struct`: sparse_spacing_struct
- `x`: Vector containing the design variables
"""
function calculate_spacing_jacobian!(spacing_struct::sparse_spacing_struct,x)
    if spacing_struct.first_opt
        return spacing_struct.spacing_vec, spacing_struct.jacobian
    end

    # actual spacing_vec value: cheap plain (non-dual) evaluation against value_context
    calculate_spacing!(spacing_struct.spacing_vec, x, spacing_struct.value_context)

    if isnothing(spacing_struct.prep[])
        spacing_struct.prep[] = prepare_jacobian(calculate_spacing!, spacing_struct.spacing_vec, spacing_struct.ad, x, Constant(spacing_struct))
    end

    jacobian!(calculate_spacing!, spacing_struct.spacing_vec, spacing_struct.jacobian, spacing_struct.prep[], spacing_struct.ad, x, Constant(spacing_struct))

    return spacing_struct.spacing_vec, spacing_struct.jacobian
end

"""
update_safe_design_variables!(spacing_struct::T,x)

Function that updates the safe design variables for the spacing constraints

# Arguments
- `spacing_struct`: sparse_spacing_struct struct
- `x`: Vector containing the design variables
"""
function update_safe_design_variables!(spacing_struct::T,x) where T <: AbstractSparseMethod
    if spacing_struct.first_opt
        spacing_struct.update_function(spacing_struct,x)
        turbine_spacing!(spacing_struct.full_spacing_vec,spacing_struct.turbine_x,spacing_struct.turbine_y)
        min_spacing = minimum(spacing_struct.full_spacing_vec) - spacing_struct.constraint_spacing

        if min_spacing >= 0.0
            spacing_struct.safe_design_variables .= x
        end
    end
    return spacing_struct.safe_design_variables
end

function update_safe_design_variables!(spacing_struct,x)
    return nothing
end

# Sparse boundary constraint methods #######################################################
"""
sparse_boundary_struct

Struct that holds all the necessary variables to calculate the boundary constraints using sparse methods

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `jacobian`: Sparse matrix containing the jacobian of the boundary constraints
- `ad`: AutoSparse(AutoForwardDiff()) object
- `prep`: Ref cell holding the DifferentiationInterface Jacobian preparation (built lazily on first use, `nothing` beforehand)
- `boundary_vec`: Vector containing the boundary constraints
- `boundary_function`: Function that calculates the boundary constraints
- `update_function`: Function that updates the boundary struct with the new design variables
- `boundary_scaling_factor`: Single float that scales the boundary constraint
- `value_context`: plain (non-dual) NamedTuple stand-in used to compute the actual boundary_vec
  value each call (mirrors `turbine_x`/`turbine_y`/`update_function`/`boundary_function`/
  `boundary_scaling_factor`, but with plain Float64 turbine positions instead of dual-typed ones)
"""
struct sparse_boundary_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10} <: AbstractSparseMethod
    turbine_x::T1
    turbine_y::T2
    jacobian::T3
    ad::T4
    prep::T5
    boundary_vec::T6
    boundary_function::T7
    update_function::T8
    boundary_scaling_factor::T9
    value_context::T10
end

"""
calculate_boundary_jacobian!(boundary_struct::sparse_boundary_struct,x)

Function that builds a sparse_boundary_struct

# Arguments
- `boundary_struct`: sparse_boundary_struct
- `x`: Vector containing the design variables
"""
function calculate_boundary_jacobian!(boundary_struct::sparse_boundary_struct,x)
    # actual boundary_vec value: cheap plain (non-dual) evaluation against value_context
    calculate_boundary!(boundary_struct.boundary_vec, x, boundary_struct.value_context)

    if isnothing(boundary_struct.prep[])
        boundary_struct.prep[] = prepare_jacobian(calculate_boundary!, boundary_struct.boundary_vec, boundary_struct.ad, x, Constant(boundary_struct))
    end

    jacobian!(calculate_boundary!, boundary_struct.boundary_vec, boundary_struct.jacobian, boundary_struct.prep[], boundary_struct.ad, x, Constant(boundary_struct))

    return boundary_struct.boundary_vec, boundary_struct.jacobian
end
