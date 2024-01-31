"""file with functions used in wind farm optimization employing sparse methods
created January 26, 2024
author: Benjamin Varela
"""

abstract type AbstractOptimizationParametersStruct end
abstract type AbstractSparseOptimizationStruct end

#TODO: documentation and examples

# ∇AEP for Yaw Optimization #######################################################
struct Sparse_yaw_AEP_gradient_struct{T1,T2,T3,T4,T5,T6} <: AbstractSparseOptimizationStruct
    jacobian::T1
    ad::T2
    cache::T3
    function_call::T4
    result_vector::T5
    gradient::T6
end

# TODO: I need to rewrite this to include everything about the farm so I can use one struct with everything and make an easier UI
# TODO: boundary and spacing structs

function create_sparse_yaw_AEP_gradient_struct(yaw,params;threshold=1E-16,unscale_function=x->x)
    results = zeros(length(yaw))
    p(y,x) = turbines_powers_yaw_optimization!(y,x,params;hours_per_year=365.25*24.0,unscale_function=unscale_function)

    # get jacobian
    jacobian = ForwardDiff.jacobian(p,results,yaw)
    jacobian[abs.(jacobian) .<= threshold] .= 0.0
    jacobian = dropzeros(sparse(jacobian))

    # set up SparseDiffTools cache
    sd = JacPrototypeSparsityDetection(; jac_prototype=jacobian)
    ad = AutoSparseForwardDiff()
    cache = sparse_jacobian_cache(ad, sd, p, results, yaw)
    gradient = zeros(length(yaw))

    return Sparse_yaw_AEP_gradient_struct(jacobian,ad,cache,p,results,gradient)
end

function turbines_powers_yaw_optimization!(wt_power,yaw,params;hours_per_year=365.25*24.0,unscale_function=x->x)
    turbine_x = params.turbine_x
    turbine_y = params.turbine_y

    rot_x, rot_y = ff.rotate_to_wind_direction(turbine_x, turbine_y, params.windresource.wind_directions[1])
    sorted_turbine_index = sortperm(rot_x)
    turbine_velocities = ff.turbine_velocities_one_direction(rot_x, rot_y, params.turbine_z, params.rotor_diameter, params.hub_height, unscale_function(yaw),
                            sorted_turbine_index, params.ct_models, params.rotor_points_y, params.rotor_points_z, params.windresource,
                            params.model_set, wind_farm_state_id=1, velocity_only=true)

    wt_power .= ff.turbine_powers_one_direction(params.generator_efficiency, params.cut_in_speed, params.cut_out_speed, params.rated_speed,
                            params.rated_power, params.rotor_diameter, turbine_velocities, x, params.windresource.air_density, params.power_models)

    wt_power .= wt_power .* hours_per_year .* params.windresource.wind_probabilities[1] .* params.obj_scale
end

function calculate_sparse_gradient_yaw!(yaw,gradient_struct)
    sparse_jacobian!(gradient_struct.jacobian, gradient_struct.ad, gradient_struct.cache, gradient_struct.function_call, gradient_struct.result_vector, yaw)
    gradient_struct.gradient .= (sum(gradient_struct.jacobian.jacobian,dims=1))[:]
    return gradient_struct.gradient
end

# ∇AEP for Layout Optimization #######################################################


# Sparse spacing constraint methods #######################################################


# Sparse boundary constraint methods #######################################################
