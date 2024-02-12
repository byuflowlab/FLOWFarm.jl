
abstract type AbstractWindFarmModel end

"""
wind_farm_struct

Unifying struct defining a wind farm

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `turbine_z`: Vector containing z positions of ground the turbines sit on
- `hub_height`: Vector containing hub heights of each turbines as measured form the ground the turbines sit on
- `rotor_diameter`: Vector containing the rotor diameter of each turbine
- `ct_models`: Vector containing ct_models for each turbine
- `generator_efficency`: Vector containing the generator efficiency of each turbine
- `cut_in_speed`: Vector containing the cut in speed of each turbine
- `cut_out_speed`: Vector containing the cut out speed of each turbine
- `rated_speed`: Vector containing the rated speed of each turbine
- `rated_power`: Vector containing the rated power of each turbine
- `wind_resource`: The windresource struct
- `power_models`: Vector containing power models of each turbine
- `model_set`: The models_set struct
- `rotor_sample_points_y`: Vector containing y sample points
- `rotor_sample_points_z`: Vector containing z sample points
- `hours_per_year`: Number of hours in a year
- `objective_scale`: Factor used to scale the objective
- `ideal_AEP`: The ideal AEP of the farm
- `boundary_struct`: Boundary struct
- `spacing_struct`: Spacing struct
- `preallocations`: preallocated space
- `turbine_x_dual`: Dual version of turbine_x
- `turbine_y_dual`: Dual version of turbine_y
- `turbine_z_dual`: Dual version of turbine_z
- `turbine_yaw_dual`: Dual version of turbine_yaw
- `preallocations_dual`: Dual version of preallocations
- `unscale_function`: function that puts the design variables back into SI units
"""
struct wind_farm_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13} <: AbstractWindFarmModel
    # design variables (user access)
    turbine_x::T1
    turbine_y::T2
    hub_height::T3
    turbine_yaw::T4
    rotor_diameter::T5

    # In place container (user access)
    AEP_gradient::T6

    constants::T7

    AEP_scale::T8
    ideal_AEP::T9
    preallocations::T10
    preallocations_dual::T11
    update_function::T12

    duals::T13
end

struct preallocations_struct{V,M}
    prealloc_turbine_velocities::V
    prealloc_turbine_ct::V
    prealloc_turbine_ai::V
    prealloc_turbine_local_ti::V
    prealloc_wake_deficits::M
    prealloc_contribution_matrix::M
    prealloc_deflections::M
    prealloc_sigma_squared::M
end

struct wind_farm_derivative_struct{T1,T2,T3,T4,T5,T6}
    turbine_x_dual::T1
    turbine_y_dual::T2
    hub_height_dual::T3
    turbine_yaw_dual::T4
    rotor_diameter_dual::T5
    forward_cfg::T6
end

struct wind_farm_constants_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12}
    turbine_z::T1
    ct_models::T2
    generator_efficiency::T3
    cut_in_speed::T4
    cut_out_speed::T5
    rated_speed::T6
    rated_power::T7
    wind_resource::T8
    power_models::T9
    model_set::T10
    rotor_sample_points_y::T11
    rotor_sample_points_z::T12
end

######### constraint structs
struct spacing_struct{T1,T2,T3,T4,T5,T6}
    constraint_spacing::T1 # Single float that defines the minimum spacing between turbines in meters
    constraint_scaling::T2 # Single float that scales the constraint
    spacing_vec::T3 # In place vector
    spacing_jacobian::T4
    using_sparsity::T5
    forward_cfg::T6
end

struct boundary_struct{T1,T2,T3,T4,T5,T6}
    boundary_scaling_factor::T1
    boundary_function::T2
    boundary_vec::T3
    boundary_jacobian::T4
    using_sparsity::T5
    forward_cfg::T6
end






############################# outdated
# """
# WindFarm(windfarm, windresource, windfarmstates)

# Struct defining a wind farm

# # Arguments
# - `turbine_x::Array{Float}(Nturbines)`: contains windturbine x coordinates in the global reference frame
# - `turbine_y::Array{Float}(Nturbines)`: contains windturbine y coordinates in the global reference frame
# - `turbine_z::Array{Float}(Nturbines)`: contains windturbine base/z coordinates in the global reference frame
# - `turbine_definition_ids::Array{Int}(Nturbines)`: contains integers for each wind turbine specifying its definition
# - `turbine_definitions::Array{AbstractTurbineDefinition}(Ntypes)`: contains structs defining each wind turbine definition (design) used in the farm
# """
# struct WindFarm{AF1,AF2,AF3,AI,AS} <: AbstractWindFarmModel

#     # farm design properties
#     turbine_x::AF1
#     turbine_y::AF2
#     turbine_z::AF3
#     turbine_definition_ids::AI
#     turbine_definitions::AS

# end


# struct SingleWindFarmState{TI,AF1,AF2,AF3,AF4,AF5,AF6,AF7,AF8,AF9,AI} <: AbstractWindFarmModel

#     # farm properties in rotated frame
#     id::TI
#     turbine_x::AF1
#     turbine_y::AF2
#     turbine_z::AF3
#     turbine_yaw::AF4
#     turbine_ct::AF5
#     turbine_ai::AF6
#     turbine_inflow_velcities::AF7
#     turbine_generators_powers::AF8
#     turbine_local_ti::AF9
#     sorted_turbine_index::AI

# end
