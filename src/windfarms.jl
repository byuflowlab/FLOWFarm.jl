
abstract type AbstractWindFarmModel end

"""
wind_farm_struct

Unifying struct defining a wind farm and all necessary variables to calculate the AEP

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `hub_height`: Vector containing hub heights of each turbines as measured form the ground the turbines sit on
- `turbine_yaw`: Vector containing yaw angle of each turbine in radians
- `rotor_diameter`: Vector containing the rotor diameter of each turbine
- `results`: DiffResults object to extract AEP when calculating AEP gradient
- `constants`: wind_farm_constants_struct
- `AEP_scale`: Scaling factor for the AEP
- `ideal_AEP`: The ideal AEP of the farm
- `preallocations`: preallocated space
- `update_function`: function that takes the design variables x and updates the farm struct
- `AEP_gradient`: The gradient of the AEP
- `AEP`: The AEP of the farm
- `config`: The ForwardDiff config object if using ForwardDiff for AEP gradient calculation, otherwise nothing
"""
struct wind_farm_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14} <: AbstractWindFarmModel
    turbine_x::T1
    turbine_y::T2
    hub_height::T3
    turbine_yaw::T4
    rotor_diameter::T5
    results::T6
    constants::T7
    AEP_scale::T8
    ideal_AEP::T9
    preallocations::T10
    update_function::T11
    AEP_gradient::T12
    AEP::T13
    config::T14
end

"""
preallocations_struct

struct that holds all the preallocated space for AEP calculation with one per thread used

# Arguments
- `prealloc_turbine_velocities`: Vector containing preallocated space for turbine velocities
- `prealloc_turbine_ct`: Vector containing preallocated space for turbine ct
- `prealloc_turbine_ai`: Vector containing preallocated space for turbine ai
- `prealloc_turbine_local_ti`: Vector containing preallocated space for turbine local ti
- `prealloc_wake_deficits`: Matrix containing preallocated space for wake deficits
- `prealloc_contribution_matrix`: Matrix containing preallocated space for contribution matrix
- `prealloc_deflections`: Matrix containing preallocated space for deflections
- `prealloc_sigma_squared`: Matrix containing preallocated space for sigma squared
"""
struct preallocations_struct{V,M,VI,Ve}
    prealloc_turbine_velocities::V
    prealloc_turbine_ct::V
    prealloc_turbine_ai::V
    prealloc_turbine_local_ti::V
    prealloc_wake_deficits::M
    prealloc_contribution_matrix::M
    prealloc_deflections::M
    prealloc_sigma_squared::M
    prealloc_rot_x::V
    prealloc_rot_y::V
    prealloc_power::V
    prealloc_sort_index::VI
    prealloc_diam::V
    prealloc_hub::V
    prealloc_yaw::V
    prealloc_state_aep::Ve
end

"""
wind_farm_constants_struct

struct that holds all the constants for the wind farm

# Arguments
- `turbine_z`: Vector containing z positions of ground the turbines sit on
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
"""
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
"""
spacing_struct

Struct defining the spacing constraints

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `constraint_spacing`: Single float that defines the minimum spacing between turbines in meters
- `constraint_scaling`: Single float that scales the constraint
- `spacing_vec`: Vector containing the spacing constraints
- `jacobian`: Matrix containing the jacobian of the spacing constraints
- `config`: The ForwardDiff config object if using ForwardDiff for jacboian calculation
- `update_function`: function that takes the design variables x and updates the spacing struct
"""
struct spacing_struct{T1,T2,T3,T4,T5,T6,T7,T8}
    turbine_x::T1
    turbine_y::T2
    constraint_spacing::T3 # Single float that defines the minimum spacing between turbines in meters
    constraint_scaling::T4 # Single float that scales the constraint
    spacing_vec::T5 # In place vector
    jacobian::T6
    config::T7
    update_function::T8
end

"""
boundary_struct

Struct defining the boundary constraints

# Arguments
- `turbine_x`: Vector containing x positions of turbines
- `turbine_y`: Vector containing y positions of turbines
- `boundary_scaling_factor`: Single float that scales the constraint
- `boundary_function`: function that takes the boundary vector and the design variables to update the boundary vector
- `boundary_vec`: Vector containing the boundary constraints
- `jacobian`: Matrix containing the jacobian of the boundary constraints
- `config`: The ForwardDiff config object if using ForwardDiff for jacboian calculation
- `update_function`: function that takes the design variables x and updates the boundary struct
"""
struct boundary_struct{T1,T2,T3,T4,T5,T6,T7,T8}
    turbine_x::T1
    turbine_y::T2
    boundary_scaling_factor::T3
    boundary_function::T4
    boundary_vec::T5
    jacobian::T6
    config::T7
    update_function::T8
end






############################# outdated
"""
WindFarm(windfarm, windresource, windfarmstates)

Struct defining a wind farm

# Arguments
- `turbine_x::Array{Float}(Nturbines)`: contains windturbine x coordinates in the global reference frame
- `turbine_y::Array{Float}(Nturbines)`: contains windturbine y coordinates in the global reference frame
- `turbine_z::Array{Float}(Nturbines)`: contains windturbine base/z coordinates in the global reference frame
- `turbine_definition_ids::Array{Int}(Nturbines)`: contains integers for each wind turbine specifying its definition
- `turbine_definitions::Array{AbstractTurbineDefinition}(Ntypes)`: contains structs defining each wind turbine definition (design) used in the farm
"""
struct WindFarm{AF1,AF2,AF3,AI,AS} <: AbstractWindFarmModel

    # farm design properties
    turbine_x::AF1
    turbine_y::AF2
    turbine_z::AF3
    turbine_definition_ids::AI
    turbine_definitions::AS

end


struct SingleWindFarmState{TI,AF1,AF2,AF3,AF4,AF5,AF6,AF7,AF8,AF9,AI} <: AbstractWindFarmModel

    # farm properties in rotated frame
    id::TI
    turbine_x::AF1
    turbine_y::AF2
    turbine_z::AF3
    turbine_yaw::AF4
    turbine_ct::AF5
    turbine_ai::AF6
    turbine_inflow_velcities::AF7
    turbine_generators_powers::AF8
    turbine_local_ti::AF9
    sorted_turbine_index::AI

end
