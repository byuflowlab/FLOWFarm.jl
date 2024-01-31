
abstract type AbstractWindFarmModel end

"""
wind_farm_struct

Struct defining a wind farm

# Arguments
- `turbine_x`: Vector containing
"""
struct wind_farm_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21,T22,T23} <: AbstractWindFarmModel
    turbine_x::T1
    turbine_y::T2
    turbine_z::T3
    rotor_diameter::T4
    hub_height::T5
    turbine_yaw::T6
    ct_models::T7
    generator_efficency::T8
    cut_in_speed::T9
    cut_out_speed::T10
    rated_speed::T11
    rated_power::T12
    wind_resource::T13
    power_models::T14
    model_set::T15
    rotor_sample_points_y::T16
    rotor_sample_points_z::T17
    hours_per_year::T18
    objective_scale::T19
    ideal_AEP::T20
    boundary_struct::T21
    spacing_struct::T22
    struct_update_function::T23
end

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

# # abstract type AbstractBoundary end

# # """
# # CircleBoundary(boundary_center, boundary_radius)

# # # Arguments
# # - `boundary_center::Float`: center of wind farm boundary
# # - 'boundary_radius::Float': radius of wind farm boundary
# # """
# # struct CircleBoundary{TF,TF} <: AbstractBoundary
# #     boundary_center::TF
# #     boundary_radius::TF
# # end

# # """
# # PolygonBoundary(boundary_center, boundary_radius)

# # # Arguments
# # - `boundary_center::Float`: center of wind farm boundary
# # - 'boundary_radius::Float': radius of wind farm boundary
# # """
# # struct PolygonBoundary{ATF} <: AbstractBoundary
# #     boundary_vertices::ATF
# #     boundary_normals::ATF
# # end

# function initialize_polygon_boundary(vertices)
#     normals = boundary_normals_calculator(vertices)
#     boundary = PolygonBoundary(vertices, normals)
#     return boundary
# end
