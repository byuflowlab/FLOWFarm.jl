
abstract type AbstractWindFarmModel end

"""
WindFarm(windfarm, windresource, windfarmstates)

Struct defining a wind farm

# Arguments
- `adding another commit`
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



