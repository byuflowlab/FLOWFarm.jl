
abstract type AbstractWindFarmModel end

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
struct WindFarm{AF,AI,AS} <: AbstractWindFarmModel
    
    # farm design properties
    turbine_x::AF
    turbine_y::AF
    turbine_z::AF
    turbine_definition_ids::AI
    turbine_definitions::AS

end


struct SingleWindFarmState{TI,AF,AI} <: AbstractWindFarmModel

    # farm properties in rotated frame
    id::TI
    turbine_x::AF
    turbine_y::AF
    turbine_z::AF
    turbine_yaw::AF
    turbine_ct::AF
    turbine_ai::AF
    sorted_turbine_index::AI
    turbine_inflow_velcities::AF
    turbine_generators_powers::AF

end

