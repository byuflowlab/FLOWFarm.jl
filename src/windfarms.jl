
abstract type AbstractWindFarmModel end

struct WindFarm{AF,AS} <: AbstractWindFarmModel
    
    # farm design properties
    turbine_x::AF
    turbine_y::AF
    turbine_z::AF
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

