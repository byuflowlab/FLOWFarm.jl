
abstract type AbstractWindFarmModel end
abstract type AbstractWindFarmProblem end
abstract type AbstractModelSet end

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

"""
    WindFarmProblemDescription(windfarm, windresource, windfarmstates)

Container for objects defining a wind farm problem

# Arguments
- `wind_farm::Array{WindFarm}(1)`: contains windturbine coordinates and definitions
- `wind_resource::Array{AbstracWindResource}(1)`: wind resource description
- `wind_farm_states::Array{SingleWindFarmState}(Nstates)`: contains turbine coordinates operational states
"""
struct WindFarmProblemDescription{AFM,AWR,AFS} <: AbstractWindFarmProblem
   
    wind_farm::AFM
    wind_resource::AWR
    wind_farm_states::AFS

end


"""
    WindFarmModelSet(wakedeficitmodel, wake_deflection_model, wake_combination_model, ti_model, wind_shear_model)

Container for objects defining models to use in wind farm calculations

# Arguments
- `wake_defiict_model::Array{AbstractWakeDeficitModel}(1)`: contains a struct defining the desired wake deficit model
- `wake_deflection_model::Array{AbstractWakeDeflectionModel}(1)`: contains a struct defining the desired wake deflection model
- `wake_combination_model::Array{AbstractWakeCombinationModel}(1)`: contains a struct defining the desired wake combination model
- `ti_model::Array{AbstractWakeDeflectionModel}(1)`: contains a struct defining the desired wake deflection model
- `wind_shear_model::Array{SingleWindFarmState}(Nstates)`: contains turbine coordinates operational states
"""
struct WindFarmModelSet{ADTM,ADNM,ACM,ATIM,ASM} <: AbstractModelSet



end
