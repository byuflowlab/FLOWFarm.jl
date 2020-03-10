include("turbines.jl")
include("wind_resources.jl")

abstract type AbstractPowerModel end

struct ConstantCp{TF} <: AbstractPowerModel
    cp::TF
    generator_efficiency::TF
end

struct PowerCurveCubic{TF} <: AbstractPowerModel
    cut_in_speed::TF
    cut_out_speed::TF
    rated_speed::TF
    rated_power::TF
end

function calculate_turbine_power(turbine::AbstractTurbine, farmstate::SingleWindFarmState, wind_model::AbstractWindResourceModel, power_model::ConstantCp)
    
    id = turbine.id
    cp = power_model.cp[1]
    generator_efficiency = power_model.generator_efficiency[1]
    air_density = wind_model.air_density
    rotor_diameter = turbine.rotor_diameter[1]
    rotor_area = pi*(rotor_diameter^2)/4.0
    wt_velocity = farmstate.turbine_inflow_velcities[turbine.id[1]]

    wt_power = generator_efficiency*(0.5*air_density*rotor_area*cp*wt_velocity^3)

    return wt_power

end

function calculate_aep(windfarm::AbstractWindFarmModel)

end