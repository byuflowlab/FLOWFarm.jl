abstract type AbstractTurbineDefinition end

struct TurbineDefinition{TI,AF,TF,CTM,PM} <: AbstractTurbineDefinition
    id::TI
    rotor_diameter::AF
    hub_height::AF
    cut_in_speed::TF
    rated_speed::TF
    cut_out_speed::TF
    rated_power::TF
    generator_efficiency::TF
    ct_model::CTM
    power_model::PM
end
