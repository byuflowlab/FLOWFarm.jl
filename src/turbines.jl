abstract type AbstractTurbineDefinition end

# not used

struct TurbineDefinition{TI,AF,CTM,PM} <: AbstractTurbineDefinition
    id::TI
    rotor_diameter::AF
    hub_height::AF
    cut_in_speed::AF
    rated_speed::AF
    cut_out_speed::AF
    rated_power::AF
    generator_efficiency::AF
    ct_model::CTM
    power_model::PM
end
