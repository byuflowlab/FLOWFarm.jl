abstract type AbstractTurbineDefinition end

struct TurbineDefinition{TI,AF,ACTM, APM} <: AbstractTurbineDefinition
    id::TI
    rotor_diameter::AF
    hub_height::AF
    ct_model::ACTM
    power_model::APM
end
