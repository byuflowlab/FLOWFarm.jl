abstract type AbstractTurbine end

struct Turbine{TI,AF,ACTM, APM} <: AbstractTurbine
    id::TI
    rotor_diameter::AF
    hub_height::AF
    ct_model::ACTM
    power_model::APM
end
