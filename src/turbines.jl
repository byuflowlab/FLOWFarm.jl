abstract type AbstractTurbine end

struct Turbine{AI,AF,AM} <: AbstractTurbine
    id::AI
    rotor_diameter::AF
    hub_height::AF
    ct_model::AM
end
