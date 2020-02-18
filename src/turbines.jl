abstract type AbstractTurbine end

struct Turbine{TI,TF} <: AbstractTurbine
    id::TI
    rotor_diameter::TF
    hub_height::TF
end
