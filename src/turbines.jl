abstract type AbstractTurbine end

struct Turbine{TI,TF,TM} <: AbstractTurbine
    id::TI
    rotor_diameter::TF
    hub_height::TF
    ct_model::TM
end
