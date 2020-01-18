abstract type AbstractTurbine end

struct Turbine{Coord,TF} <: AbstractTurbine
    coord::Coord
    rotor_diameter::TF
    hub_height::TF
    aI::TF
    yaw::TF
    ct::TF
end

struct Coord{TF}
    x::TF
    y::TF
    z::TF
end
