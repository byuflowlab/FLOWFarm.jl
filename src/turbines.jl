abstract type AbstractTurbine end

struct Turbine <: AbstractTurbine
    coord
    rotor_diameter
    hub_height
    aI
    gamma
end

struct Coord
    x
    y
    z
end
