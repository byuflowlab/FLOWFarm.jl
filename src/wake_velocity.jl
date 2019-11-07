include("turbines.jl")

abstract type AbstractWakeModel end

struct Jensen <: AbstractWakeModel
    we
end

struct Multizone <: AbstractWakeModel
    me
    we
    aU
    bU
    mU
end

struct Gauss <: AbstractWakeModel
    ka
    kb
    alpha
    beta
    ad
    bd
end

function loss(x_locations, y_locations, z_locations, turbine::Turbine,
              deflection_field, flow_field, model::Jensen)

    m = model.we
    x = x_locations .- turbine.coord.x1
    b = turbine.rotor_diameter

    boundary_line = m .* x .+ b

    y_upper = boundary_line .+ turbine.coord.x1 .+ deflection_field
    y_lower = -1 .* boundary_line .+ turbine.coord.x2 .+ deflection_field

    z_upper = boundary_line .+ turbine.hub_height
    z_lower = -1 .* boundary_line .+ turbine.hub_height

    # calculate the wake velocity
    c = (turbine.rotor_diameter ./
        (2. * model.we .* (x_locations .- turbine.coord.x1) .+ turbine.rotor_diameter)).^2

    # filter points upstream and beyond the upper and lower bounds of the wake
    # c[x_locations .- turbine.coord.x1 < 0] = 0
    # c[y_locations > y_upper] = 0
    # c[y_locations < y_lower] = 0
    # c[z_locations > z_upper] = 0
    # c[z_locations < z_lower] = 0

    return 2.0 .* turbine.aI .* c .* flow_field.u_initial, zeros(size(flow_field.u_initial)), zeros(size(flow_field.u_initial))
end


function loss(x_locations, y_locations, z_locations, turbine::Turbine,
              deflection_field, flow_field, model::Multizone)
    #return the losses from the multizone FLORIS wake model
end

function loss(x_locations, y_locations, z_locations, turbine::Turbine,
              deflection_field, flow_field, model::Gauss)
    #return the losses using the Bastankhah Gaussian wake model
end

function mult(x, y, z)
    println(x*y*z)
end

function sum2(x, y)
    println(x + y)
end
