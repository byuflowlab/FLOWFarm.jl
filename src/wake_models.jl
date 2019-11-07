abstract type AbstractWakeModel end

struct Jensen <: AbstractWakeModel
    alpha
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

function wake_model(loc, deflection, model::Jensen, turbine::Turbine)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    r0 = turbine.rotor_diameter/2.0
    del = sqrt(dy^2+dz^2)
    r = model.alpha*dx + r0
    if del > r
        loss = 0.0
    else
        loss = 2.0*turbine.aI*(r0/(r0+model.alpha*dx))^2
    end
end


function wake_model(loc, model::Jensen, turbine::Turbine)
    dx = loc[1]-turbine.coord.x
    dy = loc[2]-turbine.coord.y
    dz = loc[3]-(turbine.coord.z+turbine.hub_height)
    r0 = turbine.rotor_diameter/2.0
    del = sqrt(dy^2+dz^2)
    r = model.alpha*dx + r0
    if del > r
        loss = 0.0
    else
        loss = 2.0*turbine.aI*(r0/(r0+model.alpha*dx))^2
    end
end


function loss(x_locations, y_locations, z_locations, turbine::Turbine,
              deflection_field, flow_field, model::Multizone)
    #return the losses from the multizone FLORIS wake model
end

function loss(x_locations, y_locations, z_locations, turbine::Turbine,
              deflection_field, flow_field, model::Gauss)
    #return the losses using the Bastankhah Gaussian wake model
end
