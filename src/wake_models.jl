abstract type AbstractWakeModel end

struct Jensen <: AbstractWakeModel
    alpha
end

struct Multizone <: AbstractWakeModel
    me
    ke
    MU
    aU
    bU
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


function wake_model(loc, deflection, model::Multizone, turbine::Turbine)

    Dt = turbine.rotor_diameter
    ke = model.ke
    me = model.me
    MU = model.MU
    aU = model.aU
    bU = model.bU
    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)
    del = sqrt(dy^2+dz^2)

    Dw = zeros(3)
    for i = 1:3
        Dw[i] = max(Dt+2*ke*me[i]*(loc[1]-turbine.coord.x),0)
    end

    Rw = Dw./2

    if del < Rw[1]
        mU = MU[1]/(cos(aU+bU*turbine.gamma))
        c = (Dt/(Dt+2.0*ke*mU*dx))^2
    elseif del < Rw[2]
        mU = MU[2]/(cos(aU+bU*turbine.gamma))
        c = (Dt/(Dt+2.0*ke*mU*dx))^2
    elseif del < Rw[3]
        mU = MU[3]/(cos(aU+bU*turbine.gamma))
        c = (Dt/(Dt+2.0*ke*mU*dx))^2
    else
        c = 0.0
    end
    loss = 2.0*turbine.aI*c
end


function wake_model(loc, deflection, model::Gauss, turbine::Turbine)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

end
