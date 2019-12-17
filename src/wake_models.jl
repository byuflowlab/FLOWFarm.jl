include("combination_models.jl")
include("turbines.jl")

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
    """the original jensen wake model, from the paper: "A Note on Wind
    Generator Interaction" by N.O. Jensen (1983)"""
    # pull out the deflection distances in y (cross stream) and z (up and down)
    deflection_y = deflection[1]
    deflection_z = deflection[2]

    # find delta x, y, and z. dx is the downstream distance from the turbine to
    # the point of interest. dy and dz are the distances from the point of interest
    # and the wake center (in y and z)
    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    r0 = turbine.rotor_diameter/2.0 #turbine rotor radius
    del = sqrt(dy^2+dz^2) #distance from wake center to the point of interest
    r = model.alpha*dx + r0 #figure (1) from the paper

    if dx < 0.
        loss = 0.0
    else
        if del > r #if you're outside the wake
            loss = 0.0
        else # if you're inside the wake
            loss = 2.0*turbine.aI*(r0/(r0+model.alpha*dx))^2 #equation (2) from the paper
        end
    end
    return loss
end


function wake_model(loc, deflection, model::Multizone, turbine::Turbine)
    """The original multizone "FLORIS" wake model, from the paper:
    "Wind plant power optimization through yaw control using a parametric model
    for wake effects—a CFD simulation study" by Gebraad et al. (2014)"""

    Dt = turbine.rotor_diameter
    # extract model parameters
    ke = model.ke
    me = model.me
    MU = model.MU
    aU = model.aU
    bU = model.bU
    # pull out the deflection distances in y (cross stream) and z (up and down)
    deflection_y = deflection[1]
    deflection_z = deflection[2]

    # find delta x, y, and z. dx is the downstream distance from the turbine to
    # the point of interest. dy and dz are the distances from the point of interest
    # and the wake center (in y and z)
    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    if dx < 0.
        c = 0.0
    else
        del = sqrt(dy^2+dz^2) #distance from wake center to the point of interest

        # calculate the diameter of the wake in each of the three zones (at the specified dx)
        Dw = zeros(3)
        for i = 1:3
            Dw[i] = max(Dt+2*ke*me[i]*dx,0) # equation (13) from the paper
        end
        Rw = Dw./2 # radius of the wake zones
        # calculate the coefficient c of the appropriate wake zone:
        # equations (15, 16, and 17) from the paper
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
    end
    loss = 2.0*turbine.aI*c # calculate the wake loss. Equation (14) fro mthe paper
    return loss
end


function wake_model(loc, deflection, model::Gauss, turbine::Turbine)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    Dt = turbine.rotor_diameter

    """This is just a simple Gaussian model, needs to be replaced with the full
    Bastankhah"""
    k = 0.0325
    CT = 8.0/9.0
    sigma = k*dx + Dt/sqrt(8.0)
    loss = (1.0-sqrt(1.0-CT/(8*sigma^2/(Dt^2))))*exp(-1.0/2.0*(dy/sigma)^2)
    return loss
end


"""total loss from every turbine at a single point"""

function wake_model(loc, deflection, model::AbstractWakeModel, farm::Array{turbine::Turbine,1})
    loss_vec = zeros(length(deflection))
    for (i,turbine) in enumerate(farm)
        loss_vec[i] = wake_model(loc, deflection[i], model, turbine)
    end
    loss = combination_model(loss_vec)
    return loss
end
