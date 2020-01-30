include("combination_models.jl")
include("turbines.jl")

abstract type AbstractWakeModel end

struct JensenTopHat{TF} <: AbstractWakeModel
    alpha::TF
end

struct JensenCosine{TF} <: AbstractWakeModel
    alpha::TF
    beta::TF
end

struct Multizone{ATF, TF} <: AbstractWakeModel
    me::ATF
    ke::TF
    MU::ATF
    aU::TF
    bU::TF
end

struct GaussOriginal <: AbstractWakeModel
    k_star
end

struct GaussYaw <: AbstractWakeModel
    turbulence_intensity
    horizontal_spread_rate
    vertical_spread_rate
    alpha_star
    beta_star
end

function wake_model(loc, deflection, model::JensenTopHat, turbine::Turbine)
    """the original Jensen top hat wake model, from the paper: "A Note on Wind
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

function wake_model(loc, deflection, model::JensenCosine, turbine::Turbine)
    """the original Jensen cosine wake model, from the paper: "A Note on Wind
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
        loss = 0.0 # no loss outside the wake
    else
        d = r0/tan(model.beta) # distance from fulcrum of wake cone to wind turbine hub
        theta = atan(dy/(dx+d)) # angle from center of wake to point of interest
        if theta > model.beta # if you're outside the wake
            loss = 0.0
        else # if you're inside the wake
            n = pi / model.beta # see Jensen 1983 eq. 3. Value used for n in paper was 9, corresponding to beta = 20.0 deg.
            ftheta = (1.0 + cos(n * theta)) / 2.0 # cosine term to be applied to loss equation as per Jensen 1983
            loss = 2.0*turbine.aI*(ftheta*r0/(r0+model.alpha*dx))^2 #equation (2) from the paper
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
            mU = MU[1]/(cos(aU+bU*turbine.yaw))
            c = (Dt/(Dt+2.0*ke*mU*dx))^2
        elseif del < Rw[2]
            mU = MU[2]/(cos(aU+bU*turbine.yaw))
            c = (Dt/(Dt+2.0*ke*mU*dx))^2
        elseif del < Rw[3]
            mU = MU[3]/(cos(aU+bU*turbine.yaw))
            c = (Dt/(Dt+2.0*ke*mU*dx))^2
        else
            c = 0.0
        end
    end
    loss = 2.0*turbine.aI*c # calculate the wake loss. Equation (14) fro mthe paper
    return loss
end


function wake_model(loc, deflection, model::GaussOriginal, turbine::Turbine)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    # extract turbine properties
    Dt = turbine.rotor_diameter
    yaw = turbine.yaw
    ct = turbine.ct

    as = model.alpha_star
    bs = model.beta_star
    TI = model.turbulence_intensity

    # extract model parameters
    ks = model.k_star       # wake spread rate (k* in 2014 paper)
    TI = model.turbulence_intensity
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate
    as = model.alpha_star
    bs = model.beta_star

    """The Gaussian wake model presented by Bastankhah and Porte-Agel in
    the paper: "A new analytical model for wind-turbine wakes" (2014)"""

    # calculate beta (paper eq: 6)
    beta = 0.5*(1.0+sqrt(1.0-ct))/sqrt(1.0-ct)

    # calculate the length of the potential core (2016 paper eq: 7.3)
    x0 = ((1.0+sqrt(1.0+ct)))/(sqrt(2.0)*as*TI+bs*(1.0-sqrt(1.0-ct)))

    # calculate loss (paper eq: 23)
    enum =((dz/Dt)^2+(dy/Dt)^2)
    if dx > x0
    denom = (ks*dx/Dt+0.2*sqrt(beta))^2
    else
    denom = (ks*x0/Dt+0.2*sqrt(beta))^2
    end
    loss = (1.0 - sqrt(1.0-ct/(8.0*denom)))*exp(-enum/(2.0*denom))

end

function wake_model(loc, deflection, model::GaussYaw, turbine::Turbine)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    # extract turbine properties
    Dt = turbine.rotor_diameter
    yaw = turbine.yaw
    ct = turbine.ct

    # extract model parameters
    # ks = model.k_star       # wake spread rate (k* in 2014 paper)
    TI = model.turbulence_intensity
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate
    as = model.alpha_star
    bs = model.beta_star

    """The Gaussian wake model presented by Bastankhah and Porte-Agel in
    the paper: "Experimental and theoretical study of wind turbine wakes
    in yawed conditions" (2016)"""

    if dx > 0.0 # loss in the wake

        # calculate the length of the potential core (paper eq: 7.3)
        x0 = Dt*(cos(yaw)*(1.0+sqrt(1.0-ct)))/(sqrt(2.0)*(as*TI+bs*(1.0-sqrt(1.0-ct))))

        # calculate wake spread
        if dx > x0
            # calculate horizontal wake spread (paper eq: 7.2)
            sigma_y = Dt*(ky*(dx-x0)/Dt+cos(yaw)/sqrt(8.0))

            # calculate vertical wake spread (paper eq: 7.2)
            sigma_z = Dt*(kz*(dx-x0)/Dt+1.0/sqrt(8.0))

        else # linear interpolation in the near wakes
            # calculate horizontal wake spread
            sigma_y = Dt*(ky*(x0-x0)/Dt+cos(yaw)/sqrt(8.0))

            # calculate vertical wake spread
            sigma_z = Dt*(kz*(x0-x0)/Dt+1.0/sqrt(8.0))
        end 

        # calculate velocity deficit
        ey = exp(-0.5*(dy/sigma_y)^2.0)
        ez = exp(-0.5*(dz/sigma_z)^2.0)

        loss = (1.0-sqrt(1.0-ct*cos(yaw)/(8.0*(sigma_y*sigma_z/Dt^2.0))))*ey*ez

    else # loss upstream of turbine
        loss = 0.0

    end


end

