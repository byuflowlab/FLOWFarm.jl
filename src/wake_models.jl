include("combination_models.jl")
include("turbines.jl")

abstract type AbstractWakeModel end

struct JensenTopHat{TF} <: AbstractWakeModel
    alpha::TF
end

struct Multizone{ATF, TF} <: AbstractWakeModel
    me::ATF
    ke::TF
    MU::ATF
    aU::TF
    bU::TF
end

struct Gauss <: AbstractWakeModel
    version
    # for version==2014
    k_star
    # for version==2016
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


function wake_model(loc, deflection, model::Gauss, turbine::Turbine)

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

    if model.version == 2014
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

    elseif model.version == 2016
      """The Gaussian wake model presented by Bastankhah and Porte-Agel in
      the paper: "Experimental and theoretical study of wind turbine wakes
      in yawed conditions" (2016)"""

      # calculate the length of the potential core (paper eq: 7.3)
      x0 = (cos(yaw)*(1+sqrt(1+ct)))/(sqrt(2)*as*TI+bs*(1-sqrt(1-ct)))

      # calculate wake spread
      if dx > x0

        # calculate horizontal wake spread (paper eq: 7.2)
        sigma_y = Dt*(ky*dx/Dt+cos(yaw)/sqrt(8))

        # calculate vertical wake spread (paper eq: 7.2)
        sigma_z = Dt*(kz*dx/Dt+1/sqrt(8))

      else # linear interpolation in the near wakes

        # calculate horizontal wake spread
        sigma_y = Dt*(ky*x0/Dt+cos(yaw)/sqrt(8))

        # calculate vertical wake spread
        sigma_z = Dt*(kz*x0/Dt+1/sqrt(8))

      end

      # calculate velocity deficit
      ey = exp(-0.5*(dy/sigma_y)^2)
      ez = exp(-0.5*(dz/sigma_z)^2)

      loss = (1-sqrt(1-ct*cos(yaw)/(8*(sigma_y*sigma_z/Dt^2))))*ey*ez

    end
end
