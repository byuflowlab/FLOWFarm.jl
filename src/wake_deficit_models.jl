abstract type AbstractWakeDeficitModel end

"""
    JensenTopHat(alpha)

Container for parameters related to the Jensen Top Hat deficit model

# Arguments
- `alpha::Float`: parameter controlling the wake spreading rate and deficit decay. Default value is 0.1
"""
struct JensenTopHat{TF} <: AbstractWakeDeficitModel
    alpha::TF
end
JensenTopHat() = JensenTopHat(0.2)

"""
    JensenCosine(alpha)

Container for parameters related to the Jensen Cosine deficit model

# Arguments
- `alpha::Float`: parameter controlling the wake deficit decay rate. Default value is 0.1
- `beta::Float`: parameter controlling the width of the cosine function. Default value is 20.0 deg., given in radians.
"""
struct JensenCosine{TF} <: AbstractWakeDeficitModel
    alpha::TF
    beta::TF
end
JensenCosine() = JensenCosine(0.1, 20.0*pi/180.0)
JensenCosine(x) = JensenCosine(x, 20.0*pi/180.0)

"""
    Multizone(me, ke, MU, aU, bU)

Container for parameters related to the Multizone deficit model

# Arguments
- `me::Float`: parameter controlling general wake expansion. Default value is 0.065
- `ke::Array{Float}(3)`: parameters controlling the wake expansion of each zone respectively. Default values are [-0.5 0.22 1.0].
- `MU::Array{Float}(3)`: parameters controlling the wake deficit decay of each zone respectively. Default values are [0.5 1.0 5.5].
- `aU::Float`: parameter impacting the wake deficit decay for a constant wake deflection. Default value is 5.0.
- `bU::Float`: parameter changing the wake deficit decay under yawed conditions. Default value is 1.66.
"""
struct Multizone{ATF, TF} <: AbstractWakeDeficitModel
    me::ATF
    ke::TF
    MU::ATF
    aU::TF
    bU::TF
end
Multizone() = Multizone(0.065, [-0.5 0.22 1.0], [0.5 1.0 5.5], 5.0, 1.66)

"""
    GaussOriginal(k_star)

Container for parameters related to the origina Gaussian deficit model presented by Bastankhah and Porte-Agel 2014

# Arguments
- `k_star::Float`: parameter controlling the wake spreading rate and deficit decay. Default value is 0.075
"""
struct GaussOriginal{TF} <: AbstractWakeDeficitModel
    k_star::TF
end
GaussianOriginal() = GaussOriginal(0.075)

"""
    GaussYaw(turbulence_intensity, horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star)

Container for parameters related to the Gaussian deficit model with yaw presented by Bastankhah and Porte-Agel 2016

# Arguments
- `turbulence_intensity::Float`: the ambient turbulence intensity. No default value is provided.
- `horizontal_spread_rate::Float`: parameter controlling the horizontal spread of the deficit model. Default value is 0.022.
- `vertical_spread_rate::Float`: parameter controlling the vertical spread of the deficit model. Default value is 0.022.
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
"""
struct GaussYaw{TF} <: AbstractWakeDeficitModel
    horizontal_spread_rate::TF
    vertical_spread_rate::TF
    alpha_star::TF
    beta_star::TF
end
GaussYaw() = GaussYaw(0.022, 0.022, 2.32, 0.154)

"""
    wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::JensenTopHat, windfarmstate::SingleWindFarmState)

Computes the wake deficit according to the original Jensen top hat wake model, from the paper:
"A Note on Wind Generator Interaction" by N.O. Jensen (1983)
"""
function wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::JensenTopHat, windfarmstate::SingleWindFarmState)
    # pull out the deflection distances in y (cross stream) and z (up and down)
    deflection_y = deflection[1]
    deflection_z = deflection[2]

    # find delta x, y, and z. dx is the downstream distance from the turbine to
    # the point of interest. dy and dz are the distances from the point of interest
    # and the wake center (in y and z)
    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    dy = loc[2]-(windfarmstate.turbine_y[turbine_id]+deflection_y)
    dz = loc[3]-(windfarmstate.turbine_z[turbine_id]+turbine_definition.hub_height[1]+deflection_z)

    r0 = turbine_definition.rotor_diameter[1]/2.0 #turbine rotor radius
    del = sqrt(dy^2+dz^2) #distance from wake center to the point of interest
    r = model.alpha*dx + r0 #figure (1) from the paper

    if (dx < 0.0) || (del > r)
        loss = 0.0
    else
        loss = 2.0*windfarmstate.turbine_ai[turbine_id]*(r0/(r0+model.alpha*dx))^2 #equation (2) from the paper
    end

    return loss
end

"""
    wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::JensenCosine, windfarmstate::SingleWindFarmState)

Computes the wake deficit according to the original Jensen cosine wake model, from the paper:
"A Note on Wind Generator Interaction" by N.O. Jensen (1983)
"""
function wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::JensenCosine, windfarmstate::SingleWindFarmState)
    """the original Jensen cosine wake model, from the paper: "A Note on Wind
    Generator Interaction" by N.O. Jensen (1983)"""
    # pull out the deflection distances in y (cross stream) and z (up and down)
    deflection_y = deflection[1]
    deflection_z = deflection[2]

    # find delta x, y, and z. dx is the downstream distance from the turbine to
    # the point of interest. dy and dz are the distances from the point of interest
    # and the wake center (in y and z)
    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    dy = loc[2]-(windfarmstate.turbine_y[turbine_id]+deflection_y)
    dz = loc[3]-(windfarmstate.turbine_z[turbine_id]+turbine_definition.hub_height[1]+deflection_z)

    r0 = turbine_definition.rotor_diameter[1]/2.0 #turbine rotor radius
    del = sqrt(dy^2+dz^2) #distance from wake center to the point of interest

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
            loss = 2.0*windfarmstate.turbine_ai[turbine_id]*(ftheta*r0/(r0+model.alpha*dx))^2 #equation (2) from the paper
        end
    end

    return loss
end

"""
    wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::Multizone, windfarmstate::SingleWindFarmState)

Computes the wake deficit at a given location using the original multizone "FLORIS" wake model, from the paper:
"Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study" by Gebraad et al. (2014)
"""
function wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::Multizone, windfarmstate::SingleWindFarmState)

    dt = turbine_definition.rotor_diameter[1]
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
    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    dy = loc[2]-(windfarmstate.turbine_y[turbine_id]+deflection_y)
    dz = loc[3]-(windfarmstate.turbine_z[turbine_id]+turbine_definition.hub_height[1]+deflection_z)

    if dx < 0.
        c = 0.0
    else
        del = sqrt(dy^2+dz^2) #distance from wake center to the point of interest

        # calculate the diameter of the wake in each of the three zones (at the specified dx)
        Dw = zeros(3)
        for i = 1:3
            Dw[i] = max(dt+2*ke*me[i]*dx,0) # equation (13) from the paper
        end
        Rw = Dw./2 # radius of the wake zones

        # calculate the coefficient c of the appropriate wake zone:
        if del > Rw[3]
            c = 0.0
        else
            # equations (15, 16, and 17) from the paper
            if del < Rw[1]
                MUi = MU[1]
            elseif del < Rw[2]
                MUi = MU[2]
            elseif del < Rw[3]
                mUi = MU[3]
            end
            mU = MUi/(cos(aU+bU*windfarmstate.turbine_yaw[turbine_id]))
            c = (dt/(dt+2.0*ke*mU*dx))^2
        end

    end

    loss = 2.0*turbine_definition.aI*c # calculate the wake loss. Equation (14) fro mthe paper

    return loss
end

"""
    wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::GaussOriginal, windfarmstate::SingleWindFarmState)

Computes the wake deficit at a given location using the Gaussian wake model presented by Bastankhah and Porte-Agel in the paper: "A new analytical model for wind-turbine wakes" (2014)
"""
function wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::GaussOriginal, windfarmstate::SingleWindFarmState)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    dy = loc[2]-(windfarmstate.turbine_y[turbine_id]+deflection_y)
    dz = loc[3]-(windfarmstate.turbine_z[turbine_id]+turbine_definition.hub_height[1]+deflection_z)

    # extract turbine properties
    dt = turbine_definition.rotor_diameter[1]
    yaw = windfarmstate.turbine_yaw[turbine_id]
    ct = windfarmstate.turbine_ct[turbine_id]

    as = model.alpha_star
    bs = model.beta_star
    ti = model.turbulence_intensity

    # extract model parameters
    ks = model.k_star       # wake spread rate (k* in 2014 paper)
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate
    as = model.alpha_star
    bs = model.beta_star

    ti = windfarmstate.turbine_local_ti[turbine_id]

    # calculate beta (paper eq: 6)
    beta = 0.5*(1.0+sqrt(1.0-ct))/sqrt(1.0-ct)

    # calculate the length of the potential core (2016 paper eq: 7.3)
    x0 = ((1.0+sqrt(1.0+ct)))/(sqrt(2.0)*as*ti+bs*(1.0-sqrt(1.0-ct)))

    # calculate loss (paper eq: 23)
    enum =((dz/dt)^2+(dy/dt)^2)
    if dx > x0
        denom = (ks*dx/dt+0.2*sqrt(beta))^2
    else
        denom = (ks*x0/dt+0.2*sqrt(beta))^2
    end

    loss = (1.0 - sqrt(1.0-ct/(8.0*denom)))*exp(-enum/(2.0*denom))

end

"""
    _gauss_yaw_potential_core(dt, yaw, ct, as, ti, bs)

Helper function for wake_deficit_model when using the GaussYaw model. Computes the length of the near wake potential core.
"""
function _gauss_yaw_potential_core(d, yaw, ct, as, ti, bs)
    # from Bastankhah and Porte-Agel 2016 eqn 7.3

    x0 = d*(cos(yaw)*(1.0+sqrt(1.0-ct)))/(sqrt(2.0)*(as*ti+bs*(1.0-sqrt(1.0-ct))))

    return x0
end

"""
    _gauss_yaw_spread(dt, k, dx, x0, yaw)

Helper function for wake_deficit_model when using the GaussYaw model. Computes the standard deviation of the wake.
"""
function _gauss_yaw_spread(dt, k, dx, x0, yaw)
    # from Bastankhah and Porte-Agel 2016 eqn 7.2

    sigma = dt*(k*(dx-x0)/dt+cos(yaw)/sqrt(8.0))

    return sigma

end

"""
    wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::GaussYaw, windfarmstate::SingleWindFarmState)

Computes the wake deficit at a given location using the The Gaussian wake model presented by Bastankhah and Porte-Agel in the paper: "Experimental and theoretical study of wind turbine wakes in yawed conditions" (2016)
"""
function wake_deficit_model(loc, deflection, turbine_id, turbine_definition::TurbineDefinition, model::GaussYaw, windfarmstate::SingleWindFarmState)

    deflection_y = deflection[1]
    deflection_z = deflection[2]

    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    dy = loc[2]-(windfarmstate.turbine_y[turbine_id]+deflection_y)
    dz = loc[3]-(windfarmstate.turbine_z[turbine_id]+turbine_definition.hub_height[1]+deflection_z)

    # extract turbine properties
    dt = turbine_definition.rotor_diameter[1]
    yaw = windfarmstate.turbine_yaw[turbine_id]
    ct = windfarmstate.turbine_ct[turbine_id]

    # extract model parameters
    # ks = model.k_star       # wake spread rate (k* in 2014 paper)
    ti = windfarmstate.turbine_local_ti[turbine_id]
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate
    as = model.alpha_star
    bs = model.beta_star

    if dx > 0.0 # loss in the wake

        # calculate the length of the potential core (paper eq: 7.3)
        x0 = _gauss_yaw_potential_core(dt, yaw, ct, as, ti, bs)

        # calculate wake spread
        if dx > x0
            # calculate horizontal wake spread (paper eq: 7.2)
            sigma_y = _gauss_yaw_spread(dt, ky, dx, x0, yaw)

            # calculate vertical wake spread (paper eq: 7.2)
            sigma_z = _gauss_yaw_spread(dt, kz, dx, x0, 0.0)

        else # linear interpolation in the near wakes
            # calculate horizontal wake spread
            sigma_y = _gauss_yaw_spread(dt, ky, x0, x0, yaw)

            # calculate vertical wake spread
            sigma_z = _gauss_yaw_spread(dt, kz, x0, x0, 0.0)
        end

        # calculate velocity deficit
        ey = exp(-0.5*(dy/sigma_y)^2.0)
        ez = exp(-0.5*(dz/sigma_z)^2.0)

        loss = (1.0-sqrt(1.0-ct*cos(yaw)/(8.0*(sigma_y*sigma_z/dt^2.0))))*ey*ez

    else # loss upstream of turbine
        loss = 0.0

    end

end
