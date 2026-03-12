export NoYawDeflection, GaussYawDeflection, GaussYawVariableSpreadDeflection, JiminezYawDeflection, MultizoneDeflection
abstract type AbstractWakeDeflectionModel end

"""
    NoYawDeflection()

Allows for bypassing deflection calculations.
"""
struct NoYawDeflection <: AbstractWakeDeflectionModel
end

"""
    GaussYawDeflection(horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star, interpolation)

Container for parameters related to the Gaussian deflection model presented by Bastankhah and Porte-Agel 2016

# Arguments
- `horizontal_spread_rate::Float`: parameter controlling the horizontal spread of the deficit model. Default value is 0.022.
- `vertical_spread_rate::Float`: parameter controlling the vertical spread of the deficit model. Default value is 0.022.
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
- `interpolation::Bool`: boolean stating if the the near wake should be interpolated. Default value is true.
"""
struct GaussYawDeflection{TF1, TF2, TF3, TF4, BO} <: AbstractWakeDeflectionModel
    horizontal_spread_rate::TF1
    vertical_spread_rate::TF2
    alpha_star::TF3
    beta_star::TF4
    interpolate_sigma::BO
end
GaussYawDeflection() = GaussYawDeflection(0.022, 0.022, 2.32, 0.154, true)
GaussYawDeflection(interp) = GaussYawDeflection(0.022, 0.022, 2.32, 0.154, interp)
GaussYawDeflection(a,b,c,d) = GaussYawDeflection(a, b, c, d, true)

"""
    GaussYawDeflectionVariableSpread(alpha_star, beta_star, k1, k2, interpolation)

Container for parameters related to the Gaussian deflection model with yaw presented by Bastankhah and Porte-Agel 2016

# Arguments
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
- `k1::Float`: first parameter tuning wake spread as based on turbulence intensity
- `k2::Float`: second parameter tuning wake spread as based on turbulence intensity
- `interpolation::Bool`: boolean stating if the the near wake should be interpolated. Default value is true.
"""
struct GaussYawVariableSpreadDeflection{TF1, TF2, TF3, TF4, BO} <: AbstractWakeDeflectionModel
    alpha_star::TF1
    beta_star::TF2
    k1::TF3
    k2::TF4
    interpolate_sigma::BO
end
GaussYawVariableSpreadDeflection() = GaussYawVariableSpreadDeflection(2.32, 0.154, 0.3837, 0.003678, true)
GaussYawVariableSpreadDeflection(interp) = GaussYawVariableSpreadDeflection(2.32, 0.154, 0.3837, 0.003678, interp)
GaussYawVariableSpreadDeflection(x, y) = GaussYawVariableSpreadDeflection(x, y, 0.3837, 0.003678, true)
GaussYawVariableSpreadDeflection(x, y, interp) = GaussYawVariableSpreadDeflection(x, y, 0.3837, 0.003678, interp)

"""
    JiminezYawDeflection(horizontal_spread_rate)

Container for parameters related to the Jiminez deflection model

# Arguments
- `horizontal_spread_rate::Float`: parameter controlling the wake spreading rate and deficit decay. Default value is 0.1
"""
struct JiminezYawDeflection{TF} <: AbstractWakeDeflectionModel
    horizontal_spread_rate::TF
end
JiminezYawDeflection() = JiminezYawDeflection(0.1)

"""
    MultizoneDeflection(horizontal_spread_rate)

Container for parameters related to the Jiminez deflection model

# Arguments
- `horizontal_spread_rate::Float`: parameter controlling the wake spreading rate and deficit decay. Default value is 0.1
- `ad::Float`:Helps define the horizontal deflection of the wake at 0 deg yaw
- `bd::Float`:Helps define the horizontal deflection of the wake due to downwind distance at 0 deg yaw
"""
struct MultizoneDeflection{TF1, TF2, TF3} <: AbstractWakeDeflectionModel
    horizontal_spread_rate::TF1
    ad::TF2
    bd::TF3
end
MultizoneDeflection() = MultizoneDeflection(0.15, -4.5, -0.01)

"""
    wake_deflection_model(locx, locy, locz, turbine_id, turbine_definition::TurbineDefinition, model::NoYawDeflection, windfarmstate::SingleWindFarmState)

    Bypasses yaw deflection calculations.

"""
function wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::NoYawDeflection)

    return 0.0

end

"""
    wake_deflection_model(locx, locy, locz, turbine_id, turbine_definition::TurbineDefinition, model::JiminezYawDeflection)

    Calculates the horizontal deflection of the wind turbine wake

    Based on:
    [1] Jiminez 2010 "Wake defl ection of a wind turbine in yaw"
    [2] Gebraad 2014 "Wind plant optimization by yaw control using a parametric wake model"
    this version ignores the corrections made to the yaw model for rotor rotation as described in [2] and
    [3] Thomas 2017 "Improving the FLORIS wind plant model for compatibility with gradient-based optimization"
"""
function wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::JiminezYawDeflection)

    dx = locx-turbine_x[turbine_id]
    yaw = -turbine_yaw[turbine_id] # Jiminez used opposite rotation convention, hence (-) sign
    ct = turbine_ct[turbine_id]
    diam = rotor_diameter[turbine_id]

    kd = model.horizontal_spread_rate

    initial_yaw_angle = 0.5*((cos(yaw))^2)*sin(yaw)*ct  # [1] eq. 20, [2] eq. 8

    # [2] eq. 10
    a = 2.0*kd*dx/diam + 1.0
    b = initial_yaw_angle*(15.0*a^4+initial_yaw_angle^2)
    c = (30.0*kd/diam)*a^5
    d = initial_yaw_angle*diam*(15.0 + initial_yaw_angle^2)
    e = 30.0*kd

    y_deflection = b/c - d/e

    return y_deflection

end

"""
    wake_deflection_model(locx, locy, locz, turbine_id, turbine_definition::TurbineDefinition, model::MultizoneDeflection, windfarmstate::SingleWindFarmState)

    Calculates the horizontal deflection of the wind turbine wake accounting for both yaw and rotational deflection

    Based on:
    [1] Jiminez 2010 "Wake defl ection of a wind turbine in yaw"
    [2] Gebraad 2014 "Wind plant optimization by yaw control using a parametric wake model"
    this version ignores the corrections made to the yaw model for rotor rotation as described in [2] and
    [3] Thomas 2017 "Improving the FLORIS wind plant model for compatibility with gradient-based optimization"
"""

function wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::MultizoneDeflection)

    dx = locx-turbine_x[turbine_id]
    yaw = -turbine_yaw[turbine_id] # Jiminez used opposite rotation convention, hence (-) sign
    ct = turbine_ct[turbine_id]
    diam = rotor_diameter[turbine_id]

    kd = model.horizontal_spread_rate
    ad = model.ad
    bd = model.bd

    initial_yaw_angle = 0.5*((cos(yaw))^2)*sin(yaw)*ct  # [1] eq. 20, [2] eq. 8

    # [2] eq. 10
    a = 2.0*kd*dx/diam + 1.0
    b = initial_yaw_angle*(15.0*a^4+initial_yaw_angle^2)
    c = (30.0*kd/diam)*a^5
    d = initial_yaw_angle*diam*(15.0 + initial_yaw_angle^2)
    e = 30.0*kd

    # [2] eq. 10, 11, and 12 define deflection
    yaw_deflection = b/c - d/e
    rotation_deflection = ad + bd*(dx)
    y_deflection = -1*(yaw_deflection + rotation_deflection)

    return y_deflection

end

function _bpa_theta_0(yaw, ct)

    theta0 = (0.3*yaw/cos(yaw))*(1.0-sqrt(1.0-ct*cos(yaw)))

    return theta0
end

function _bpa_deflection(diam, ct, yaw, ky, kz, sigmay, sigmaz, theta0, x0)
    diam_inv = 1.0/diam
    a = theta0*x0*diam_inv
    b1 = cos(yaw)
    b = (theta0/14.7)*sqrt(b1/(ky*kz*ct))*(2.9-1.3*sqrt(1.0-ct)-ct)

    c1 = 1.6*sqrt(8.0*sigmay*sigmaz * diam_inv * diam_inv / b1)
    c2 = sqrt(ct)
    
    # c = (1.6+c2)*(c1-ct)
    # d = (1.6-c2)*(c1+ct)
    g1 = 1.6*c1-c2*ct
    g2 = -1.6*ct+c1*c2
    c = g1+g2
    d = g1-g2

    y_deflection = diam*(a+b*log(c/d))
    return y_deflection
end

"""
    wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::GaussYawDeflection)

    Calculates the horizontal deflection of the wind turbine wake

    Based on:
    [1] Bastankhah and Porte-Agel 2016 "Experimental and theoretical study of
    wind turbine wakes in yawed conditions"
"""
function wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::GaussYawDeflection)

    dx = locx-turbine_x[turbine_id]
    yaw = turbine_yaw[turbine_id]
    ct = turbine_ct[turbine_id]
    diam = rotor_diameter[turbine_id]
    ti = turbine_local_ti[turbine_id]

    as = model.alpha_star
    bs = model.beta_star
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate

    # [1] eqn 6.12
    theta0 = _bpa_theta_0(yaw, ct)

    # [1] eqn 7.4
    x0 = _gauss_yaw_potential_core(diam, yaw, ct, as, ti, bs)

    # calculate the discontinuity point of the gauss yaw model
    xd = _gauss_yaw_discontinuity(diam, x0, ky, kz, yaw, ct)

    # calculate horizontal wake spread (paper eq: 7.2)
    sigmay = _gauss_yaw_spread_interpolated(diam, ky, dx, x0, yaw, xd)

    # calculate vertical wake spread (paper eq: 7.2)
    sigmaz = _gauss_yaw_spread_interpolated(diam, kz, dx, x0, 0.0, xd)


    y_deflection = _bpa_deflection(diam, ct, yaw, ky, kz, sigmay, sigmaz, theta0, x0)

    return y_deflection
end

"""
    wake_deflection_model(oc, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::GaussYawVariableSpreadDeflection)

    Calculates the horizontal deflection of the wind turbine wake. Varies based on local turbulence intensity.

    Based on:
    [1] Bastankhah and Porte-Agel 2016 "Experimental and theoretical study of
    wind turbine wakes in yawed conditions"
    [2] Niayifar and Porte-Agel 2016 "Analytical Modeling of Wind Farms:
    A New Approach for Power Prediction"
"""
function wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::GaussYawVariableSpreadDeflection)

    dx = locx-turbine_x[turbine_id]
    yaw = turbine_yaw[turbine_id]
    ct = turbine_ct[turbine_id]
    dt = rotor_diameter[turbine_id]
    ti = turbine_local_ti[turbine_id]

    as = model.alpha_star
    bs = model.beta_star

    # [2] calculate wake spread based on local turbulence intensity
    ky = kz = _k_star_func(ti, model.k1, model.k2)

    # [1] eqn 6.12 initial wake angle
    theta0 = _bpa_theta_0(yaw, ct)

    # [1] eqn 7.4
    x0 = _gauss_yaw_potential_core(dt, yaw, ct, as, ti, bs)

    # calculate the discontinuity point of the gauss yaw model
    xd = _gauss_yaw_discontinuity(dt, x0, ky, kz, yaw, ct)

    # calculate horizontal wake spread (paper eq: 7.2)
    sigma_y = _gauss_yaw_spread_interpolated(dt, ky, dx, x0, yaw, xd)

    # calculate vertical wake spread (paper eq: 7.2)
    sigma_z = _gauss_yaw_spread_interpolated(dt, kz, dx, x0, 0.0, xd)

    # finally, calculate deflection
    y_deflection = _bpa_deflection(dt, ct, yaw, ky, kz, sigma_y, sigma_z, theta0, x0)

    return y_deflection
end
