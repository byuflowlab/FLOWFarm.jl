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
struct GaussYawDeflection{TF, BO} <: AbstractWakeDeflectionModel
    horizontal_spread_rate::TF
    vertical_spread_rate::TF
    alpha_star::TF
    beta_star::TF
    interpolate_sigma::BO
end
GaussYawDeflection() = GaussYawDeflection(0.022, 0.022, 2.32, 0.154, true)
GaussYawDeflection(interp) = GaussYawDeflection(0.022, 0.022, 2.32, 0.154, interp)
GaussYawDeflection(a,b,c,d) = GaussYawDeflection(a, b, c, d, true)
GaussYawDeflection(a,b,c,d,interp) = GaussYawDeflection(a, b, c, d, interp)

"""
    GaussTiltDeflection(turbulence_intensity, horizontal_spread_rate_tilt_1, horizontal_spread_rate_tilt_2, vertical_spread_rate_tilt_upper_1, vertical_spread_rate_tilt_upper_2, vertical_spread_rate_tilt_lower_1, vertical_spread_rate_tilt_lower_2, sigma_z0_upper_1, sigma_z0_upper_2, sigma_z0_upper_3, sigma_z0_lower_1, sigma_z0_lower_2, sigma_z0_lower_3, sigma_y0_1, sigma_y0_2, sigma_y0_3)

Container for parameters related to the Bastankhah Wake Model modified to include tilt

# Arguments
- `turbulence_intensity::Float`: Turbulence intensity usually set to around 0.09 (somewhat dependent on wind speed)
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
- `horizontal_spread_rate_tilt_1::Float`: ky is a function of tilt and this is the slope of the linear fit
- `horizontal_spread_rate_tilt_2::Float`: ky is a function if tilt and this is the intercept of the linear fit
- `vertical_spread_rate_tilt_upper_1::Float`: slope for the linear fit for kz of the upper portion of the wake
- `vertical_spread_rate_tilt_upper_2::Float`: intercept for the linear fit for kz of the upper portion of the wake
- `vertical_spread_rate_tilt_lower_1::Float`: slope for the linear fit for kz of the lower portion of the wake
- `vertical_spread_rate_tilt_lower_2::Float`: intercept for the linear fit for kz of the lower portion of the wake
- 'sigma_z0_upper_1::Float': used to define sigma_z0 as a function of tilt angle for the upper portion of the wake
- 'sigma_z0_upper_2::Float': used to define sigma_z0 as a function of tilt angle for the upper portion of the wake
- 'sigma_z0_upper_3::Float': used to define sigma_z0 as a function of tilt angle for the upper portion of the wake
- 'sigma_z0_lower_1::Float': used to define sigma_z0 as a function of tilt angle for the lower portion of the wake
- 'sigma_z0_lower_2::Float': used to define sigma_z0 as a function of tilt angle for the lower portion of the wake
- 'sigma_z0_lower_3::Float': used to define sigma_z0 as a function of tilt angle for the lower portion of the wake
- 'sigma_y0_1::Float': used to define sigma_y0 as a function of tilt angle
- 'sigma_y0_2::Float': used to define sigma_y0 as a function of tilt angle
- 'sigma_y0_3::Float': used to define sigma_y0 as a function of tilt angle
- `interpolation::Bool`: boolean stating if the the near wake should be interpolated. Default value is true.
"""
struct GaussTiltDeflection{TF,ATF,BO} <: AbstractWakeDeficitModel
    alpha_star::TF
    beta_star::TF
    ky1::TF
    ky2::TF
    kz1_up::TF
    kz2_up::TF
    kz1::TF
    kz2::TF
    sigy1::TF
    sigy2::TF
    sigy3::TF
    sigz1_up::TF
    sigz2_up::TF
    sigz3_up::TF
    sigz1::TF
    sigz2::TF
    sigz3::TF
    wec_factor::ATF
    interpolate_sigma::BO
end
GaussTiltDeflection() = GaussTiltDeflection(2.32, 0.154, 0.04666, 0.02229, -0.0352, 0.02071, 0.00655, 0.00029, 0.2608, -0.4913, 2.6534, 0.3536, 0.1766, -3.8565, 0.2473, -0.1921, 0.9547, [1.0], true)
GaussTiltDeflection(interp) = GaussTiltDeflection(2.32, 0.154, 0.04666, 0.02229, -0.0352, 0.02071, 0.00655, 0.00029, 0.2608, -0.4913, 2.6534, 0.3536, 0.1766, -3.8565, 0.2473, -0.1921, 0.9547, [1.0], interp)
GaussTiltDeflection(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) = GaussTiltDeflection(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,true)
GaussTiltDeflection(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,interp) = GaussTiltDeflection(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,interp)


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
struct GaussYawVariableSpreadDeflection{TF, BO} <: AbstractWakeDeflectionModel
    alpha_star::TF
    beta_star::TF
    k1::TF
    k2::TF
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
struct MultizoneDeflection{TF} <: AbstractWakeDeflectionModel
    horizontal_spread_rate::TF
    ad::TF
    bd::TF
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
    a = theta0*x0/diam
    b = (theta0/14.7)*sqrt(cos(yaw)/(ky*kz*ct))*(2.9-1.3*sqrt(1.0-ct)-ct)
    c = (1.6+sqrt(ct))*(1.6*sqrt(8.0*sigmay*sigmaz/(cos(yaw)*diam^2))-ct)
    d = (1.6-sqrt(ct))*(1.6*sqrt(8.0*sigmay*sigmaz/(cos(yaw)*diam^2))+ct)
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

"""
    wake_deflection_model(locx, locy, locz, turbine_x, turbine_tilt, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::GaussTiltDeflection)

    Calculates the vertical deflection of the wind turbine wake using a surrogate model
"""
function wake_deflection_model(locx, locy, locz, turbine_x, turbine_tilt, turbine_ct, turbine_id, rotor_diameter, turbine_local_ti, model::GaussTiltDeflection)

    dx = locx-turbine_x[turbine_id]
    tilt = turbine_tilt[turbine_id]
    ct = turbine_ct[turbine_id]
    diam = rotor_diameter[turbine_id]
    ti = turbine_local_ti[turbine_id]

    # extract model parameters
    # ks = model.k_star       # wake spread rate (k* in 2014 paper)
    as = model.alpha_star
    bs = model.beta_star
    ky1 = model.ky1
    ky2 = model.ky2
    kz1_up = model.kz1_up
    kz2_up = model.kz2_up
    kz1 = model.kz1
    kz2 = model.kz2
    sigy1 = model.sigy1
    sigy2 = model.sigy2
    sigy3 = model.sigy3
    sigz1_up = model.sigz1_up
    sigz2_up = model.sigz2_up
    sigz3_up = model.sigz3_up
    sigz1 = model.sigz1
    sigz2 = model.sigz2
    sigz3 = model.sigz3
    wec_factor = model.wec_factor[1]

    # [1] eqn 7.4 (currently the same for tilt and yaw)
    # There should also potentially be a surrogate model for x0
    x0 = _gauss_yaw_potential_core(diam, tilt, ct, as, ti, bs)

    #### Put Tilt Deflection surrogate model here 
    #### Will be a function of tilt and downstream distance
    #### Someday it should also be a function of turbulence intensity
    z_deflection = _gauss_tilt_surrogate(tilt, x0, diam, locx)

    # # Need to find ky and kz, but they rely on z_deflection

    # # calculate the discontinuity point (currently the same for tilt and yaw)
    # xd = _gauss_yaw_discontinuity(diam, x0, ky, kz, tilt, ct)
    
    # # calculate horizontal wake spread (paper eq: 7.2)
    # sigmay = _gauss_yaw_spread_interpolated(diam, ky, dx, x0, yaw, xd)

    # # calculate vertical wake spread (paper eq: 7.2)
    # sigmaz = _gauss_yaw_spread_interpolated(diam, kz, dx, x0, 0.0, xd)

    
    # z_deflection = _bpa_deflection(diam, ct, yaw, ky, kz, sigmay, sigmaz, theta0, x0)

    return z_deflection
end