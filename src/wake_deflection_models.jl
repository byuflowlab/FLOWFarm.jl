abstract type AbstractWakeDeflectionModel end

"""
    GaussYawDeflection(turbulence_intensity, horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star)

Container for parameters related to the Gaussian deflection model presented by Bastankhah and Porte-Agel 2016

# Arguments
- `turbulence_intensity::Float`: the ambient turbulence intensity. No default value is provided.
- `horizontal_spread_rate::Float`: parameter controlling the horizontal spread of the deficit model. Default value is 0.022.
- `vertical_spread_rate::Float`: parameter controlling the vertical spread of the deficit model. Default value is 0.022.
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
"""
struct GaussYawDeflection{TF} <: AbstractWakeDeflectionModel
    turbulence_intensity::TF
    horizontal_spread_rate::TF
    vertical_spread_rate::TF
    alpha_star::TF
    beta_star::TF
end
GaussYawDeflection(x) = GaussYawDeflection(x, 0.022, 0.022, 2.32, 0.154)

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
    wake_deflection_model(loc, turbine_id, turbine_definition::TurbineDefinition, model::JiminezYawDeflection, windfarmstate::SingleWindFarmState)
    
    Calculates the horizontal deflection of the wind turbine wake

    Based on:
    [1] Jiminez 2010 "Wake defl ection of a wind turbine in yaw"
    [2] Gebraad 2014 "Wind plant optimization by yaw control using a parametric wake model"
    this version ignores the corrections made to the yaw model for rotor rotation as described in [2] and 
    [3] Thomas 2017 "Improving the FLORIS wind plant model for compatibility with gradient-based optimization"
"""
function wake_deflection_model(loc, turbine_id, turbine_definition::TurbineDefinition, model::JiminezYawDeflection, windfarmstate::SingleWindFarmState)

    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    yaw = -windfarmstate.turbine_yaw[turbine_id] # Jiminez used opposite rotation convention, hence (-) sign
    ct = windfarmstate.turbine_ct[turbine_id]
    diam = turbine_definition.rotor_diameter[1]

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

function _bpa_theta_0(yaw, ct)

    theta0 = (0.3*yaw/cos(yaw))*(1.0-sqrt(1.0-ct*cos(yaw)))

    return theta0
end

function _bpa_deflection(d, ct, yaw, ky, kz, sigmay, sigmaz, theta0, x0)
    a = theta0*x0/d
    b = (theta0/14.7)*sqrt(cos(yaw)/(ky*kz*ct))*(2.9-1.3*sqrt(1.0-ct)-ct)
    c = (1.6+sqrt(ct))*(1.6*sqrt(8.0*sigmay*sigmaz/(cos(yaw)*d^2))-ct)
    d = (1.6-sqrt(ct))*(1.6*sqrt(8.0*sigmay*sigmaz/(cos(yaw)*d^2))+ct)
    y_deflection = d*(a+b*log(c/d))
    return y_deflection
end

"""
    wake_deflection_model(loc, turbine_id, turbine_definition::TurbineDefinition, model::GaussYawDeflection, windfarmstate::SingleWindFarmState)
    
    Calculates the horizontal deflection of the wind turbine wake

    Based on:
    [1] Bastankhah and Porte-Agel 2016 "Experimental and theoretical study of
    wind turbine wakes in yawed conditions"
"""
function wake_deflection_model(loc, turbine_id, turbine_definition::TurbineDefinition, model::GaussYawDeflection, windfarmstate::SingleWindFarmState)

    dx = loc[1]-windfarmstate.turbine_x[turbine_id]
    yaw = windfarmstate.turbine_yaw[turbine_id]
    ct = windfarmstate.turbine_ct[turbine_id]
    diam = turbine_definition.rotor_diameter[1]

    as = model.alpha_star
    bs = model.beta_star
    ti = model.turbulence_intensity
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate

    # [1] eqn 6.12
    theta0 = _bpa_theta_0(yaw, ct)

    # [1] eqn 7.4
    x0 = _gauss_yaw_potential_core(diam, yaw, ct, as, ti, bs)
    sigmay = _gauss_yaw_spread(diam, ky, dx, x0, yaw)
    sigmaz = _gauss_yaw_spread(diam, kz, dx, x0, 0.0)
    y_deflection = _bpa_deflection(d, ct, yaw, ky, kz, sigmay, sigmaz, theta0, x0)

    return y_deflection
end

# function deflection_model(loc, turbine::Turbine)

#     dx = loc[1]-turbine.coord.x
#     dy = loc[2]-(turbine.coord.y+deflection_y)
#     dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

#     del = sqrt(dy^2+dz^2)

# end
