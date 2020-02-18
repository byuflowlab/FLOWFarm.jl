abstract type AbstractWakeDeflectionModel end

#TODO add Jiminez deflection

struct GaussYawDeflection{TF} <: AbstractWakeDeflectionModel
    turbulence_intensity::TF
    horizontal_spread_rate::TF
    vertical_spread_rate::TF
    alpha_star::TF
    beta_star::TF
end

struct JiminezYawDeflection{TF} <: AbstractWakeDeflectionModel
    horizontal_spread_rate::TF
end

function wake_deflection_model(loc, model::JiminezYawDeflection, turbine::Turbine)
    # based on:
    # [1] Jiminez 2010 "Wake defl ection of a wind turbine in yaw"
    # [2] Gebraad 2014 "Wind plant optimization by yaw control using a parametric wake model"
    # this version ignores the corrections made to the yaw model for rotor rotation as described in [2] and 
    # [3] Thomas 2017 "Improving the FLORIS wind plant model for compatibility with gradient-based optimization"

    ct = turbine.ct     # thurst coefficient
    yaw = -turbine.yaw   # yaw missalignment (rad)
    kd = model.horizontal_spread_rate
    dx = loc[1]-turbine.coord.x
    diam = turbine.rotor_diameter

    initial_yaw_angle = 0.5*((cos(yaw))^2)*sin(yaw)*ct  # [1] eq. 20, [2] eq. 8

    # [2] eq. 10
    a = 2.0*kd*dx/diam + 1.0
    b = initial_yaw_angle*(15.0*a^4+initial_yaw_angle^2)
    c = (30.0*kd/diam)*a^5
    d = initial_yaw_angle*diam*(15.0 + initial_yaw_angle^2)
    e = 30.0*kd

    y_defflection = b/c - d/e

    return y_defflection

end
function wake_deflection_model(loc, model::GaussYawDeflection, turbine::Turbine)
    # [1] Bastankhah and Porte-Agel 2016

    dx = loc[1]-turbine.coord.x
    yaw = turbine.yaw
    ct = turbine.ct
    dt = turbine.rotor_diameter

    as = model.alpha_star
    bs = model.beta_star
    ti = model.turbulence_intensity
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate

    # [1] eqn 6.12
    theta0 = (0.3*yaw/cos(yaw))*(1.0-sqrt(1.0-ct*cos(yaw)))

    # [1] eqn 7.4
    x0 = _gauss_yaw_potential_core(dt, yaw, ct, as, ti, bs)
    sigmay = _gauss_yaw_spread(dt, ky, dx, x0, yaw)
    sigmaz = _gauss_yaw_spread(dt, kz, dx, x0, 0.0)
    a = theta0*x0/dt
    b = (theta0/14.7)*sqrt(cos(yaw)/(ky*kz*ct))*(2.9-1.3*sqrt(1.0-ct)-ct)
    c = (1.6+sqrt(ct))*(1.6*sqrt(8.0*sigmay*sigmaz/(cos(yaw)*dt^2))-ct)
    d = (1.6-sqrt(ct))*(1.6*sqrt(8.0*sigmay*sigmaz/(cos(yaw)*dt^2))+ct)
    y_defflection = dt*(a+b*log(c/d))

    return y_defflection
end

function deflection_model(loc, turbine::Turbine)

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    del = sqrt(dy^2+dz^2)

end
