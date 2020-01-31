abstract type AbstractDeflectionModel end

struct GaussYawDeflection <: AbstractDeflectionModel
    turbulence_intensity
    horizontal_spread_rate
    vertical_spread_rate
    alpha_star
    beta_star
end

function deflection_model(loc, model::GaussYawDeflection, turbine::Turbine)
    # from on Bastankhah and Porte-Agel 2016 [1]

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
    dy = dt*(a+b*log(c/d))

    return dy
end

function deflection_model(loc, turbine::Turbine)

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    del = sqrt(dy^2+dz^2)

end
