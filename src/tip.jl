"""
Convert angle from cartesian coordinate system to meterological coordinate system
"""

function tip_cart2met(wind_direction_cart)
    # convert from cart. polar system (CCW, 0 rad.=E) to meteorological polar system (CW, 0 rad.=N) 
    wind_direction_met = 3.0*pi/2.0 - wind_direction_cart

    # adjust to be between 0 and 2pi
    if wind_direction_met < 0.0
        wind_direction_met += 2.0*pi
    end
    
    return wind_direction_met
end

function tip_wake_combination_model(deficit, wind_speed, old_deficit_sum)
    # Katic et al. 1986
    new_deficit_sum = sqrt(old_deficit_sum^2 + (wind_speed*deficit)^2)
    if new_deficit_sum > wind_speed
        new_deficit_sum = wind_speed
    end

    return new_deficit_sum

end

# get angle for two given turbines
function tip_get_angle(dx, dy)
    angle = atan(dy,dx)
end

# get probability for a given angle
function tip_get_probability_spline(wind_rose, step)
    pre_spline = Akima(wind_rose[:,1], wind_rose[:,3])
    p_spline(x) = pre_spline(x)*length(wind_rose[:,1])/(2*pi/step)
    return p_spline
end

# get speed for a given angle
function tip_get_speed_spline(wind_rose)
    speed_spline = Akima(wind_rose[:,1], wind_rose[:,2])
    return speed_spline
end

# get effective velocity for a given turbine
function tip_get_effvelocity(turbid, x, y, uave, probability_spline, speed_spline, ct_spline, pstep; r0=0.5, alpha=0.1)
    nturbines = length(x)
    deficit_sum = 0.0
    for j = 1:nturbines
        dx = x[j] - x[turbid]
        dy = y[j] - y[turbid]
        dr = sqrt(dx^2+dy^2)
        angle = tip_cart2met(tip_get_angle(dx, dy)) # radians
        
        wake_range = angle-atan(alpha):pstep:angle+atan(alpha)
        p = sum(probability_spline(wake_range))
        wind_speed = sum(speed_spline(wake_range))/length(wake_range)
        # counter = 0
        # for ang in (angle-beta*0.0):1:(angle+beta*0.0)
        #     p += probability_spline(ang)
        #     wind_speed += speed_spline(ang)
        #     counter += 1
        # end
        # wind_speed /= counter
        a = tip_get_ai(wind_speed, ct_spline)
        deficit = p*2.0*a*(r0/(r0+alpha*dr))^2
        
        deficit_sum = tip_wake_combination_model(deficit, wind_speed, deficit_sum)
    end
    
    ueff = uave - deficit_sum
    return ueff
end

# get effective cp for a turbine in a given direction
function tip_get_cp_spline(cpdata)
    cp_spline = Akima(cpdata[:,1], cpdata[:,2])
    return cp_spline
end
    
# create spline for thrust coefficient
function tip_get_ct_spline(ctdata)
    ct_spline = Akima(ctdata[:,1], ctdata[:,2])
    return ct_spline
end

# calculate ai given wind speed and a ct curve
function tip_get_ai(speed, ct_spline)
    ct = ct_spline(speed)
    if (ct > 0.96) # Glauert condition
        ai = 0.143 + sqrt(0.0203-0.6427*(0.889 - ct))
    else
        ai = 0.5(1 - sqrt(1-ct))
    end
    return ai
end

# get average wind speed
function tip_get_ave_wind_speed(wind_rose)
    wind_speeds = wind_rose[:,2]
    probabilities = wind_rose[:,3]
    ave_speed = sum(wind_speeds.*probabilities)
    return ave_speed
end

# get average power for a given turbine
function tip_get_effpower(turbid, x, y, uave, probability_spline, speed_spline, cp_spline, ct_spline, pstep; prated=5E6, u_cut_in=3.0, u_cut_out=25.0, density=1.225, r0=0.5, alpha=0.1)
    area = pi*r0^2
    ueff = tip_get_effvelocity(turbid, x, y, uave, probability_spline, speed_spline, ct_spline, pstep, r0=r0, alpha=alpha)
    peff = 0.0
    if ueff < u_cut_in
        return peff
    elseif ueff > u_cut_out
        return peff
    else 
        cp = cp_spline(ueff)
        peff = 0.5*density*area*cp*ueff^3
    end
        
    if peff > prated
        peff = prated
    end
    return peff
end

# get annual energy production
function tip_get_aep(x, y, wind_rose, cpdata, ctdata; pstep=0.01, prated=5E6, u_cut_in=3.0, u_cut_out=25.0, density=1.225, r0=0.5, alpha=0.1)
    nturbines = length(x)
    uave = tip_get_ave_wind_speed(wind_rose)
    power = zeros(nturbines)
    probability_spline = tip_get_probability_spline(wind_rose,pstep)
    speed_spline = tip_get_speed_spline(wind_rose)
    cp_spline = tip_get_cp_spline(cpdata)
    ct_spline = tip_get_ct_spline(ctdata)
    for i = 1:nturbines
        power[i] = tip_get_effpower(i, x, y, uave, probability_spline, speed_spline, cp_spline, ct_spline, pstep, prated=prated, u_cut_in=u_cut_in, u_cut_out=u_cut_out, density=density, r0=r0, alpha=alpha)
    end
    aep = (365.25)*(24.0)*sum(power)
    return aep
end