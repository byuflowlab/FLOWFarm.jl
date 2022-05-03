abstract type AbstractPowerModel end

"""
    PowerModelConstantCp(cp)

Models will assume a constant cp value as provided

# Arguments
- `cp::Float`: constant power coefficient value
- 'pp::TI': exponent for adjusting power for wind turbine yaw
"""
struct PowerModelConstantCp{TF,TI} <: AbstractPowerModel
    cp::TF
    pp::TI
end
PowerModelConstantCp(x) = PowerModelConstantCp(x,2)

"""
    PowerModelCpPoints(vel_points, cp_points)

Models will use adjust cp based on cp curve using linear interpolation of provided points

# Arguments
- `vel_points::Array{N,Float}`: wind speed values in m/s
- `cp_points::Array{N,Float}`: power coefficient values corresponding to the provided speeds
- 'pp::TF': exponent for adjusting power for wind turbine yaw
"""
struct PowerModelCpPoints{ATF,TF} <: AbstractPowerModel
    vel_points::ATF
    cp_points::ATF
    pp::TF
end
PowerModelCpPoints(x,y) = PowerModelCpPoints(x, y, 2)

"""
    PowerModelPowerPoints(vel_points, cp_points)

Models will use adjust wind turbine power based on power curve using linear interpolation of 
provided points

# Arguments
- `vel_points::Array{N,Float}`: wind speed values in m/s
- `power_points::Array{N,Float}`: power values corresponding to the provided speeds
- 'pp::TF': exponent for adjusting power for wind turbine yaw
"""
struct PowerModelPowerPoints{ATF, TF} <: AbstractPowerModel
    vel_points::ATF
    power_points::ATF
    pp::TF
end
PowerModelPowerPoints(x,y) = PowerModelPowerPoints(x, y, 2)

"""
    PowerModelPowerCurveCubic()

Power will be calculated based on turbine specifications assuming a cubic power curve. Note
that this method is inherently incorrect and should only be used for theoretical purposes 
or after careful validation.

# Arguments
- 'pp::TF': exponent for adjusting power for wind turbine yaw

"""
struct PowerModelPowerCurveCubic{TF} <: AbstractPowerModel
    pp::TF
end
PowerModelPowerCurveCubic() = PowerModelPowerCurveCubic(2)
"""
    calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)

Calculate the power for a wind turbine based on standard theory for region 2

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `cp::Float`: Power coefficient of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
"""
function calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity, wt_yaw; pp=2)
    power = generator_efficiency*0.5*air_density*rotor_area*cp*(cos(wt_yaw)^pp)*wt_velocity^3
    return power
end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model)

Calculate the power for a wind turbine based on standard theory for region 2 using a constant cp

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelConstantCp`: Struct containing the cp value to be used in region 2
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, wt_yaw, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelConstantCp)

    if wt_velocity < cut_in_speed
        power = 0.0
    elseif wt_velocity < rated_speed
        # extract cp_value
        cp = power_model.cp
        pp = power_model.pp
        power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity, wt_yaw, pp=pp)
        if power > rated_power
            power = rated_power
        end
    elseif wt_velocity < cut_out_speed
        power = rated_power
    elseif wt_velocity > cut_out_speed
        power = 0.0
    end

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

Calculate the power for a wind turbine based on a cp curve with linear interpolation

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelCpPoints`: Struct containing the velocity and cp values defining the cp curve
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, wt_yaw,
    cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelCpPoints)

    # obtain exponent for adjusting cp to yaw
    pp = power_model.pp

    # initialize power to zero
    power = 0.0

    # use specs if inflow wind speed is less than the wind speeds provided in the power curve
    if wt_velocity < power_model.vel_points[1]

        # calculated wind turbine power
        if wt_velocity < cut_in_speed
            power = 0.0
        elseif wt_velocity < rated_speed
            # use cp value corresponding to lowest provided velocity point
            cp = power_model.cp_points[1]
            # println(cp)
            # calculate power
            power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity, wt_yaw, pp=pp)
        end

    # use cp points where provided
    elseif wt_velocity < power_model.vel_points[end]

        # estimate cp_value using linear interpolation
        cp = linear(power_model.vel_points, power_model.cp_points, wt_velocity)
        # println(cp)
        # calculate power
        power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity, wt_yaw, pp=pp)

    # use specs if above vel_points max
    else

        if wt_velocity <= cut_out_speed
            # use cp value corresponding to highest provided velocity point
            cp = power_model.cp_points[end]
            # println(cp)
            # calculate power
            power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity, wt_yaw, pp=pp)
        end

    end

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

Calculate the power for a wind turbine based on a pre-determined power curve with linear
    interpolation

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelPowerPoints`: Struct containing the velocity and power values
    defining the power curve
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, wt_yaw,
    cut_in_speed, rated_speed, cut_out_speed, rated_power, 
    power_model::PowerModelPowerPoints; return_derivative=false)

    # get exponent for yaw adjustment
    pp = power_model.pp

    power = 0.0
    dpower_dvelocity = 0.0
    # use specs if inflow wind speed is less than the wind speeds provided in the power curve
    if wt_velocity < power_model.vel_points[1]

        # calculated wind turbine power
        if wt_velocity < rated_speed
            # use power value corresponding to lowest provided velocity point
            power = (cos(wt_yaw)^pp)*linear([cut_in_speed, power_model.vel_points[1]], [0.0, power_model.power_points[1]], wt_velocity)
        end

    # use power points where provided
    elseif wt_velocity < power_model.vel_points[end]

        # calculate power
        power = (cos(wt_yaw)^pp)*linear(power_model.vel_points, power_model.power_points, wt_velocity)
        if return_derivative
            dpower_dvelocity = (cos(wt_yaw)^pp)*gradient(power_model.vel_points, power_model.power_points, wt_velocity)
        end
    # use specs if above vel_points max
    else

        if wt_velocity <= cut_out_speed
            # use power corresponding to highest wind speed provided
            power = (cos(wt_yaw)^pp)*power_model.power_points[end]
        end
    end

    if !return_derivative
        return power
    else
        return power, dpower_dvelocity
    end
end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

Calculates wind turbine power using a cubic estimation based on turbine specifications
    as defined in https://github.com/byuflowlab/iea37-wflo-casestudies/blob/master/cs3-4/iea37-cs3-announcement.pdf

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelPowerCurveCubic`: Empty struct
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, wt_yaw,
    cut_in_speed, rated_speed, cut_out_speed, rated_power, 
    power_model::PowerModelPowerCurveCubic)

    pp = power_model.pp

    if wt_velocity < cut_in_speed
        power = 0.0
    elseif wt_velocity < rated_speed
        power = rated_power*((wt_velocity - cut_in_speed)/(rated_speed - cut_in_speed))^3
        # power = rated_power*((wt_velocity)/(rated_speed))^3
    elseif wt_velocity < cut_out_speed
        power = rated_power
    elseif wt_velocity > cut_out_speed
        power = 0.0
    end

    return power*(cos(wt_yaw)^pp)

end

"""
    calculate_turbine_power(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, 
    rated_power, rotor_diameter, wt_velocity, power_model::AbstractPowerModel, air_density)

Calculate the power for all wind turbines. Dispaches to desired power model.

# Arguments
- `generator_efficiency::Array{Float,nTurbines}`
- `cut_in_speed::Array{Float,nTurbines}` 
- `cut_out_speed::Array{Float,nTurbines}`
- `rated_speed::Array{Float,nTurbines}`
- `rated_power::Array{Float,nTurbines}`
- `rotor_diameter::Array{Float,nTurbines}`
- `wt_velocity::Array{Float,nTurbines}`: turbine effective wind speeds for current state only
- `power_model::AbstractPowerModel)
- `air_density::Float`
"""
function calculate_turbine_power(generator_efficiency, cut_in_speed, cut_out_speed, 
    rated_speed, rated_power, rotor_diameter, wt_velocity, wt_yaw, power_model::AbstractPowerModel, 
    air_density; return_derivatives=false)

    # calculated wind turbine rotor-swept area
    rotor_area = pi*(rotor_diameter^2)/4.0

    if !return_derivatives
        wt_power = calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, wt_yaw, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)
        return wt_power
    else
        wt_power, dp_du = calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, wt_yaw, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model, return_derivative=return_derivatives)
        return wt_power, dp_du
    end
end


