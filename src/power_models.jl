abstract type AbstractPowerModel end

struct PowerModelConstantCp{TF} <: AbstractPowerModel
    cp::TF
end

struct PowerModelCpPoints{ATF} <: AbstractPowerModel
    vel_points::ATF
    cp_points::ATF
end

struct PowerModelPowerPoints{ATF} <: AbstractPowerModel
    vel_points::ATF
    power_points::ATF
end

struct PowerModelPowerCurveCubic{} <: AbstractPowerModel
end

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
function calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)
    power = generator_efficiency*(0.5*air_density*rotor_area*cp*wt_velocity^3)
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
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelConstantCp)

    if wt_velocity < cut_in_speed
        power = 0.0
    elseif wt_velocity < rated_speed
        # extract cp_value
        cp = power_model.cp
        power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)
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
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelCpPoints)

    # use specs if inflow wind speed is less than the wind speeds provided in the power curve
    if wt_velocity < power_model.vel_points[1]

        # calculated wind turbine power
        if wt_velocity < cut_in_speed
            power = 0.0
        elseif wt_velocity < rated_speed
            # use cp value corresponding to lowest provided velocity point
            cp = power_model.cp_points[1]
            # calculate power
            power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)
        end

    # use cp points where provided
    elseif wt_velocity < power_model.vel_points[end]

        # estimate cp_value using linear interpolation
        cp = linear(power_model.vel_points, power_model.cp_points, wt_velocity)

        # calculate power
        power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)

    # use specs if above vel_points max
    else

        if wt_velocity <= cut_out_speed
            # use cp value corresponding to highest provided velocity point
            cp = power_model.cp_points[end]
            # calculate power
            power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)
        elseif wt_velocity > cut_out_speed
            power = 0.0
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
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelPowerPoints)

    # use specs if inflow wind speed is less than the wind speeds provided in the power curve
    if wt_velocity < power_model.vel_points[1]

        # calculated wind turbine power
        if wt_velocity < cut_in_speed
            power = 0.0
        elseif wt_velocity < rated_speed
            # use power value corresponding to lowest provided velocity point
            power = linear([cut_in_speed, power_model.vel_points[1]], [0.0, power_model.power_points[1]], wt_velocity)
        end

    # use power points where provided
    elseif wt_velocity < power_model.vel_points[end]

        # calculate power
        power = linear(power_model.vel_points, power_model.power_points, wt_velocity)

    # use specs if above vel_points max
    else

        if wt_velocity <= cut_out_speed
            # use power corresponding to highest wind speed provided
            power = power_model.power_points[end]
        elseif wt_velocity > cut_out_speed
            power = 0.0
        end

    end

    return power

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
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelPowerCurveCubic)


    if wt_velocity < cut_in_speed
        power = 0.0
    elseif wt_velocity < rated_speed
        # power = rated_power*((wt_velocity - cut_in_speed)/(rated_speed - cut_in_speed))^3
        power = rated_power*((wt_velocity)/(rated_speed))^3
    elseif wt_velocity < cut_out_speed
        power = rated_power
    elseif wt_velocity > cut_out_speed
        power = 0.0
    end

    return power

end

"""
    calculate_turbine_power(turbine_id, turbine_definition::TurbineDefinition,
        farmstate::SingleWindFarmState, wind_model::AbstractWindResourceModel)

Calculate the power for a wind turbine based on a pre-determined power curve with linear
    interpolation

# Arguments
- `turbine_id::Int`: Efficiency of the turbine generator
- `turbine_definition::TurbineDefinition`: Struct containing the relevant wind turbine deffinition
- `farmstate::SingleWindFarmState`: Struct contatining the current wind farm state, including correct inflow velocities
- `wind_model::AbstractWindResourceModel`: Struct defining the wind resource
"""
function calculate_turbine_power(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, rotor_diameter, wt_velocity, power_model::AbstractPowerModel, air_density)

    # calculated wind turbine rotor-swept area
    rotor_area = pi*(rotor_diameter^2)/4.0

    wt_power = calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

    return wt_power
end


function turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, air_density, power_model::AbstractPowerModel)

    # get number of turbines and rotor sample point
    nturbines = length(rotor_diameter)
    wt_power = zeros(nturbines)

    for d=1:nturbines
        wt_power[d] = calculate_turbine_power(generator_efficiency[d], cut_in_speed[d], cut_out_speed[d], rated_speed[d], rated_power[d], rotor_diameter[d], turbine_inflow_velcities[d], power_model, air_density)
    end
    return wt_power
end


"""
    calculate_aep(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem;
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0])

Calculate wind farm AEP


# Arguments
- `turbine_id::Int`: Efficiency of the turbine generator
- `turbine_definition::TurbineDefinition`: Struct containing the relevant wind turbine deffinition
- `farmstate::SingleWindFarmState`: Struct contatining the current wind farm state, including correct inflow velocities
- `wind_model::AbstractWindResourceModel`: Struct defining the wind resource
"""


function calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
            hub_height, turbine_yaw, turbine_ai, ct_model, generator_efficiency, cut_in_speed,
            cut_out_speed, rated_speed, rated_power, wind_resource, power_model::AbstractPowerModel, model_set::AbstractModelSet;
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], k1=0.2, k2=0.003)

    wind_probabilities = wind_resource.wind_probabilities

    nstates = length(wind_probabilities)
    hours_per_year = 365.25*24.0

    # state_energy = Vector{typeof(wind_farm.turbine_x[1])}(undef,nstates)
    state_energy = zeros(nstates)
    for i = 1:nstates

        rot_x, rot_y = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[i])

        sorted_turbine_index = sortperm(rot_x)

        turbine_velocities, turbine_ct, turbine_local_ti = turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                            turbine_ai, sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                            model_set; wind_farm_state_id=i, k1=k1, k2=k2)

        wt_power = turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                            rated_power, rotor_diameter, turbine_velocities, wind_resource.air_density, power_model)

        state_power = sum(wt_power)
        state_energy[i] = state_power*hours_per_year*wind_probabilities[i]
    end

    return sum(state_energy)
end
