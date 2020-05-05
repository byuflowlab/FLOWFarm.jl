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

# struct PowerModelPowerCurveCubic{} <: AbstractPowerModel
# end

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
- `power_model::PowerModelConstantCp`: Struct containing the cp value to be used in region 2
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model::PowerModelConstantCp)

    # extract cp_value
    cp = power_model.cp

    power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model)

Calculate the power for a wind turbine based on standard theory for region 2 using a cp 
    curve with linear interpolation

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `power_model::PowerModelCpPoints`: Struct containing the velocity and cp values defining the cp curve
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model::PowerModelCpPoints)

    # estimate cp_value using linear interpolation
    cp = linear(power_model.vel_points, power_model.cp_points, wt_velocity)

    # calculate power
    power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model)

Calculate the power for a wind turbine based on a pre-determined power curve with linear 
    interpolation

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `power_model::PowerModelPowerPoints`: Struct containing the velocity and power values 
    defining the power curve
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model::PowerModelPowerPoints)

    # estimate power using linear interpolation
    power = linear(power_model.vel_points, power_model.power_points, wt_velocity)

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
function calculate_turbine_power(turbine_id, turbine_definition::TurbineDefinition, farmstate::SingleWindFarmState, wind_model::AbstractWindResourceModel)

    # extract turbine design information
    generator_efficiency = turbine_definition.generator_efficiency[1]
    cut_in_speed = turbine_definition.cut_in_speed
    cut_out_speed = turbine_definition.cut_out_speed
    rated_speed = turbine_definition.rated_speed
    rated_power = turbine_definition.rated_power
    rotor_diameter = turbine_definition.rotor_diameter[1]

    # calculated wind turbine rotor-swept area
    rotor_area = pi*(rotor_diameter^2)/4.0

    # extract wind resource information
    air_density = wind_model.air_density
    
    # extract wind turbine inflow velocity
    wt_velocity = farmstate.turbine_inflow_velcities[turbine_id]

    # extract wind turbine power model
    power_model = turbine_definition.power_model

    # calculated wind turbine power
    if wt_velocity < cut_in_speed[1]
        wt_power = 0.0
    elseif wt_velocity < rated_speed[1]
        wt_power = calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model)
    elseif wt_velocity < cut_out_speed[1]
        wt_power = rated_power[1]
    elseif wt_velocity > cut_out_speed[1]
        wt_power = 0.0
    end

    # adjust calculated power to not get higher than rated power
    if wt_power > rated_power[1]
        wt_power = 0.0
    end

    return wt_power

end

# function calculate_turbine_power(turbine::AbstractTurbine, farmstate::SingleWindFarmState, wind_model::AbstractWindResourceModel, power_model::PowerCurveCubic)
#
#     id = turbine.id
#
#     cut_in_speed = power_model.cut_in_speed
#     cut_out_speed = power_model.cut_out_speed
#     rated_speed = power_model.rated_speed
#     rated_power = power_model.rated_power
#
#     generator_efficiency = power_model.generator_efficiency[1]
#     rotor_diameter = turbine.rotor_diameter[1]
#     rotor_area = pi*(rotor_diameter^2)/4.0
#     wt_velocity = farmstate.turbine_inflow_velcities[turbine.id[1]]
#
#     wt_power = generator_efficiency*(0.5*air_density*rotor_area*cp*wt_velocity^3)
#
#     return wt_power
#
# end

function calculate_aep(windfarm::AbstractWindFarmModel)

end
