export PowerModelConstantCp, PowerModelCpPoints, PowerModelPowerPoints, PowerModelPowerCurveCubic
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
struct PowerModelCpPoints{ATF1,ATF2,TF} <: AbstractPowerModel
    vel_points::ATF1
    cp_points::ATF2
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
struct PowerModelPowerPoints{ATF1,ATF2, TF} <: AbstractPowerModel
    vel_points::ATF1
    power_points::ATF2
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

"""
    turbine_powers_one_direction((generator_efficiency, cut_in_speed, cut_out_speed,
        rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, air_density, power_model::AbstractPowerModel)

Calculate the power for all wind turbines for a given state

# Arguments
- `generator_efficiency::Array{Float,nTurbines}`
- `cut_in_speed::Array{Float,nTurbines}`
- `cut_out_speed::Array{Float,nTurbines}`
- `rated_speed::Array{Float,nTurbines}`
- `rated_power::Array{Float,nTurbines}`
- `rotor_diameter::Array{Float,nTurbines}`
- `turbine_inflow_velcities::Array{Float,nTurbines}`: for current state only
- `air_density::Float`
- `power_models::Array{nturbines})` elements of array should be be of sub-types or AbstractPowerModel
"""
function turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed,
    rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, turbine_yaw, air_density,
    power_models; jac=nothing)

    # get number of turbines and rotor sample point
    nturbines = length(turbine_inflow_velcities)

    arr_type = promote_type(typeof(generator_efficiency[1]),typeof(cut_in_speed[1]),typeof(cut_out_speed[1]),typeof(rated_speed[1]),
                            typeof(rated_power[1]),typeof(rotor_diameter[1]),typeof(turbine_inflow_velcities[1]),typeof(turbine_yaw[1]))
    wt_power = zeros(arr_type, nturbines)

    if jac === nothing
        for d=1:nturbines
            wt_power[d] = calculate_turbine_power(generator_efficiency[d], cut_in_speed[d], cut_out_speed[d], rated_speed[d], rated_power[d], rotor_diameter[d], turbine_inflow_velcities[d], turbine_yaw[d], power_models[d], air_density)
        end
    else
        for d=1:nturbines
            wt_power[d], jac[d] = calculate_turbine_power(generator_efficiency[d], cut_in_speed[d], cut_out_speed[d], rated_speed[d], rated_power[d], rotor_diameter[d], turbine_inflow_velcities[d], turbine_yaw[d], power_models[d], air_density; return_derivatives=true)
        end
    end
    return wt_power
end

"""
    calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0])

Calculate power for each turbine for all states respectively

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine east-west locations in the global
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine north-south locations in the global
    reference frame
- `turbine_z::Array{Float,nTurbines}`: turbine base height in the global reference frame
- `rotor_diameter::Array{Float,nTurbines}`
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in
    radians
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes
    with state etc
- `generator_efficiency::Array{Float,nTurbines}`
- `cut_in_speed::Array{Float,nTurbines}`
- `cut_out_speed::Array{Float,nTurbines}`
- `rated_speed::Array{Float,nTurbines}`
- `rated_power::Array{Float,nTurbines}`
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds,
    frequencies, etc)
- `power_models::Array{nTurbines}`: elemenst of array should be sub-types of AbstractPowerModel
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across
    the rotor swept area when calculating the effective wind speed for the wind turbine.
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
"""
function calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0])

    nturbines = length(turbine_x)

    wind_probabilities = wind_resource.wind_probabilities

    nstates = length(wind_probabilities)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),typeof(hub_height[1]),typeof(turbine_yaw[1]),
            typeof(generator_efficiency[1]),typeof(cut_in_speed[1]),typeof(cut_out_speed[1]),typeof(rated_speed[1]),typeof(rated_power[1]))

    turbine_powers_by_direction = zeros(arr_type,(nstates, nturbines))

    for i = 1:nstates
        rot_x, rot_y = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[i])

        sorted_turbine_index = sortperm(rot_x)

        turbine_velocities = turbine_velocities_one_direction(rot_x, rot_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                            model_set, wind_farm_state_id=i, velocity_only=true)

        wt_power = turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                            rated_power, rotor_diameter, turbine_velocities, turbine_yaw, wind_resource.air_density, power_models)
        turbine_powers_by_direction[i,:] = wt_power
    end

    return turbine_powers_by_direction
end

"""
    calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models::Array{AbstractPowerModel}, model_set::AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0])

Calculate AEP for each requested state respectively

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine east-west locations in the global
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine north-south locations in the global
    reference frame
- `turbine_z::Array{Float,nTurbines}`: turbine base height in the global reference frame
- `rotor_diameter::Array{Float,nTurbines}`
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in
    radians
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes
    with state etc
- `generator_efficiency::Array{Float,nTurbines}`
- `cut_in_speed::Array{Float,nTurbines}`
- `cut_out_speed::Array{Float,nTurbines}`
- `rated_speed::Array{Float,nTurbines}`
- `rated_power::Array{Float,nTurbines}`
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds,
    frequencies, etc)
- `power_model::Array{nTurbines}`: elemenst of array should be sub-types of AbstractPowerModel
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across
    the rotor swept area when calculating the effective wind speed for the wind turbine.
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `hours_per_year::Float`: hours per year (averaged for leap year by default)
"""
function calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
            hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
            cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::AbstractModelSet;
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0, weighted=true, wind_farm_state_id=1)

    # get number of states
    nstates = length(wind_resource.wind_probabilities)

    # set array type
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),typeof(hub_height[1]),typeof(turbine_yaw[1]),
                typeof(generator_efficiency[1]),typeof(cut_in_speed[1]),typeof(cut_out_speed[1]),typeof(rated_speed[1]),typeof(rated_power[1]))

    # initialize state energy
    state_energy = zeros(arr_type,nstates)

    # pre-allocate arrays for later calculations
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_yaw[1]))
    n_turbines = length(turbine_x)
    prealloc_turbine_velocities = zeros(arr_type, n_turbines)
    prealloc_turbine_ct = zeros(arr_type, n_turbines)
    prealloc_turbine_ai = zeros(arr_type, n_turbines)
    prealloc_turbine_local_ti = zeros(arr_type, n_turbines)

    # loop over all states
    for i = 1:nstates

        state_energy[i] = calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
            turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
            rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
            model_set; wind_farm_state_id=i, hours_per_year=hours_per_year, weighted=weighted, prealloc_turbine_velocities=prealloc_turbine_velocities,
            prealloc_turbine_ct=prealloc_turbine_ct, prealloc_turbine_ai=prealloc_turbine_ai, prealloc_turbine_local_ti=prealloc_turbine_local_ti)

    end

    return state_energy
end

function calculate_state_aep(turbine_x::T0, turbine_y::T1, turbine_z::T2, rotor_diameter::T3, hub_height::T4,
    turbine_yaw::T5, ct_model::Vector{<:AbstractThrustCoefficientModel}, generator_efficiency::Vector{T6}, cut_in_speed::Vector{T6}, cut_out_speed::Vector{T6}, rated_speed::Vector{T6},
    rated_power::Vector{T6}, power_models::Vector{<:AbstractPowerModel}, rotor_sample_points_y::Vector{T6}, rotor_sample_points_z::Vector{T6}, wind_resource,
    model_set; wind_farm_state_id=1, hours_per_year=365.25*24.0, weighted=true, wind_speed_ids=nothing, prealloc_turbine_velocities=nothing,
    prealloc_turbine_ct=nothing, prealloc_turbine_ai=nothing, prealloc_turbine_local_ti=nothing, prealloc_wake_deficits=nothing, prealloc_contribution_matrix=nothing,
    prealloc_deflections=nothing, prealloc_sigma_squared=nothing) where {T0, T1, T2, T3, T4, T5, T6}

    # rotate turbine locations to match the direction of the current state
    rot_x, rot_y = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # get turbine indices in sorted order from upstream to downstream
    sorted_turbine_index = sortperm(rot_x)

    if wind_speed_ids === nothing
        # calculate wind turbine velocities for given state
        prealloc_turbine_velocities = turbine_velocities_one_direction(rot_x, rot_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                            model_set; turbine_velocities=prealloc_turbine_velocities, turbine_ct=prealloc_turbine_ct, turbine_ai=prealloc_turbine_ai, turbine_local_ti=prealloc_turbine_local_ti,
                            wind_farm_state_id=wind_farm_state_id, velocity_only=true)

        # calculate wind turbine powers for given state
        wt_power = turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                            rated_power, rotor_diameter, prealloc_turbine_velocities, turbine_yaw, wind_resource.air_density, power_models)

        # calculate wind farm power for given state
        state_power = sum(wt_power)

        # calculate aep for given state
        state_aep = state_power*hours_per_year*wind_resource.wind_probabilities[wind_farm_state_id]

    else # wind_farm_state_id carries the wind direction and wind speed information

        # calculate wind turbine velocities for given direction and moderate speed
        prealloc_turbine_velocities = turbine_velocities_one_direction(rot_x, rot_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                            model_set; turbine_velocities=prealloc_turbine_velocities, turbine_ct=prealloc_turbine_ct, turbine_ai=prealloc_turbine_ai, turbine_local_ti=prealloc_turbine_local_ti,
                            wind_farm_state_id=wind_farm_state_id, velocity_only=true, wake_deficits=prealloc_wake_deficits,
                            contribution_matrix=prealloc_contribution_matrix,deflections=prealloc_deflections,sigma_squared=prealloc_sigma_squared)

        # back out turbine deficits from the turbine velocities
        turbine_deficits = prealloc_turbine_velocities./wind_resource.wind_speeds[wind_farm_state_id]

        # initialize state aep (which is actually the directional aep in this case)
        state_aep = 0.0

        # loop over all speeds for the given direction
        for i = eachindex(wind_speed_ids)

            # calculate turbine velocities for this wind speed based on the deficits
            prealloc_turbine_velocities = turbine_deficits.*wind_resource.wind_speeds[wind_speed_ids[i]]

            # calculate the power of the turbines at each wind speed for the given direction
            wt_power = turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                        rated_power, rotor_diameter, prealloc_turbine_velocities, turbine_yaw, wind_resource.air_density, power_models)

            # sum the turbine powers to get the powers for the current speed/direction combination
            state_power = sum(wt_power)

            # add the weighted state power to the state (directional) aep
            state_aep += state_power*hours_per_year*wind_resource.wind_probabilities[wind_speed_ids[i]]

        end

        # calculate weighted average for state power
        state_power = state_aep/hours_per_year

    end

    # return desired quantity
    if weighted
        # return state aep or directional aep
        return state_aep
    else
        # return state power or average directional power
        return state_power
    end

end


"""
    calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models::Array{AbstractPowerModel}, model_set::AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0)

Calculate wind farm AEP

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine east-west locations in the global
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine north-south locations in the global
    reference frame
- `turbine_z::Array{Float,nTurbines}`: turbine base height in the global reference frame
- `rotor_diameter::Array{Float,nTurbines}`
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in
    radians
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes
    with state etc
- `generator_efficiency::Array{Float,nTurbines}`
- `cut_in_speed::Array{Float,nTurbines}`
- `cut_out_speed::Array{Float,nTurbines}`
- `rated_speed::Array{Float,nTurbines}`
- `rated_power::Array{Float,nTurbines}`
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds,
    frequencies, etc)
- `power_model::Array{)`: elements of array should be sub types of AbstractPowerModel
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across
    the rotor swept area when calculating the effective wind speed for the wind turbine.
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `hours_per_year::Float`: hours per year (averaged for leap year by default)
"""
function calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
            hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
            cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::AbstractModelSet;
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0, prealloc_turbine_velocities=nothing,
            prealloc_turbine_ct=nothing, prealloc_turbine_ai=nothing, prealloc_turbine_local_ti=nothing, prealloc_wake_deficits=nothing,
            prealloc_contribution_matrix=nothing, prealloc_deflections=nothing, prealloc_sigma_squared=nothing,
            force_single_thread=false)

    # find how many wind states are being calculated
    nstates = length(wind_resource.wind_directions)
    n_turbines = length(turbine_x)

    # check if we can run a single wind speed for each direction to save time
    if typeof(model_set.wake_combination_model) == SumOfSquaresFreestreamSuperposition
        # find unique directions
        unique_directions = unique(wind_resource.wind_directions)

        # find how many unique directions there are
        ndirections = length(unique_directions)
    end

    # state_energy = Vector{typeof(wind_farm.turbine_x[1])}(undef,nstates)
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),typeof(hub_height[1]),typeof(turbine_yaw[1]),
                typeof(generator_efficiency[1]),typeof(cut_in_speed[1]),typeof(cut_out_speed[1]),typeof(rated_speed[1]),typeof(rated_power[1]))

    if arr_type == ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}
        reverse_diff = true
    else
        reverse_diff = false
    end

    AEP = arr_type(0.0)

    # calculate AEP in parallel using multi-threading
    if Threads.nthreads() > 1 && !reverse_diff && !force_single_thread
        n_threads = Threads.nthreads()
        if prealloc_turbine_velocities === nothing
            prealloc_turbine_velocities = zeros(arr_type,n_turbines,n_threads)
        end
        if prealloc_turbine_ct === nothing
            prealloc_turbine_ct = zeros(arr_type,n_turbines,n_threads)
        end
        if prealloc_turbine_ai === nothing
            prealloc_turbine_ai = zeros(arr_type,n_turbines,n_threads)
        end
        if prealloc_turbine_local_ti === nothing
            prealloc_turbine_local_ti = zeros(arr_type,n_turbines,n_threads)
        end
        if prealloc_wake_deficits === nothing
            prealloc_wake_deficits = zeros(arr_type,n_turbines,n_turbines,n_threads)
        end
        if prealloc_contribution_matrix === nothing
            prealloc_contribution_matrix = zeros(arr_type,n_turbines,n_turbines,n_threads)
        end
        if prealloc_deflections === nothing
            prealloc_deflections = zeros(arr_type,n_turbines,n_turbines,n_threads)
        end
        if prealloc_sigma_squared === nothing
            prealloc_sigma_squared = zeros(arr_type,n_turbines,n_turbines,n_threads)
        end

        if typeof(model_set.wake_combination_model) == SumOfSquaresFreestreamSuperposition
            state_aep = zeros(arr_type,ndirections)

            n_per_thread, rem = divrem(ndirections,n_threads)
            rem > 0 && (n_per_thread += 1)
            assignments = 1:n_per_thread:ndirections

            Threads.@threads for i_assignment in eachindex(assignments)

                i_start = assignments[i_assignment]
                i_stop = min(i_start+n_per_thread-1, ndirections)

                for i = i_start:i_stop
                    # get indices to all speeds corresponding to this unique direction
                    wind_speed_ids = findall(wind_resource.wind_directions .== unique_directions[i])

                    # take a speed in the middle (so it is not zero)
                    middle_id = wind_speed_ids[max(cld(length(wind_speed_ids),2),1)]

                    # get direction aep including all wind speeds for that direction
                    state_aep[i] = calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
                        turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                        rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                        model_set; wind_farm_state_id=middle_id, hours_per_year=hours_per_year, wind_speed_ids=wind_speed_ids,
                        prealloc_turbine_velocities=view(prealloc_turbine_velocities,:,i_assignment), prealloc_turbine_ct=view(prealloc_turbine_ct,:,i_assignment),
                        prealloc_turbine_ai=view(prealloc_turbine_ai,:,i_assignment), prealloc_turbine_local_ti=view(prealloc_turbine_local_ti,:,i_assignment),
                        prealloc_wake_deficits=view(prealloc_wake_deficits,:,:,i_assignment), prealloc_contribution_matrix=view(prealloc_contribution_matrix,:,:,i_assignment),
                        prealloc_deflections=view(prealloc_deflections,:,:,i_assignment), prealloc_sigma_squared=view(prealloc_sigma_squared,:,:,i_assignment))
                end
            end
        else
            state_aep = zeros(arr_type,nstates)
            n_per_thread, rem = divrem(nstates,n_threads)
            rem > 0 && (n_per_thread += 1)
            assignments = 1:n_per_thread:nstates

            Threads.@threads for i_assignment in eachindex(assignments)

                i_start = assignments[i_assignment]
                i_stop = min(i_start+n_per_thread-1, nstates)

                for i = i_start:i_stop

                    state_aep[i] = calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
                        turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                        rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                        model_set; wind_farm_state_id=i, hours_per_year=hours_per_year,
                        prealloc_turbine_velocities=view(prealloc_turbine_velocities,:,i_assignment), prealloc_turbine_ct=view(prealloc_turbine_ct,:,i_assignment),
                        prealloc_turbine_ai=view(prealloc_turbine_ai,:,i_assignment), prealloc_turbine_local_ti=view(prealloc_turbine_local_ti,:,i_assignment),
                        prealloc_wake_deficits=view(prealloc_wake_deficits,:,:,i_assignment), prealloc_contribution_matrix=view(prealloc_contribution_matrix,:,:,i_assignment),
                        prealloc_deflections=view(prealloc_deflections,:,:,i_assignment), prealloc_sigma_squared=view(prealloc_sigma_squared,:,:,i_assignment))
                end
            end
        end

        AEP = sum(state_aep)

    # calculate AEP using distributed processing
    elseif Distributed.nworkers() > 1 && !reverse_diff && !force_single_thread
        # if possible, avoid recalculating wakes for more than one speed in each direction
        if typeof(model_set.wake_combination_model) == SumOfSquaresFreestreamSuperposition
            state_aep = zeros(arr_type,ndirections)
            AEP = @sync @distributed (+) for i = 1:ndirections

                # get indices to all speeds corresponding to this unique direction
                wind_speed_ids = findall(wind_resource.wind_directions .== unique_directions[i])

                # take a speed in the middle (so it is not zero)
                middle_id = wind_speed_ids[max(cld(length(wind_speed_ids),2),1)]

                # get direction aep including all wind speeds for that direction
                calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=middle_id, hours_per_year=hours_per_year, wind_speed_ids=wind_speed_ids)
            end
        else
            AEP = @sync @distributed (+) for i = 1:nstates

                calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=i, hours_per_year=hours_per_year)
            end
        end
    else
        state_aep = arr_type(0.0)
        if prealloc_turbine_velocities === nothing
            prealloc_turbine_velocities = zeros(arr_type, n_turbines)
        end
        if prealloc_turbine_ct === nothing
            prealloc_turbine_ct = zeros(arr_type, n_turbines)
        end
        if prealloc_turbine_ai === nothing
            prealloc_turbine_ai = zeros(arr_type, n_turbines)
        end
        if prealloc_turbine_local_ti === nothing
            prealloc_turbine_local_ti = zeros(arr_type, n_turbines)
        end
        if prealloc_wake_deficits === nothing
            prealloc_wake_deficits = zeros(arr_type,n_turbines,n_turbines)
        end
        if prealloc_contribution_matrix === nothing
            prealloc_contribution_matrix = zeros(arr_type,n_turbines,n_turbines)
        end
        if prealloc_deflections === nothing
            prealloc_deflections = zeros(arr_type,n_turbines,n_turbines)
        end
        if prealloc_sigma_squared === nothing
            prealloc_sigma_squared = zeros(arr_type,n_turbines,n_turbines)
        end
        if typeof(model_set.wake_combination_model) == SumOfSquaresFreestreamSuperposition
            for i = 1:ndirections
                # get indices to all speeds corresponding to this unique direction
                wind_speed_ids = findall(wind_resource.wind_directions .== unique_directions[i])

                # take a speed in the middle (so it is not zero)
                middle_id = wind_speed_ids[max(cld(length(wind_speed_ids),2),1)]

                # get direction aep including all wind speeds for that direction
                state_aep += calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=middle_id, hours_per_year=hours_per_year, wind_speed_ids=wind_speed_ids,
                    prealloc_turbine_velocities=view(prealloc_turbine_velocities,:,1), prealloc_turbine_ct=view(prealloc_turbine_ct,:,1),
                    prealloc_turbine_ai=view(prealloc_turbine_ai,:,1), prealloc_turbine_local_ti=view(prealloc_turbine_local_ti,:,1),
                    prealloc_wake_deficits=view(prealloc_wake_deficits,:,:,1), prealloc_contribution_matrix=view(prealloc_contribution_matrix,:,:,1),
                    prealloc_deflections=view(prealloc_deflections,:,:,1), prealloc_sigma_squared=view(prealloc_sigma_squared,:,:,1))
            end
        else
            for i = 1:nstates
                state_aep += calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height,
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=i, hours_per_year=hours_per_year,
                    prealloc_turbine_velocities=view(prealloc_turbine_velocities,:,1), prealloc_turbine_ct=view(prealloc_turbine_ct,:,1),
                    prealloc_turbine_ai=view(prealloc_turbine_ai,:,1), prealloc_turbine_local_ti=view(prealloc_turbine_local_ti,:,1),
                    prealloc_wake_deficits=view(prealloc_wake_deficits,:,:,1), prealloc_contribution_matrix=view(prealloc_contribution_matrix,:,:,1),
                    prealloc_deflections=view(prealloc_deflections,:,:,1), prealloc_sigma_squared=view(prealloc_sigma_squared,:,:,1))
            end
        end

        AEP = state_aep
    end

    return AEP
end

"""
    calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models::Array{AbstractPowerModel}, model_set::AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0)

Calculate ideal wind farm AEP (AEP with no wake loss)

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine east-west locations in the global
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine north-south locations in the global
    reference frame
- `turbine_z::Array{Float,nTurbines}`: turbine base height in the global reference frame
- `rotor_diameter::Array{Float,nTurbines}`
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in
    radians
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes
    with state etc
- `generator_efficiency::Array{Float,nTurbines}`
- `cut_in_speed::Array{Float,nTurbines}`
- `cut_out_speed::Array{Float,nTurbines}`
- `rated_speed::Array{Float,nTurbines}`
- `rated_power::Array{Float,nTurbines}`
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds,
    frequencies, etc)
- `power_model::Array{)`: elements of array should be sub types of AbstractPowerModel
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across
    the rotor swept area when calculating the effective wind speed for the wind turbine.
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `hours_per_year::Float`: hours per year (averaged for leap year by default)
"""
function calculate_ideal_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0)

    no_wake_model_set = WindFarmModelSet(NoWakeDeficit(),model_set.wake_deflection_model,
        model_set.wake_combination_model,model_set.local_ti_model,model_set.point_velocity_average_factor)

    AEP = calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, zeros(length(turbine_yaw)), ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, no_wake_model_set;
    rotor_sample_points_y, rotor_sample_points_z, hours_per_year)

    return AEP
end
