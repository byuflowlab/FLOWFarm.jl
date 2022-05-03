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

function calculate_state_aep(turbine_x::Vector{T0}, turbine_y::Vector{T1}, turbine_z::Vector{T2}, rotor_diameter::Vector{T3}, hub_height::Vector{T4}, 
    turbine_yaw::Vector{T5}, ct_model::Vector{<:AbstractThrustCoefficientModel}, generator_efficiency::Vector{T6}, cut_in_speed::Vector{T6}, cut_out_speed::Vector{T6}, rated_speed::Vector{T6},
    rated_power::Vector{T6}, power_models::Vector{<:AbstractPowerModel}, rotor_sample_points_y::Vector{T6}, rotor_sample_points_z::Vector{T6}, wind_resource,
    model_set; wind_farm_state_id=1, hours_per_year=365.25*24.0, weighted=true, wind_speed_ids=nothing, prealloc_turbine_velocities=nothing,
    prealloc_turbine_ct=nothing, prealloc_turbine_ai=nothing, prealloc_turbine_local_ti=nothing) where {T0, T1, T2, T3, T4, T5, T6}

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
                            wind_farm_state_id=wind_farm_state_id, velocity_only=true)

        # back out turbine deficits from the turbine velocities
        turbine_deficits = prealloc_turbine_velocities./wind_resource.wind_speeds[wind_farm_state_id]

        # initialize state aep (which is actualy the directional aep in this case)
        state_aep = 0.0

        # loop over all speeds for the given direction
        for i = 1:length(wind_speed_ids)

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
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0, distributed=false)

    # find how many wind states are being calculated
    nstates = length(wind_resource.wind_directions)

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

    # calculate AEP in parallel using multi-threading
    if Threads.nthreads() > 1
        if typeof(model_set.wake_combination_model) == SumOfSquaresFreestreamSuperposition
            state_aep = zeros(arr_type,ndirections)
            Threads.@threads for i = 1:ndirections

                # get indices to all speeds corresponding to this unique direction
                wind_speed_ids = findall(wind_resource.wind_directions .== unique_directions[i])

                # take a speed in the middle (so it is not zero)
                middle_id = wind_speed_ids[Int(length(wind_speed_ids)/2.0)]
                
                # get direction aep including all wind speeds for that direction
                state_aep[i] = calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, 
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=middle_id, hours_per_year=hours_per_year, wind_speed_ids=wind_speed_ids)
            end
        else
            state_aep = zeros(arr_type,nstates)
            Threads.@threads for i = 1:nstates

                state_aep[i] = calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, 
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=i, hours_per_year=hours_per_year)
            end
        end

        AEP = sum(state_aep)
        
    # calculate AEP in serial or in parallel using distributed processing
    else
        # if possible, avoid recalculating wakes for more than one speed in each direction
        if typeof(model_set.wake_combination_model) == SumOfSquaresFreestreamSuperposition
            state_aep = zeros(arr_type,ndirections)
            AEP = @sync @distributed (+) for i = 1:ndirections

                # get indices to all speeds corresponding to this unique direction
                wind_speed_ids = findall(wind_resource.wind_directions .== unique_directions[i])

                # take a speed in the middle (so it is not zero)
                middle_id = wind_speed_ids[max(Int(round(length(wind_speed_ids)/2.0)),1)]
                
                # get direction aep including all wind speeds for that direction
                calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, 
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=middle_id, hours_per_year=hours_per_year, wind_speed_ids=wind_speed_ids)
            end
        else
            AEP = @sync @distributed (+) for i = 1:nstates
            
                state_aep = calculate_state_aep(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, 
                    turbine_yaw, ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                    rated_power, power_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set; wind_farm_state_id=i, hours_per_year=hours_per_year)
            end
        end
    end

    return AEP
end