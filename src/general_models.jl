abstract type AbstractModelSet end
# using CSV
# using DataFrames

"""
    WindFarmModelSet(wakedeficitmodel, wake_deflection_model, wake_combination_model, local_ti_model)

Container for objects defining models to use in wind farm calculations

# Arguments
- `wake_defiict_model::AbstractWakeDeficitModel`: contains a struct defining the desired wake deficit model
- `wake_deflection_model::AbstractWakeDeflectionModel`: contains a struct defining the desired wake deflection model
- `wake_combination_model::AbstractWakeCombinationModel`: contains a struct defining the desired wake combination model
- `local_ti_model::AbstractTurbulenceIntensityModel`: contains a struct defining the desired turbulence intensity model
"""
struct WindFarmModelSet{DTM,DNM,CM,TIM} <: AbstractModelSet

    wake_deficit_model::DTM
    wake_deflection_model::DNM
    wake_combination_model::CM
    local_ti_model::TIM

end

"""
    point_velocity(loc, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource, model_set::AbstractModelSet;
    wind_farm_state_id=1, downwind_turbine_id=0)

Calculates the wind speed at a given point for a given state

# Arguments
- `loc::Array{TF,3}`: Location of interest
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_tilt::Array{TF,nTurbines}`: turbine tilt for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: turbine thrust coefficients for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: array containing indices of wind turbines 
    from most upwind to most downwind turbine in the given state
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: contains wind resource discreption (directions,
    speeds, frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `downwind_turbine_id::Int`: index of wind turbine of interest (if any). If not a point for
    calculating effective wind speed of a turbine, then provide 0 (default)
"""
function point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set::AbstractModelSet;
                    wind_farm_state_id=1, downwind_turbine_id=0)
    
    # extract flow information
    wind_speed = wind_resource.wind_speeds[wind_farm_state_id]
    reference_height = wind_resource.measurement_heights[wind_farm_state_id]

    # set ground height 
    ground_height = wind_resource.wind_shear_model.ground_height    # TODO: allow topology to be given

    # find order for wind shear and deficit calculations
    shear_order = wind_resource.wind_shear_model.shear_order

    # adjust wind speed for wind shear
    if shear_order == "nothing"
        wind_speed_internal = wind_speed
    elseif shear_order == "first"
        wind_speed_internal = adjust_for_wind_shear(locz, wind_speed, reference_height, ground_height, wind_resource.wind_shear_model)
    else
        wind_speed_internal = wind_speed
    end

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize deficit summation term to zero
    deficit_sum = 0.0

    # loop through all turbines
    for u=1:nturbines

        # get index of upstream turbine
        upwind_turb_id = Int(sorted_turbine_index[u])

        # don't allow turbine to impact itself
        if upwind_turb_id == downwind_turbine_id; continue; end

        # downstream distance between upstream turbine and point
        x = locx - turbine_x[upwind_turb_id]

        # check turbine relative locations
        if x > 1E-6
            # skip this loop if it would include a turbine's impact on itself)
            if upwind_turb_id==downwind_turbine_id; continue; end

            # Check if tilt or yaw is being evaluated with
            if sum(turbine_yaw) == 0.0;
                vertical_deflection = wake_deflection_model(locx, locy, locz, turbine_x, turbine_tilt, turbine_ct,
                upwind_turb_id, rotor_diameter, turbine_local_ti, model_set.wake_deflection_model)

                horizontal_deflection = 0.0

                # velocity difference in the wake
                deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                turbine_local_ti, turbine_ct, turbine_tilt, model_set.wake_deficit_model)

            elseif sum(turbine_tilt) == 0.0;
                # calculate wake deflection of the current wake at the point of interest
                horizontal_deflection = wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct,
                                upwind_turb_id, rotor_diameter, turbine_local_ti, model_set.wake_deflection_model)

                vertical_deflection = 0.0

                # velocity difference in the wake
                deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                turbine_local_ti, turbine_ct, turbine_yaw, model_set.wake_deficit_model)
            else
                print("ERROR: This model only works for either Tilt or Yaw, not both.")
            end

            # combine deficits according to selected wake combination method
            deficit_sum = wake_combination_model(deltav, wind_speed_internal, wtvelocities[upwind_turb_id], deficit_sum, model_set.wake_combination_model)
            # println(deficit_sum, " ", downwind_turbine_id, " ", upwind_turb_id)
        end
    end

    # find velocity at point without shear
    point_velocity = wind_speed_internal - deficit_sum

    if shear_order == "nothing"
        point_velocity_out = point_velocity
    elseif shear_order == "first"
        point_velocity_out = point_velocity
    else
        point_velocity_out = adjust_for_wind_shear(locz, point_velocity, reference_height, ground_height, wind_resource.wind_shear_model)        
    end

    return point_velocity_out

end

"""
    turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set::AbstractModelSet; wind_farm_state_id::Int=1, velocity_only::Bool=true, turbine_velocities=nothing,
    turbine_ct=nothing, turbine_ai=nothing, turbine_local_ti=nothing)

# Arguments
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine sorted order upstream to downstream 
    for given state
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes 
    with state etc
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across 
    the rotor swept area when calculating the effective wind speed for the wind turbine. 
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades) 
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the 
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw, turbine_tilt,
                    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set::AbstractModelSet; wind_farm_state_id::Int=1, velocity_only::Bool=true, turbine_velocities=nothing,
                    turbine_ct=nothing, turbine_ai=nothing, turbine_local_ti=nothing)
    
    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)

    # initialize correct array types
    if turbine_velocities===nothing || turbine_ai === nothing || turbine_ct === nothing || turbine_local_ti === nothing
        arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                                typeof(hub_height[1]),typeof(turbine_yaw[1]))
        # initialize arrays
        if turbine_velocities === nothing
            turbine_velocities = zeros(arr_type, n_turbines)
        end
        if turbine_ct === nothing
            turbine_ct = zeros(arr_type, n_turbines)
        end
        if turbine_ai === nothing
            turbine_ai = zeros(arr_type, n_turbines)
        end
        if turbine_local_ti === nothing
            turbine_local_ti = zeros(arr_type, n_turbines)
        end
    end

    if typeof(model_set.wake_deficit_model) == CumulativeCurl{Float64,Vector{Float64}}
        turbine_velocities_one_direction_CC!(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
        sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
        model_set, turbine_velocities,
        turbine_ct, turbine_ai, turbine_local_ti; wind_farm_state_id=wind_farm_state_id, velocity_only=velocity_only)
    else
        turbine_velocities_one_direction!(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw, turbine_tilt,
        sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
        model_set, turbine_velocities,
        turbine_ct, turbine_ai, turbine_local_ti; wind_farm_state_id=wind_farm_state_id, velocity_only=velocity_only)
    end

    

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end

end

function turbine_velocities_one_direction!(turbine_x::Vector{T0}, turbine_y::Vector{T1}, turbine_z::Vector{T2}, rotor_diameter::Vector{T3}, hub_height::Vector{T4}, turbine_yaw::Vector{T5}, turbine_tilt::Vector{T6},
    sorted_turbine_index::Vector{Int}, ct_model::Vector{<:AbstractThrustCoefficientModel}, rotor_sample_points_y::Vector{T7}, rotor_sample_points_z::Vector{T7}, wind_resource,
    model_set::AbstractModelSet, turbine_velocities::Vector{T8},
    turbine_ct::Vector{T8}, turbine_ai::Vector{T8}, turbine_local_ti::Vector{T8}; wind_farm_state_id::Int=1, velocity_only::Bool=true) where {T0, T1, T2, T3, T4, T5, T6, T7, T8}

    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    # loop over all turbines
    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        wind_turbine_velocity = 0.0

        # initialize point vel with shear 
        point_velocity_with_shear = 0.0

        # loop over all rotor sample points to approximate the effective inflow velocity
        for p=1:n_rotor_sample_points

            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            # put rotor sample points in wind direction coordinate system, and account for tilt
            if sum(turbine_yaw) == 0.0
                locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_*sin(turbine_yaw[downwind_turbine_id])
                locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
                locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z
            # put rotor sample points in wind direction coordinate system, and account for yaw
            elseif sum(turbine_tilt) == 0.0;
                locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_z*sin(turbine_tilt[downwind_turbine_id])
                locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y
                locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z*cos(turbine_tilt[downwind_turbine_id])
            else
                print("ERROR: This model only works for either Tilt or Yaw, not both.")
            end

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                    wind_resource, model_set,
                                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id])

        # get local turbulence intensity for this wind state
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        
        # update local turbulence intensity for downstream turbine

        if sum(turbine_yaw) == 0.0
            turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_tilt, turbine_local_ti, sorted_turbine_index,
                                turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)
        elseif sum(turbine_tilt) == 0.0
            turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                                turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)
        end
    end

end

# Calculates the wind speeds for the farm with the Cumulative Curl model as defined in https://doi.org/10.5194/wes-2022-17
function turbine_velocities_one_direction_CC!(turbine_x::Vector{T0}, turbine_y::Vector{T1}, turbine_z::Vector{T2}, rotor_diameter::Vector{T3}, hub_height::Vector{T4}, turbine_yaw::Vector{T5},
    sorted_turbine_index::Vector{Int}, ct_model::Vector{<:AbstractThrustCoefficientModel}, rotor_sample_points_y::Vector{T6}, rotor_sample_points_z::Vector{T6}, wind_resource,
    model_set::AbstractModelSet, turbine_velocities::Vector{T7},
    turbine_ct::Vector{T7}, turbine_ai::Vector{T7}, turbine_local_ti::Vector{T7}; wind_farm_state_id::Int=1, velocity_only::Bool=true) where {T0, T1, T2, T3, T4, T5, T6, T7}

    @inbounds begin
        arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                                    typeof(hub_height[1]),typeof(turbine_yaw[1]))

        U_inf = wind_resource.wind_speeds[wind_farm_state_id]

        no_yaw = false
        if iszero(turbine_yaw)
            no_yaw = true
        end

        # get number of turbines, rotor sample points, and initialize contribution vector (see CC model)
        n_turbines = length(turbine_x)
        n_rotor_sample_points = length(rotor_sample_points_y)

        C = zeros(arr_type,n_turbines,n_turbines)
        sigma2 = zeros(arr_type,n_turbines,n_turbines)
        deflections = zeros(arr_type,n_turbines,n_turbines)

        deficits = zeros(arr_type,n_turbines,n_rotor_sample_points)
        point_velocities = zeros(arr_type,n_turbines,n_rotor_sample_points)
        point_velocities .= U_inf

        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]

        a_f = model_set.wake_deficit_model.a_f
        b_f = model_set.wake_deficit_model.b_f
        c_f = model_set.wake_deficit_model.c_f
        wec_factor = model_set.wake_deficit_model.wec_factor[1]

        zPos = zeros(n_turbines)

        for i = 1:n_turbines
            zPos[i] = turbine_z[i] + hub_height[i]
        end

        C_temp = zeros(arr_type,n_rotor_sample_points)

        # loop over all turbines (n)
        for n=1:n_turbines
            current_turbine_id = Int(sorted_turbine_index[n])
            x_n = turbine_x[current_turbine_id]
            y_n = turbine_y[current_turbine_id]
            z_n = zPos[current_turbine_id]

            # update current turbine velocity from sample points
            tot = 0
            for t = 1:n_rotor_sample_points
                tot += (point_velocities[current_turbine_id,t]^3)
            end
            turbine_velocities[current_turbine_id] = cbrt((tot / n_rotor_sample_points))

            # update coefficient of thrust for current turbine
            turbine_ct[current_turbine_id] = calculate_ct(turbine_velocities[current_turbine_id], ct_model[current_turbine_id])

            # update local TI for current turbine
            turbine_local_ti[current_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                                            turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=current_turbine_id, tol=1E-6)                     

            for d = n+1:n_turbines
                downwind_turbine_id = Int(sorted_turbine_index[d])
                for p = 1:n_rotor_sample_points
                    # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
                    local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
                    local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

                    # put rotor sample points in wind direction coordinate system, and account for yaw
                    x = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
                    y = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
                    z = zPos[downwind_turbine_id] + local_rotor_sample_point_z

                    # find order for wind shear and deficit calculations
                    shear_order = wind_resource.wind_shear_model.shear_order
                    # adjust wind speed for wind shear
                    if shear_order == "nothing"
                        wind_speed_internal = wind_speed
                    elseif shear_order == "first"
                        wind_speed_internal = adjust_for_wind_shear(z, wind_resource.wind_speeds[wind_farm_state_id], wind_resource.measurement_heights[wind_farm_state_id], wind_resource.wind_shear_model.ground_height, wind_resource.wind_shear_model)
                    else
                        wind_speed_internal = wind_speed
                    end

                    # CC model
                    @fastmath x_tilde_n = abs(x - x_n) / rotor_diameter[current_turbine_id]
                    @fastmath m = a_f*exp(b_f*x_tilde_n)+c_f
                    @fastmath a1 = 2^(2/m - 1)
                    @fastmath a2 = a1^2

                    if p == 1
                        if no_yaw == false
                            deflections[current_turbine_id,downwind_turbine_id] = wake_deflection_model(x, y, z, turbine_x, turbine_yaw, turbine_ct, current_turbine_id, 
                                    rotor_diameter, turbine_local_ti, model_set.wake_deflection_model)
                        end
                        sigma2[current_turbine_id,downwind_turbine_id] = wake_expansion(turbine_ct[current_turbine_id],turbine_local_ti[current_turbine_id],x_tilde_n,model_set.wake_deficit_model)
                    end

                    dy = deflections[current_turbine_id,downwind_turbine_id]
                    sigma_n = sigma2[current_turbine_id,downwind_turbine_id]
                    sum_C = 0.0

                    for i = 1:n-1
                        other_turbine_id = Int(sorted_turbine_index[i])
                        y_i = turbine_y[other_turbine_id]
                        z_i = zPos[other_turbine_id]
                        sigma_i = sigma2[other_turbine_id,downwind_turbine_id]
                        dy_i = deflections[other_turbine_id,downwind_turbine_id]
                        @fastmath lambda_n_i = sigma_n/(sigma_n+sigma_i) * exp(-((y_n-y_i-dy_i)^2 + (z_n-z_i)^2)/(2.0*(sigma_n+sigma_i)))
                        @fastmath sum_C += lambda_n_i * C[other_turbine_id,downwind_turbine_id]
                    end

                    @fastmath calc = abs_smooth(a2 - (m*turbine_ct[current_turbine_id]*cos(turbine_yaw[current_turbine_id]))/(16.0*gamma(2/m)*(sigma_n^(2/m))*(1-sum_C/U_inf)^2),0.1)
                    @fastmath C_point = (1-sum_C/U_inf) * (a1-sqrt(calc))
                    @fastmath r_tilde = (sqrt((y-y_n-dy)^2 + (z-z_n)^2)/rotor_diameter[current_turbine_id])
                    @fastmath velDef = C_point*exp(-1 * (r_tilde^m)/(2.0*sigma_n*wec_factor))
                    @fastmath deficits[downwind_turbine_id,p] += velDef * turbine_velocities[current_turbine_id]
                    if d == n+1
                        point_velocities[downwind_turbine_id,p] = U_inf - deficits[downwind_turbine_id,p]
                        # find order for wind shear and deficit calculations
                        shear_order = wind_resource.wind_shear_model.shear_order
                        # adjust wind speed for wind shear
                        if shear_order == "nothing"
                            point_velocities[downwind_turbine_id,p] = point_velocities[downwind_turbine_id,p]
                        elseif shear_order == "first"
                            point_velocities[downwind_turbine_id,p] = point_velocities[downwind_turbine_id,p]
                        else
                            point_velocities[downwind_turbine_id,p] = adjust_for_wind_shear(z, point_velocities[downwind_turbine_id,p], wind_resource.measurement_heights[wind_farm_state_id], wind_resource.wind_shear_model.ground_height, wind_resource.wind_shear_model)
                        end
                    end
                    C_temp[p] = C_point
                    C_avg = 0
                    for s = 1:p
                        C_avg += C_temp[s]
                    end
                    C_avg /= p

                    C[current_turbine_id,downwind_turbine_id] = C_avg
                end
            end
        end
    end
end

# Helper function for turbine_velocities_one_direction_CC!
function wake_expansion(Ct,TI,x_tilde,model)
    @fastmath beta = 0.5*(1.0+sqrt(1.0-Ct))/sqrt(1.0-Ct)
    epsilon = (model.c_s1*Ct+model.c_s2)*sqrt(beta)
    k = (model.a_s*TI+model.b_s)
    sigma = k*x_tilde+epsilon

    return sigma^2
end

# function turbine_velocities_one_direction(x, turbine_z, rotor_diameter, hub_height, turbine_yaw,
#     sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
#     model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)

#     n_turbines = Int(length(x)/2)
#     # println(typeof(x), n_turbines)
#     turbine_x = x[1:n_turbines] 
#     turbine_y = x[n_turbines+1:end]
#     # println(turbine_x)
#     # println("turbine_x type ", typeof(turbine_x))
#     # println("type of x ", typeof(x))

#     # get number of turbines and rotor sample point
#     # n_turbines = length(turbine_x)
#     n_rotor_sample_points = length(rotor_sample_points_y)

#     arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
#                 typeof(hub_height[1]),typeof(turbine_yaw[1]))
#     turbine_velocities = zeros(arr_type, n_turbines)
#     turbine_ct = zeros(arr_type, n_turbines)
#     turbine_ai = zeros(arr_type, n_turbines)
#     turbine_local_ti = zeros(arr_type, n_turbines)

#     for d=1:n_turbines

#         # get index of downstream turbine
#         downwind_turbine_id = Int(sorted_turbine_index[d])

#         # initialize downstream wind turbine velocity to zero
#         # println("start array: ", turbine_velocities[downwind_turbine_id])
#         # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
#         wind_turbine_velocity = 0.0
#         # turbine_velocities[downwind_turbine_id] = 0.0

#         for p=1:n_rotor_sample_points


#             # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
#             local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
#             local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

#             locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
#             locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
#             locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

#             # calculate the velocity at given point
#             point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
#                                 rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
#                                 wind_resource, model_set,
#                                 wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

#             # add sample point velocity to turbine velocity to be averaged later
#             wind_turbine_velocity += point_velocity_with_shear

#         end

#         # final velocity calculation for downstream turbine (average equally across all points)
#         wind_turbine_velocity /= n_rotor_sample_points

#         turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

#         # update thrust coefficient for downstream turbine
#         turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

#         # update axial induction for downstream turbine
#         turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id])

#         # update local turbulence intensity for downstream turbine
#         ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
#         turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
#                     turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

#     end

#     if velocity_only
#         return turbine_velocities 
#     else
#         return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
#     end

# end

# turbine_velocities_one_direction!(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1) = turbine_velocities_one_direction!([0.0], [0.0],
# model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1)

"""
calculate_flow_field(xrange, yrange, zrange, model_set::AbstractModelSet, turbine_x, 
    turbine_y, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai, rotor_diameter, hub_height, 
    turbine_local_ti, sorted_turbine_index, wtvelocities, wind_resource; wind_farm_state_id=1)

Generates a flow field for a given state and cross section

# Arguments
- `xrange::Range`: range defining east-west locations to sample in global reference frame
- `yrange::Range`: range defining north-west locations to sample in global reference frame
- `zrange::Range`: range defining vertical locations to sample in global reference frame
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the global reference frame
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_tilt::Array{TF,nTurbines}`: turbine tilt for the given wind direction in 
radians
- `turbine_ct::Array{TF,nTurbines}`: thrust coefficient of each turbine for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource; wind_farm_state_id=1)

    xlen = length(xrange)
    ylen = length(yrange)
    zlen = length(zrange)
    npoints = xlen*ylen*zlen
    point_velocities = zeros(npoints)
    point_velocities = reshape(point_velocities, (zlen, ylen, xlen))

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                locx = xrange[xi]
                locy = yrange[yi]
                locz = zrange[zi]
                locx, locy = rotate_to_wind_direction(locx, locy, wind_resource.wind_directions[wind_farm_state_id])

                point_velocities[zi, yi, xi] = point_velocity(locx, locy, locz, rot_tx, rot_ty, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set,
                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=0)

            end
        end
    end

    

    # if zlen == 1
    #     return point_velocities[1,1:ylen,1:xlen]
    # elseif ylen == 1
    #     return point_velocities[1:zlen,1,1:xlen]
    # elseif xlen == 1
    #     return point_velocities[1:zlen,1:ylen,1]
    # else
    return point_velocities[1:zlen,1:ylen,1:xlen]
    # end

end

function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_tilt,
    rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource; wind_farm_state_id=1)

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = turbine_velocities_one_direction(rot_tx, rot_ty, turbine_z, rotor_diameter, hub_height, turbine_yaw, turbine_tilt,
    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set, wind_farm_state_id=wind_farm_state_id, velocity_only=false)

    return calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_tilt, turbine_ct, turbine_ai,
        rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
        wind_resource, wind_farm_state_id=wind_farm_state_id)

end