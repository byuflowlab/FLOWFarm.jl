export WindFarmModelSet
abstract type AbstractModelSet end

"""
    WindFarmModelSet(wakedeficitmodel, wake_deflection_model, wake_combination_model, local_ti_model)

Container for objects defining models to use in wind farm calculations

# Arguments
- `wake_defiict_model::AbstractWakeDeficitModel`: contains a struct defining the desired wake deficit model
- `wake_deflection_model::AbstractWakeDeflectionModel`: contains a struct defining the desired wake deflection model
- `wake_combination_model::AbstractWakeCombinationModel`: contains a struct defining the desired wake combination model
- `local_ti_model::AbstractTurbulenceIntensityModel`: contains a struct defining the desired turbulence intensity model
- `point_velocity_average_factor::Float`: factor used to determine how point velocity is averaged across the turbine (3 would result in a cubic mean)
"""
struct WindFarmModelSet{DTM,DNM,CM,TIM,F} <: AbstractModelSet

    wake_deficit_model::DTM
    wake_deflection_model::DNM
    wake_combination_model::CM
    local_ti_model::TIM
    point_velocity_average_factor::F

end

WindFarmModelSet(a,b,c,d) = WindFarmModelSet(a,b,c,d,1.)

"""
    point_velocity(loc, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
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
function point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set::AbstractModelSet, wake_deficits, contribution_matrix, deflections, sigma_squared;
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
    deficit_sum = eltype(hub_height)(0.0)

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

            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct,
                            upwind_turb_id, rotor_diameter, turbine_local_ti, model_set.wake_deflection_model)

            if !isempty(deflections)
                deflections[upwind_turb_id, downwind_turbine_id] = horizontal_deflection
            end

            vertical_deflection = 0.0

            # velocity difference in the wake
            deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                            upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                            turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, u, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model_set.wake_deficit_model)

            if !isempty(wake_deficits)
                wake_deficits[upwind_turb_id,downwind_turbine_id] = deltav
            end

            # combine deficits according to selected wake combination method
            deficit_sum = wake_combination_model(deltav, wind_speed_internal, wtvelocities[upwind_turb_id], deficit_sum, model_set.wake_combination_model)
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
function turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set::AbstractModelSet; wind_farm_state_id::Int=1, velocity_only::Bool=true, turbine_velocities=nothing,
                    turbine_ct=nothing, turbine_ai=nothing, turbine_local_ti=nothing, using_sparsity::Bool=false,
                    wake_deficits=nothing, contribution_matrix=nothing, deflections=nothing, sigma_squared=nothing)

    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)

    # initialize correct array types
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                                typeof(hub_height[1]),typeof(turbine_yaw[1]))

    # initialize arrays
    if turbine_velocities === nothing
        turbine_velocities = zeros(arr_type, n_turbines)
    else
        turbine_velocities .= 0
    end
    if turbine_ct === nothing
        turbine_ct = zeros(arr_type, n_turbines)
    else
        turbine_velocities .= 0
    end
    if turbine_ai === nothing
        turbine_ai = zeros(arr_type, n_turbines)
    else
        turbine_ai .= 0
    end
    if turbine_local_ti === nothing
        turbine_local_ti = zeros(arr_type, n_turbines)
    else
        turbine_local_ti .= 0
    end
    if wake_deficits === nothing
        wake_deficits = zeros(arr_type,n_turbines,n_turbines)
    else
        wake_deficits .= 0
    end
    if contribution_matrix === nothing
        contribution_matrix = zeros(arr_type,n_turbines,n_turbines)
    else
        contribution_matrix .= 0
    end
    if deflections === nothing
        deflections = zeros(arr_type,n_turbines,n_turbines)
    else
        deflections .= 0
    end
    if sigma_squared === nothing
        sigma_squared = zeros(arr_type,n_turbines,n_turbines)
    else
        sigma_squared .= 0
    end

    turbine_velocities_one_direction!(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
        sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource, model_set, turbine_velocities,
        turbine_ct, turbine_ai, turbine_local_ti, wake_deficits, contribution_matrix, deflections, sigma_squared;
        wind_farm_state_id=wind_farm_state_id, velocity_only=velocity_only)

    wake_deficits .= wake_deficits'

    if using_sparsity
        return turbine_velocities, wake_deficits
    elseif velocity_only
        return turbine_velocities
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end

end

function turbine_velocities_one_direction!(turbine_x::T0, turbine_y::T1, turbine_z::T2, rotor_diameter::T3, hub_height::T4, turbine_yaw::T5,
    sorted_turbine_index::Vector{Int}, ct_model::Vector{<:AbstractThrustCoefficientModel}, rotor_sample_points_y::Vector{T6}, rotor_sample_points_z::Vector{T6}, wind_resource,
    model_set::AbstractModelSet, turbine_velocities::T7,
    turbine_ct::T7, turbine_ai::T7, turbine_local_ti::T7, wake_deficits::T8, contribution_matrix::T9, deflections::T10, sigma_squared::T11; wind_farm_state_id::Int=1, velocity_only::Bool=true) where {T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11}

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

            # put rotor sample points in wind direction coordinate system, and account for yaw
            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                    wind_resource, model_set, wake_deficits, contribution_matrix, deflections, sigma_squared,
                                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear^model_set.point_velocity_average_factor
        end


        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        wind_turbine_velocity = (wind_turbine_velocity)^(1.0/model_set.point_velocity_average_factor)

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id])

        # get local turbulence intensity for this wind state
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]

        # update local turbulence intensity for downstream turbine
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                            turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)
    end

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
    turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, rotor_diameter, hub_height,
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
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
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

    wake_deficits = []
    contribution_matrix = []
    deflections = []
    sigma_squared = []

    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                locx = xrange[xi]
                locy = yrange[yi]
                locz = zrange[zi]
                locx, locy = rotate_to_wind_direction(locx, locy, wind_resource.wind_directions[wind_farm_state_id])

                point_velocities[zi, yi, xi] = point_velocity(locx, locy, locz, rot_tx, rot_ty, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set, wake_deficits, contribution_matrix, deflections, sigma_squared,
                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=0)

            end
        end
    end

    return point_velocities[1:zlen,1:ylen,1:xlen]

end

function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw,
    rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource; wind_farm_state_id=1)

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = turbine_velocities_one_direction(rot_tx, rot_ty, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set, wind_farm_state_id=wind_farm_state_id, velocity_only=false)

    return calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
        rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
        wind_resource, wind_farm_state_id=wind_farm_state_id)

end
