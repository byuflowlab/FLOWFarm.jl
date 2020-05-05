abstract type AbstractWindFarmProblem end
abstract type AbstractModelSet end

"""
    WindFarmProblemDescription(windfarm, windresource, windfarmstates)

Container for objects defining a wind farm problem

# Arguments
- `wind_farm::AbstractWindFarm`: contains windturbine coordinates and definitions
- `wind_resource::AbstracWindResource`: wind resource description
- `wind_farm_states::Array{SingleWindFarmState}(Nstates)`: contains turbine coordinates operational states
"""
struct WindFarmProblemDescription{FM,WR,AFS} <: AbstractWindFarmProblem

    wind_farm::FM
    wind_resource::WR
    wind_farm_states::AFS

end


"""
    WindFarmModelSet(wakedeficitmodel, wake_deflection_model, wake_combination_model, ti_model, wind_shear_model)

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

function rotate_to_wind_direction(xlocs, ylocs, wind_direction_met)
    # use radians

    # convert from meteorological polar system (CW, 0 rad.=N) to standard polar system (CCW, 0 rad.=E)
    wind_direction_cart = (3*pi/2 - wind_direction_met)

    if wind_direction_cart < 0.0
        wind_direction_cart += 2*pi
    end

    cos_wdr = cos(-wind_direction_cart)
    sin_wdr = sin(-wind_direction_cart)

    # convert to cartesian coordinates with wind to positive x
    x_cart = xlocs*cos_wdr - ylocs*sin_wdr
    y_cart = xlocs*sin_wdr + ylocs*cos_wdr

    return x_cart, y_cart
end

function point_velocity(loc, model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem;
    wind_farm_state_id=1, downwind_turbine_id=0)

    wakedeficitmodel = model_set.wake_deficit_model
    wakedeflectionmodel = model_set.wake_deflection_model
    wakecombinationmodel = model_set.wake_combination_model

    windfarm = problem_description.wind_farm
    windresource = problem_description.wind_resource
    windfarmstate = problem_description.wind_farm_states[wind_farm_state_id]

    # extract turbine locations in rotated reference frame
    turbine_x = windfarmstate.turbine_x
    turbine_y = windfarmstate.turbine_y
    turbine_z = windfarm.turbine_z
    turbine_definition_ids = windfarm.turbine_definition_ids

    # extract turbine definitions
    turbine_definitions = windfarm.turbine_definitions

    # get sorted wind turbine index in currect direction
    sorted_turbine_index = windfarmstate.sorted_turbine_index

    # get current inflow velocities at each turbine
    wtvelocities = windfarmstate.turbine_inflow_velcities

    # extract flow information
    wind_speed = windresource.wind_speeds[wind_farm_state_id]
    reference_height = windresource.measurement_heights[wind_farm_state_id]
    wind_shear_model = windresource.wind_shear_model[1]
    shear_exponent = wind_shear_model.shear_exponent

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize deficit summation term to zero
    deficit_sum = 0.0

    # initialize point velocity with shear to zero
    point_velocity_with_shear = 0.0

    # loop through all turbines
    for u=1:nturbines

        # get index of upstream turbine
        upwind_turb_id = sorted_turbine_index[u]

        # skip this loop if it would include a turbine's impact on itself)
        if upwind_turb_id==downwind_turbine_id; continue; end

        # get turbine definition
        upwind_turbine = turbine_definitions[turbine_definition_ids[upwind_turb_id]]

        # downstream distance between upstream turbine and point
        x = loc[1] - turbine_x[upwind_turb_id]

        # set this iterations velocity deficit to 0
        deltav = 0.0

        # check turbine relative locations
        if x > 0.0

            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = wake_deflection_model(loc, upwind_turb_id, upwind_turbine, wakedeflectionmodel, windfarmstate)
            vertical_deflection = 0.0

            # velocity difference in the wake
            deltav = wake_deficit_model(loc, [horizontal_deflection, vertical_deflection], upwind_turb_id, upwind_turbine, wakedeficitmodel, windfarmstate)

            # combine deficits according to selected wake combination method
            turb_inflow = wtvelocities[upwind_turb_id]
            deficit_sum = wake_combination_model(deltav, wind_speed, turb_inflow, deficit_sum, wakecombinationmodel)

        end

        # find velocity at point without shear
        point_velocity_without_shear = wind_speed - deficit_sum

        # adjust sample point velocity for shear
        point_velocity_with_shear = adjust_for_wind_shear(loc, point_velocity_without_shear, reference_height, turbine_z[upwind_turb_id], wind_shear_model)

    end

    return point_velocity_with_shear

end


function turbine_velocities_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
    model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1)

    windfarm = problem_description.wind_farm
    windfarmstate = problem_description.wind_farm_states[wind_farm_state_id]

    # get number of turbines and rotor sample point
    n_turbines = length(windfarmstate.turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = windfarmstate.sorted_turbine_index[d]

        # get turbine definition of downstream turbine
        turbine_definition_id = windfarm.turbine_definition_ids[downwind_turbine_id]
        downwind_turbine = windfarm.turbine_definitions[turbine_definition_id]

        # initialize downstream wind turbine velocity to zero
        wind_turbine_velocity = 0.0

        for p=1:n_rotor_sample_points

            loc = zeros(3)
            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*downwind_turbine.rotor_diameter[1]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*downwind_turbine.rotor_diameter[1]

            # move sample points to correct height and yaw location in wind farm state reference frame
            loc[1] = windfarmstate.turbine_x[downwind_turbine_id] + local_rotor_sample_point_y*sin(windfarmstate.turbine_yaw[downwind_turbine_id])
            loc[2] = windfarmstate.turbine_y[downwind_turbine_id] + local_rotor_sample_point_y*cos(windfarmstate.turbine_yaw[downwind_turbine_id])
            loc[3] = windfarmstate.turbine_z[downwind_turbine_id] + downwind_turbine.hub_height[1] + local_rotor_sample_point_z

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(loc, model_set, problem_description,
                wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        # update wind farm state with new velocity for downstream turbine
        windfarmstate.turbine_inflow_velcities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        windfarmstate.turbine_ct[downwind_turbine_id] = calculate_ct(wind_turbine_velocity, downwind_turbine.ct_model[1])

        # update local turbulence intensity for downstream turbine
        ambient_ti = problem_description.wind_resource.ambient_tis[wind_farm_state_id]
        windfarmstate.turbine_local_ti[downwind_turbine_id] = calculate_local_ti(ambient_ti, model_set.local_ti_model)

    end

end

turbine_velocities_one_direction!(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1) = turbine_velocities_one_direction!([0.0], [0.0],
model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1)


function calculate_flow_field(direction_id, xrange, yrange, zrange, rotor_sample_points_y, rotor_sample_points_z,
    model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem;
    wind_farm_state_id=1)

    windresource = problem_description.wind_resource

    xlen = length(xrange)
    ylen = length(yrange)
    zlen = length(zrange)
    npoints = xlen*ylen*zlen
    point_velocities = zeros(npoints)
    point_velocities = reshape(point_velocities, (zlen, ylen, xlen))

    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                loc = [xrange[xi], yrange[yi], zrange[zi]]
                loc[1], loc[2] = rotate_to_wind_direction(loc[1], loc[2], windresource.wind_directions[direction_id])

                point_velocities[zi, yi, xi] = point_velocity(loc, model_set, problem_description, wind_farm_state_id=wind_farm_state_id)

            end
        end
    end

    if zlen == 1
        return point_velocities[1,1:ylen,1:xlen]
    elseif ylen == 1
        return point_velocities[1:zlen,1,1:xlen]
    elseif xlen == 1
        return point_velocities[1:zlen,1:ylen,1]
    end

end


function hermite_spline(x, x0, x1, y0, dy0, y1, dy1)
       """This function produces the y and dy values for a hermite cubic spline
       interpolating between two end points with known slopes

       :param x: x position of output y
       :param x0: x position of upwind endpoint of spline
       :param x1: x position of downwind endpoint of spline
       :param y0: y position of upwind endpoint of spline
       :param dy0: slope at upwind endpoint of spline
       :param y1: y position of downwind endpoint of spline
       :param dy1: slope at downwind endpoint of spline

       :return: y: y value of spline at location x"""

    # initialize coefficients for parametric cubic spline
    c3 = (2.0*(y1))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) - \
         (2.0*(y0))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) + \
         (dy0)/(x0^2 - 2.0*x0*x1 + x1^2) + \
         (dy1)/(x0^2 - 2.0*x0*x1 + x1^2)

    c2 = (3.0*(y0)*(x0 + x1))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) - \
         ((dy1)*(2.0*x0 + x1))/(x0^2 - 2.0*x0*x1 + x1^2) - ((dy0)*(x0 +
         2.0*x1))/(x0^2 - 2.0*x0*x1 + x1^2) - (3.0*(y1)*(x0 + x1))/(x0^3 -
         3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3)

    c1 = ((dy0)*(x1^2 + 2.0*x0*x1))/(x0^2 - 2.0*x0*x1 + x1^2) + ((dy1)*(x0^2 +
         2.0*x1*x0))/(x0^2 - 2.0*x0*x1 + x1^2) - (6.0*x0*x1*(y0))/(x0^3 -
         3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) + (6.0*x0*x1*(y1))/(x0^3 -
         3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3)

    c0 = ((y0)*(- x1^3 + 3.0*x0*x1^2))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 -
         x1^3) - ((y1)*(- x0^3 + 3.0*x1*x0^2))/(x0^3 - 3.0*x0^2*x1 +
         3.0*x0*x1^2 - x1^3) - (x0*x1^2*(dy0))/(x0^2 - 2.0*x0*x1 + x1^2) - \
         (x0^2*x1*(dy1))/(x0^2 - 2.0*x0*x1 + x1^2)

    # Solve for y and dy values at the given point
    y = c3*x^3 + c2*x^2 + c1*x + c0
    # dy_dx = c3*3*x^2 + c2*2*x + c1

    # return y, dy_dx
    return y
end
