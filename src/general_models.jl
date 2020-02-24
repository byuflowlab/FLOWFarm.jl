function rotate_to_wind_direction(xlocs, ylocs, wind_direction)

    nTurbines = length(xlocs)
    wd = -pi*(270. - wind_direction)/180. #shift to traditional wind direction coords and to radians
    xw = xlocs.*cos(wd)-ylocs.*sin(wd)
    yw = xlocs.*sin(wd)+ylocs.*cos(wd)
    return xw, yw

end

function point_velocity(loc,
    windfarm::WindFarm,
    windfarmstate::SingleWindFarmState,
    windresource::AbstractWindResourceModel,
    wakedeficitmodel::AbstractWakeDeficitModel, 
    wakedeflectionmodel::AbstractWakeDeflectionModel, 
    wakecombinationmodel::AbstractWakeCombinationModel,
    windshearmodel::AbstractWindShearModel,
    turbine_id=0)

    # get state id
    direction_id = windfarmstate.id

    # extract turbine locations in rotated reference frame
    turbine_x = windfarmstate.turbine_x
    turbine_y = windfarmstate.turbine_y
    turbine_z = windfarm.turbine_z

    # extract turbine definitions
    turbines = windfarm.turbine_definitions

    # get sorted wind turbine index in currect direction
    sorted_turbine_index = windfarmstate.sorted_turbine_index

    # get current inflow velocities at each turbine
    wtvelocities = windfarmstate.turbine_inflow_velcities

    # extract flow information
    wind_speed = windresource.wind_speeds[direction_id]
    reference_height = windresource.measurement_heights[direction_id]
    shear_exponent = windresource.shear_exponent

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize deficit summation term to zero
    deficit_sum = 0.0

    # initialize point velocity with shear to zero
    point_velocity_with_shear = 0.0
    
    # loop through all turbines
    for u=1:nturbines 
        
        # get index of upstream turbine
        turb = sorted_turbine_index[u]
        
        # skip this loop if it would include a turbine's impact on itself)
        if turb==turbine_id; continue; end

        # get turbine definition
        turbine = turbines[turb]
        
        # downstream distance between upstream turbine and point
        x = loc[1] - turbine_x[turb]
    
        # set this iterations velocity deficit to 0
        deltav = 0.0
        
        # check turbine relative locations
        if x > 0.0
        
            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = wake_deflection_model(loc, wakedeflectionmodel, turbine, windfarmstate)
            vertical_deflection = 0.0
            
            # velocity difference in the wake
            deltav = wake_deficit_model(loc, [horizontal_deflection, vertical_deflection], wakedeficitmodel, turbine, windfarmstate)
            
            # combine deficits according to selected wake combination method 
            turb_inflow = wtvelocities[turb]
            deficit_sum += wake_combination_model(deltav, wind_speed, turb_inflow, deficit_sum, wakecombinationmodel)

        end
        
        # find velocity at point without shear
        point_velocity_without_shear = wind_speed - deficit_sum
        
        # adjust sample point velocity for shear
        point_velocity_with_shear = adjust_for_wind_shear(loc, point_velocity_without_shear, reference_height, turbine_z[turb], windshearmodel)

    end

    return point_velocity_with_shear

end

function turbine_velocities_one_direction!(rotor_sample_points_y, rotor_sample_points_z, 
    windfarm::AbstractWindFarmModel, 
    windfarmstate::SingleWindFarmState,
    windresource::AbstractWindResourceModel, 
    windshearmodel::AbstractWindShearModel, 
    wakedeficitmodel::AbstractWakeDeficitModel, 
    wakedeflectionmodel::AbstractWakeDeflectionModel, 
    wakecombinationmodel::AbstractWakeCombinationModel)

    # get number of turbines and rotor sample point
    n_turbines = length(windfarmstate.turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    for d=1:n_turbines
    
        # get index of downstream turbine
        downstream_turb_index = windfarmstate.sorted_turbine_index[d]

        # get turbine definition of downstream turbine
        downstream_turbine = windfarm.turbine_definitions[downstream_turb_index]

        # initialize downstream wind turbine velocity to zero
        wind_turbine_velocity = 0.0

        for p=1:n_rotor_sample_points
        
            loc = [0.0 0.0 0.0]
            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*downstream_turbine.rotor_diameter
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*downstream_turbine.rotor_diameter
            
            # move sample points to correct height and yaw location in wind farm state reference frame
            loc[1] = windfarmstate.turbine_x[downstream_turb_index] + local_rotor_sample_point_y*sin(windfarmstate.turbine_yaw[downstream_turb_index]) 
            loc[2] = windfarmstate.turbine_y[downstream_turb_index] + local_rotor_sample_point_y*cos(windfarmstate.turbine_yaw[downstream_turb_index]) 
            loc[3] = windfarmstate.turbine_z[downstream_turb_index] + downstream_turbine.hub_height + local_rotor_sample_point_z
            
            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(loc,
                windfarm::WindFarm,
                windfarmstate::SingleWindFarmState,
                windresource::AbstractWindResourceModel,
                wakedeficitmodel::AbstractWakeDeficitModel, 
                wakedeflectionmodel::AbstractWakeDeflectionModel, 
                wakecombinationmodel::AbstractWakeCombinationModel,
                windshearmodel::AbstractWindShearModel,
                downstream_turb_index)
            
            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear
        
        end
    
        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        # update wind farm state with new velocity for downstream turbine
        windfarmstate.turbine_inflow_velcities[downstream_turb_index] = wind_turbine_velocity
        
        # update thrust coefficient for downstream turbine
        windfarmstate.turbine_ct[downstream_turb_index] = calculate_ct(downstream_turbine.ct_model)

        # TODO add local turbulence intensity calculations

    end

end