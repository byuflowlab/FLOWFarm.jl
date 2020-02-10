# include("wake_models.jl")
# include("deflection_models.jl")

# abstract type AbstractWindFarmModel end
# abstract type AbstractWindResourceModel end

# struct WindFarm{AF} <: AbstractWindFarmModel
    
#     turbinex::AF
#     turbiney::AF
#     turbinez::AF
#     hub_height::AF

# end

# struct DiscretizedWindResource{AF} <: AbstractWindResourceModel
    
#     winddirections::AF
#     windspeeds::AF
#     windpropabilities::AF

# end

# function _wind_frame(xlocs, ylocs, wind_direction)

#     nTurbines = length(xlocs)
#     wd = -pi*(270. - wind_direction)/180. #shift to traditional wind direction coords and to radians
#     xw = xlocs.*cos(wd)-ylocs.*sin(wd)
#     yw = xlocs.*sin(wd)+ylocs.*cos(wd)
#     return xw, yw

# end

# function _point_velocity(loc, turb, turbI, sortedturbinexindex, wakemodel::AbstractWakeModel, windfarm::WindFarm, deflectionmodel::AbstractDeflectionModel)

#     # extract turbine locations
#     turbinex = windfarm.turbinex
#     turbineY = windfarm.turbiney

#     # get number of turbines
#     nturbines = length(turbinex)

#     # get hub height and vertical location
#     hub_height = windfarm.hub_height
#     turbinez = windfarm.turbinez

#     # extract turbine types
#     turbines = windfarm.turbines

#     # initialize deficit summation term to zero
#     deficitsum = 0.0
    
#     # loop through all turbines
#     for u=1:nTurbines 
        
#         # get index of upstream turbine
#         turb = sortedturbinexindex[u] + 1.0
        
#         # skip this loop if turb = turbI (turbines impact on itself)
#         if turb==turbI; continue; end
    
#         # downstream distance between upstream turbine and point
#         x = loc[1] - turbinexw[turb]
    
#         # set this iterations velocity deficit to 0
#         deltav = 0.0
        
#         # check turbine relative locations
#         if x > 0.0
        
#             # calculate wake deflection of the current wake at the point of interest
#             horizontal_deflection = deflection_model(loc, wake_model, turbines[turb])
            
#             # cross wind distance from point location to upstream turbine wake center
#             deltay = pointY - (turbineyw[turb] + horizontal_deflection)

#             # vertical distance from upstream hub height to height of point of interest
#             deltaz = pointZ - (turbinez + hub_height)

#             #TODO left off here


#     #         if (x > x0) then
#     #             ! velocity difference in the wake
#     #             call deltav_func(deltay, deltaz, Ct_local(turb), yaw(turb), &
#     #                             & sigmay, sigmaz, rotorDiameter(turb), & 
#     #                             & wake_model_version, kz_local(turb), x, &
#     #                             & wec_factor, sigmay_spread, sigmaz_spread, deltav)
                                
#     #         else
#     #             ! velocity deficit in the nearwake (linear model)
#     #             call deltav_near_wake_lin_func(deltay, deltaz, &
#     #                              & Ct_local(turb), yaw(turb), sigmay_0, sigmaz_0, x0, & 
#     #                              & rotorDiameter(turb), x, discontinuity_point, & 
#     #                              & sigmay_d, sigmaz_d, wake_model_version, &
#     #                              & kz_local(turb), x0, sigmay_spread, &
#     #                              & sigmaz_spread, sigmay_0_spread, &
#     #                              & sigmaz_0_spread, wec_factor, WECH, deltav)
                                 
#     #         end if

#     #         ! save deficit sum in holder for AD purposes
#     #         old_deficit_sum = deficit_sum
            
#     #         ! combine deficits according to selected wake combination method
#     #         call wake_combination_func(wind_speed, wtVelocity(turb), deltav,         &
#     #                                    wake_combination_method, old_deficit_sum, deficit_sum)
        
#     #     end if
    
#     # end do
    
#     # ! find velocity at point without shear
#     # point_velocity = wind_speed - deficit_sum
    
#     # ! adjust sample point velocity for shear
#     # call wind_shear_func(pointZ, point_velocity, z_ref, z_0, shear_exp, point_velocity_with_shear)
#     # 
#     end
# # end subroutine point_velocity_with_shear_func



# end
