abstract type AbstractWindShearModel end

struct PowerLawWindShear{TF} <: AbstractWindShearModel
    
    # model parameter
    shear_exponent::TF

end

#TODO add log shear

function adjust_for_wind_shear(loc, reference_velocity, reference_height, ground_height, model::PowerLawWindShear)

    # initialize adjusted wind speed to zero
    adjusted_wind_speed = 0.0
    shear_exp = model.shear_exponent[1]

    # check that the point of interest is above ground level
    if loc[3] >= ground_height
        # adjusted wind speed for wind shear if point is above ground
        adjusted_wind_speed = reference_velocity*((loc[3]-ground_height)/(reference_height-ground_height))^shear_exp
    end

    return adjusted_wind_speed
end

# function adjust_for_wind_shear(loc, point_velocity_no_shear, reference_height, ground_height, shear_model::CubicWindShear)

#     # initialize adjusted wind speed to zero
#     adjusted_wind_speed = 0.0

#     # check that the point of interest is above ground level
#     if loc[3] >= ground_height
#         # adjusted wind speed for wind shear if point is above ground
#         adjusted_wind_speed = point_velocity_no_shear*((loc[3]-ground_height)/(reference_height-ground_height))^shear_exp
#     else 
#         # if the point of interest is below ground, set the wind speed to 0.0
#         adjusted_wind_speed = 0.0
#     end

#     return adjusted_wind_speed
# end