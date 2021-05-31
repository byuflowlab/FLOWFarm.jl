abstract type AbstractWindShearModel end

struct PowerLawWindShear{TF} <: AbstractWindShearModel
    
    # model parameter
    shear_exponent::TF

end

#TODO add log shear
"""
    adjust_for_wind_shear(locz, reference_velocity, reference_height, ground_height, model::PowerLawWindShear)

Uses provided velocity at a given height to estimate the velocity at
a different height due to wind shear.

# Arguments
- `locz::Float`: height of desired velocity
- `reference_velocity::Float`: known velocity at reference_height
- `reference_height::Float`: height of known velocity 
- `ground_height::Float`: height of the ground (typically zero)
- `model::AbstractWindShearModel`: wind shear model to use for calculations
"""
function adjust_for_wind_shear(locz, reference_velocity, reference_height, ground_height, model::PowerLawWindShear)

    # initialize adjusted wind speed to zero
    adjusted_wind_speed = 0.0
    shear_exp = model.shear_exponent
    
    # check that the point of interest is above ground level
    if locz >= ground_height
        # adjusted wind speed for wind shear if point is above ground
        adjusted_wind_speed = reference_velocity*((locz-ground_height)/(reference_height-ground_height))^shear_exp
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