"""
    hermite_spline(x, x0, x1, y0, dy0, y1, dy1)
    
Produces the y and (optionally) dy values for a hermite cubic spline
interpolating between two end points with known slopes

# Arguments
- `x::Float`: x position of output y
- `x0::Float`: x position of upwind endpoint of spline
- `x1::Float`: x position of downwind endpoint of spline
- `y0::Float`: y position of upwind endpoint of spline
- `dy0::Float`: slope at upwind endpoint of spline
- `y1::Float`: y position of downwind endpoint of spline
- `dy1::Float`: slope at downwind endpoint of spline
"""
function hermite_spline(x, x0, x1, y0, dy0, y1, dy1, return_deriv=false)

    # initialize coefficients for parametric Hermite cubic spline
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
    dy_dx = c3*3*x^2 + c2*2*x + c1

    if return_deriv
        return y, dy_dx
    else
        return y
end

"""
    overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y,
    wake_center_z, wake_diameter; tol=1E-6)
    
Produces the y and (optionally) dy values for a hermite cubic spline
interpolating between two end points with known slopes

All dimensions should be in meters

# Arguments
- `turbine_y::Float`: cross wind location of turbine hub
- `turbine_z::Float`: vertical location of turbine hub
- `rotor_diameter::Float`
- `wake_center_y::Float`: cross wind location of wake center
- `wake_center_z::Float`: vertical location of wake center
- `wake_diameter::Float`
- `tol::Float`: extra distance for comparisons. Default 1E-6
"""
# calculates the overlap area between a given wake and a rotor area
function overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y,
    wake_center_z, wake_diameter; tol=1E-6)

    # distance between wake center and rotor center
    if (wake_center_z > (turbine_z + tol)) || (wake_center_z < (turbine_z - tol))
        OVdYd = sqrt((wake_center_y-turbine_y)^2 + (wake_center_z - turbine_z)^2)
    elseif (wake_center_y > (turbine_y + tol))
        OVdYd = wake_center_y - turbine_y
    elseif (turbine_y > (wake_center_y + tol))
        OVdYd = turbine_y - wake_center_y
    else
        OVdYd = 0.0
    end

    # find rotor radius
    OVr = rotor_diameter/2.0

    # find wake radius
    OVRR = wake_diameter/2.0

    # make sure the distance from wake center to turbine hub is positive
    OVdYd = abs(OVdYd)

    # determine if there is overlap
    if (OVdYd < (OVr+OVRR)) # if the rotor overlaps the wake zone

        # check that turbine and wake centers are not perfectly aligned
        if (OVdYd > (0.0 + tol))

            # check if the rotor is wholly contained in the wake
            if ((OVdYd + OVr) < OVRR + tol)
                # wake_overlap = pi*OVr*OVr
                wake_overlap = 3.1415926535897*OVr*OVr
            elseif ((OVdYd + OVRR) < OVr + tol)
                # wake_overlap = pi*OVRR*OVRR
                wake_overlap = 3.1415926535897*OVRR*OVRR
            else
                # calculate the distance from the wake center to the chord connecting the lens cusps
                OVL = (-OVr*OVr+OVRR*OVRR+OVdYd*OVdYd)/(2.0*OVdYd)

                OVz = sqrt(OVRR*OVRR-OVL*OVL)
                OVz2 = sqrt(OVr*OVr-(OVdYd-OVL)*(OVdYd-OVL))

                wake_overlap = OVRR*OVRR*acos(OVL/OVRR) + OVr*OVr*acos((OVdYd-OVL)/OVr) - OVL*OVz - (OVdYd-OVL)*OVz2
            end

        # perfect overlap case where the wake is larger than the rotor
        elseif (OVRR > OVr)
            wake_overlap = pi*OVr*OVr
        # perfect overlap case where the rotor is larger than the wake
        else
            wake_overlap = pi*OVRR*OVRR
        end

    # case with no overlap
    else
        wake_overlap = 0.0
    end

end

"""
    smooth_max_ndim(x; s=100.0)

Calculate the smoothmax (a.k.a. softmax or LogSumExponential) of the elements in x.

Based on John D. Cook's writings at 
(1) https://www.johndcook.com/blog/2010/01/13/soft-maximum/
and
(2) https://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/

# Arguments
- `x::Float`: first value for comparison
- `y::Float`: second value for comparison
- `s::Float` : controls the level of smoothing used in the smooth max
"""
function smooth_max(x, y; s=10.0)

    # LogSumExponential Method - used this in the past
    # g = (x*exp(s*x)+y*exp(s*y))/(exp(s*x)+exp(s*y))

    # non-overflowing version of Smooth Max function (see ref 2 above)
    max_val = max(x, y)
    min_val = min(x, y)
    r = (log(1.0 + exp(s*(min_val - max_val))) + s*max_val)/s

    return r

end


"""
    smooth_max_ndim(x; s=100.0)

Calculate the smoothmax (a.k.a. softmax or LogSumExponential) of the elements in x.

Based on John D. Cook's writings at 
(1) https://www.johndcook.com/blog/2010/01/13/soft-maximum/
and
(2) https://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/

And based on article in FeedlyBlog
(3) https://blog.feedly.com/tricks-of-the-trade-logsumexp/

# Arguments
- `x::Array{Float,1}` : vector with all the input values
- `s::Float` : controls the level of smoothing used in the smooth max
"""
function smooth_max_ndim(x; s=100.0)

    # non-overflowing version of Smooth Max function (see ref 2 and 3 above)
    
    # get the maximum value and the index of maximum value
    max_val, max_ind = findmax(x)

    # get the indices of x
    indices = collect(1:length(x))

    # remove the index of the maximum value
    splice!(indices, max_ind)

    # LogSumExp with smoothing factor s
    r = (log(sum([1.0; exp.(s*(x[indices].-max_val))])) + s*max_val)/s

    return r

end