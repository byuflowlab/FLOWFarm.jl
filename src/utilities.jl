
using DelimitedFiles 
"""
    latlong_to_xy(latitude, longitude, utm_zone; isnorth=true, units="m")

Converts arrays of points from latitude and longitude to x and y in meters in a local
coordinate frame based on the point with the lowest magnitude latitude,

# Arguments
- `latitude::Array{Float,N}`
- `longitude::Array{Float,N}`
- `isnorth::Float`: specifies if the point is in the northern hemisphere (defaul: true)
"""
function latlong_to_xy(latitude, longitude, utm_zone; isnorth=true)
    
    # get number of points
    npoints = length(latitude)

    # set up utm region
    utmregion = gd.UTMfromLLA(utm_zone, isnorth, wgs84)

    # get zero point index
    zp = argmin(latitude)

    # get zero point lat long
    minlla = gd.LLA(latitude[zp], longitude[zp])
    
    # get zero point utm
    minutm = utmregion(minlla)

    # get x and y for all points
    x = zeros(npoints)
    y = zeros(npoints)
    for i in 1:npoints
        # lat and long for current point
        lla = gd.LLA(latitude[i], longitude[i])

        # convert to utm coordinates
        utm = utmregion(lla)

        # get x location
        x[i] = utm.x - minutm.x

        # get y location
        y[i] = utm.y - minutm.y

    end
    
    return x, y
end

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
function hermite_spline(x, x0, x1, y0, dy0, y1, dy1; return_deriv=false)

    # initialize coefficients for parametric Hermite cubic spline
    c3 = (2.0*(y1))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) - 
        (2.0*(y0))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) + 
        (dy0)/(x0^2 - 2.0*x0*x1 + x1^2) + 
        (dy1)/(x0^2 - 2.0*x0*x1 + x1^2)

    c2 = (3.0*(y0)*(x0 + x1))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) - 
        ((dy1)*(2.0*x0 + x1))/(x0^2 - 2.0*x0*x1 + x1^2) - ((dy0)*(x0 +
        2.0*x1))/(x0^2 - 2.0*x0*x1 + x1^2) - (3.0*(y1)*(x0 + x1))/(x0^3 -
        3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3)

    c1 = ((dy0)*(x1^2 + 2.0*x0*x1))/(x0^2 - 2.0*x0*x1 + x1^2) + ((dy1)*(x0^2 +
        2.0*x1*x0))/(x0^2 - 2.0*x0*x1 + x1^2) - (6.0*x0*x1*(y0))/(x0^3 -
        3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3) + (6.0*x0*x1*(y1))/(x0^3 -
        3.0*x0^2*x1 + 3.0*x0*x1^2 - x1^3)

    c0 = ((y0)*(- x1^3 + 3.0*x0*x1^2))/(x0^3 - 3.0*x0^2*x1 + 3.0*x0*x1^2 -
        x1^3) - ((y1)*(- x0^3 + 3.0*x1*x0^2))/(x0^3 - 3.0*x0^2*x1 +
        3.0*x0*x1^2 - x1^3) - (x0*x1^2*(dy0))/(x0^2 - 2.0*x0*x1 + x1^2) - 
        (x0^2*x1*(dy1))/(x0^2 - 2.0*x0*x1 + x1^2)

    # Solve for y and dy values at the given point
    y = c3*x^3 + c2*x^2 + c1*x + c0
    dy_dx = c3*3*x^2 + c2*2*x + c1

    if return_deriv
        return y, dy_dx
    else
        return y
    end
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
    smooth_max(x; s=10.0)

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
function smooth_max(x; s=10.0)

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

function return_model_set(filename)

    # import FlowFarm; const ff = FlowFarm

    # set initial turbine x and y locations
    turbine_x = [-3.0, 0.0, 3.0, 0.0, 0.0, -1.5, 0.0, 1.5, 0.0].*80.0
    turbine_y = [0.0, 3.0, 0.0, -3.0, 0.0, 0.0, 1.5, 0.0, -1.5].*80.0

    # calculate the number of turbines
    nturbines = length(turbine_x)

    # set turbine base heights
    turbine_z = zeros(nturbines)

    # set turbine yaw values
    turbine_yaw = zeros(nturbines)

    # set turbine design parameters
    rotor_diameter = zeros(nturbines) .+ 80.0 # m
    hub_height = zeros(nturbines) .+ 70.0   # m
    cut_in_speed = zeros(nturbines) .+4.  # m/s
    cut_out_speed = zeros(nturbines) .+25.  # m/s
    rated_speed = zeros(nturbines) .+16.  # m/s
    rated_power = zeros(nturbines) .+2.0E6  # W
    generator_efficiency = zeros(nturbines) .+0.944

    # rotor swept area sample points (normalized by rotor radius)
    rotor_points_y = [0.0]
    rotor_points_z = [0.0]

    # set flow parameters
    wind_speed = 8.0
    air_density = 1.1716  # kg/m^3
    ambient_ti = 0.077
    shearexponent = 0.15
    winddirections = [275.0*pi/180.0, 0.0, pi]
    windspeeds = [wind_speed, wind_speed, wind_speed]
    windprobabilities = [1.0/3.0,1.0/3.0,1.0/3.0]
    ambient_tis = [ambient_ti, ambient_ti, ambient_ti]
    measurementheight = [hub_height[1], hub_height[1], hub_height[1]]

    # load power curve
    powerdata = readdlm("inputfiles/niayifar_vestas_v80_power_curve_observed.txt",  ',', skipstart=1)
    pvelpoints = powerdata[:,1]
    powerpoints = powerdata[:,2]*1E6

    # initialize power model
    power_model = PowerModelPowerPoints(pvelpoints, powerpoints)
    power_models = Vector{typeof(power_model)}(undef, nturbines)
    for i = 1:nturbines
        power_models[i] = power_model
    end

    # load thrust curve
    ctdata = readdlm("inputfiles/predicted_ct_vestas_v80_niayifar2016.txt",  ',', skipstart=1)
    ctvelpoints = ctdata[:,1]
    ctpoints = ctdata[:,2]

    # initialize thurst model
    ct_model1 = ThrustModelCtPoints(ctvelpoints, ctpoints)
    ct_models = Vector{typeof(ct_model1)}(undef, nturbines)
    for i = 1:nturbines
        ct_models[i] = ct_model1
    end

    # initialize wind shear model
    wind_shear_model = PowerLawWindShear(shearexponent)

    # get sorted indecies 
    sorted_turbine_index = sortperm(turbine_x)

    # initialize the wind resource definition
    windresource = DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

    # set up wake and related models
    wakedeficitmodel = GaussYaw()
    wakedeflectionmodel = GaussYawDeflection()
    wakecombinationmodel = LinearLocalVelocitySuperposition()
    localtimodel = LocalTIModelMaxTI()

    # initialize model set
    model_set = WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    return model_set, windresource, power_models, ct_models, turbine_x, turbine_y, turbine_z, turbine_yaw, 
        rotor_diameter, hub_height, cut_in_speed, cut_out_speed, rated_speed, rated_power,
        generator_efficiency, rotor_points_y, rotor_points_z
    
end