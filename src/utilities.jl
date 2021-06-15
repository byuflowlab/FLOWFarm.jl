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
        OVdYd = norm([wake_center_y-turbine_y, wake_center_z - turbine_z])
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

"""
    closeBndryLists(bndryPts_x, bndryPts_y)

Appends the 1st element to the end of each array for a closed boundary.
Note, this will not function properly if there is only one region.
For only one region, use `closeBndryList(bndryPts_x, bndryPts_y)` (note the
singular, not plural 'List' in the function title)

# Arguments
- `bndryPts_x::Array{Float,1}` : N-D array of x-coordinates for the vertices
        around N-many closed boundaries
- `bndryPts_y::Array{Float,1}` : N-D array of y-coordinates for the vertices
        around N-many closed boundaries
"""
function closeBndryLists(region_bndry_x, region_bndry_y)
    # Determine how many regions and points per region were passed
    nRegions = length(region_bndry_x)
    # For every region
    for i in 1:nRegions
        # Append the initial points to the end of that row (if needed)
        region_bndry_x[i], region_bndry_y[i] = closeBndryList(region_bndry_x[i], region_bndry_y[i])
    end
    
    return region_bndry_x, region_bndry_y
end

"""
    closeBndryList(bndryPts_x, bndryPts_y)

Appends the 1st element to the end of the coordinate arrays if it is not already
repeated. Note, this will only work on 1-D arrays. For an array of 1-D arrays,
use `closeBndryLists(bndryPts_x, bndryPts_y)` (note the
plural, not singular 'Lists' in the function title)

# Arguments
- `bndryPts_x::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndryPts_y::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
"""
function closeBndryList(bndry_x, bndry_y)
    # If the end point isn't already a repeat of the staring point
    if !((bndry_x[1] == bndry_x[end]) && (bndry_y[1] == bndry_y[end]))
        # Append the first point to the end of the array
        bndryPts_x_clsd = push!(bndry_x, bndry_x[1])
        bndryPts_y_clsd = push!(bndry_y, bndry_y[1])
    end

    return bndryPts_x_clsd, bndryPts_y_clsd
end

"""
    PointsOnCircum(center_x, center_y, r, n = 100)

Given a circle center, radius, and number of discrete points, returns an array
of discrete points along the circle's circumference

# Arguments
- `center_x::Float64` : cartesian x-coordinate for the center of the circle
- `center_y::Float64` : cartesian y-coordinate for the center of the circle
- `r::Float64` : distance from circle's center to the circumference points
- `n::Float64` : defaults to 100, is the number of discrete evenly-spaced points
        that will be returned along the circle's circumference
"""
function DiscreteCircum(center_x, center_y, r, n = 100)
    bndry_x = zeros(n+1)
    bndry_y = zeros(n+1)
    for i in 0:n
        bndry_x[i+1] = center_x + cos(2*pi/n*i)*r
        bndry_y[i+1] = center_y + sin(2*pi/n*i)*r
    end
    return bndry_x, bndry_y
end

"""
    calcMinorAngle(bndry_x, bndry_y, bndry_z=[0,0,0])

Given three points in space, calculates the magnitude of the non-reflex angle
formed at the center point. Created to be used in VR_bounary_startup()

# Arguments
- `bndry_x::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndry_y::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
- `bndry_z::Array{Float,1}` : 1-D array of z-coordinates for the vertices
        around a singlar closed boundary. Default to [0,0,0] for X-Y plane
"""
function calcMinorAngle(bndry_x, bndry_y, bndry_z=[0,0,0])
    # Calculates the magnitude of the non-reflex angle formed at the center point
    ABx = bndry_x[2]-bndry_x[1]
    ABy = bndry_y[2]-bndry_y[1]
    ABz = bndry_z[2]-bndry_z[1]
    
    BCx = bndry_x[2]-bndry_x[3]
    BCy = bndry_y[2]-bndry_y[3]
    BCz = bndry_z[2]-bndry_z[3]
    
    Num = (ABx*BCx) + (ABy*BCy) + (ABz*BCz) # Dot Product
    
    Denom = norm([ABx, ABy, ABz]) * norm([BCx, BCy, BCz]) # Multiplication of magnitudes
    Theta = acosd(Num/Denom) # Get the angle formed

    # If it's greater than 180, get the 
    if (Theta > 180.0)
        Theta = 360.0 - Theta    # Get the explementary angle
    end
    
    return Theta
end

"""
    calcSmallestAngle(bndry_x_clsd, bndry_y_clsd)

Given a 1-D closed array of boundary verticies (with first point repeated at the
end) it determines the smallest non-reflex angle created by any three
consecutive verticies along the boundary. Created to be used in
VR_bounary_startup()

# Arguments
- `bndry_x::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndry_y::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
"""
function calcSmallestAngle(bndry_x_clsd, bndry_y_clsd)
    num_angles = length(bndry_x_clsd)-1
    # Loop the second point back on so we get the angle at the beginning
    if !(((bndry_x_clsd[1] == bndry_x_clsd[end-1]) && (bndry_y_clsd[1] == bndry_y_clsd[end-1])) && ((bndry_x_clsd[2] == bndry_x_clsd[end]) && (bndry_y_clsd[2] == bndry_y_clsd[end])))
        bndryPts_x_loopd = vcat(bndry_x_clsd, bndry_x_clsd[2])
        bndryPts_y_loopd = vcat(bndry_y_clsd, bndry_y_clsd[2])
    end
    
    #- Calculate the smallest angle -#
    smallest_angle = 360
    for i in 1:num_angles
        temp_angle = calcMinorAngle(bndryPts_x_loopd[i:i+2], bndryPts_y_loopd[i:i+2])

        if (temp_angle < smallest_angle)
            smallest_angle = temp_angle
        end
    end
    
    return smallest_angle
end

"""
    getPerimeterLength(bndry_x_clsd, bndry_y_clsd)

Given a 1-D closed array of boundary verticies (with first point repeated at the
end) returns the length along the perimeter. Created to be used in
VR_bounary_startup()

# Arguments
- `bndry_x::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndry_y::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
"""
function getPerimeterLength(bndry_x_clsd, bndry_y_clsd)
    num_bndry_pts = length(bndry_x_clsd)-1
    nLength = zeros(num_bndry_pts)
    for i in 1:num_bndry_pts
        nLength[i] = norm([bndry_x_clsd[i+1]-bndry_x_clsd[i], bndry_y_clsd[i+1]-bndry_y_clsd[i]])
    end
    
    return nLength
end

"""
    coordDist(x1, y1, x2, y2)

Given a two points (x1, y1) and (x2, y2), returns the euclidean distance between
them

# Arguments
- `x1::Float64` : x-coord of the first point
- `y1::Float64` : y-coord of the first point
- `x2::Float64` : x-coord of the second point
- `y2::Float64` : y-coord of the second point
"""
function coordDist(x1, y1, x2, y2)
    xDiff = x1 - x2
    yDiff = y1 - y2
    return norm([xDiff, yDiff])
end

"""

    boundary_normals_calculator(boundary_vertices)

Outputs the unit vectors perpendicular to each boundary of a shape, given the Cartesian coordinates for the shape's vertices.

# Arguments
- `boundary_vertices::Array{Float,1}` : n-by-2 array containing all the boundary vertices, counterclockwise

"""

function boundary_normals_calculator(boundary_vertices)

    # get number of vertices in shape
    nvertices = length(boundary_vertices[:, 1])

    # add the first vertex to the end of the array to form a closed-loop
    boundary_vertices = [boundary_vertices; boundary_vertices[1,1] boundary_vertices[1,2]]

    # initialize array to hold boundary normals
    boundary_normals = zeros(nvertices, 2)

    # iterate over each boundary
    for i = 1:nvertices

        # create a vector normal to the boundary
        boundary_normals[i, :] = [ -(boundary_vertices[i+1, 2] - boundary_vertices[i, 2]) ; boundary_vertices[i+1, 1] - boundary_vertices[i, 1] ]
        
        # normalize the vector
        boundary_normals[i, :] = boundary_normals[i, :] / norm(boundary_normals[i, :])
    
    end

    return boundary_normals

end

function sunflower_points(n; alpha=0.0)
    # this function generates n points within a circle in a sunflower seed pattern
    # the code is based on the example found at
    # https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle

    function radius(k, n, b)
        if k > n - b
            r = 1 # put on the boundary
        else
            r = sqrt(k - 1.0/2.0)/sqrt(n - (b + 1.0)/2.0)  # apply squareroot
        end
        return r
    end

    x = zeros(n)
    y = zeros(n)

    b = round(alpha*sqrt(n)) # number of boundary points

    phi = (sqrt(5.0) + 1.0)/2.0  # golden ratio

    for k in 1:n
        r = radius(k, n, b)

        theta = 2.0*pi*k/phi^2

        x[k] = r*cos(theta)
        y[k] = r*sin(theta)
    end

    return x, y
end

"""
    rotor_sample_points(nsamplepoints=1)

Initializes the sampling locations in the rotor-swept-area. Returns values such that
zero is at the turbine hub location and 1 is at the tip of the blades. If a single
sample is requested, it will be at the hub location. Otherwise, the points will be located
using the sunflower packcing algorithm.

# Arguments
- `nsamplepoints::Int`: controlls how many sample points to generate
"""
function rotor_sample_points(nsamplepoints=1)

    if nsamplepoints > 1
        rotor_sample_points_y, rotor_sample_points_z = ff.sunflower_points(nsamplepoints)
    else
        rotor_sample_points_y = rotor_sample_points_z = [0.0]
    end

    return rotor_sample_points_y, rotor_sample_points_z
end

"""
    print_state_layouts_in_cartesian_frame(turbinex, turbiney, winddirections)

Given a wind farm layout in the global reference frame, print the layout rotated to the 
cartesian frame with wind to the positive x axis (right) for all wind directions.

# Arguments
- `turbinex::Array{T,1}`: x locations of turbines in global reference frame 
- `turbiney::Array{T,1}`: y locations of turbines in global reference frame
- `winddirections::Array{T,1}`: all wind directions in radians in meteorological coordinates (0 rad. = from North)
"""
function print_layout_in_cartesian_frame_excel(turbinex, turbiney, winddirections; center=[0.0,0.0], plot_layouts=false, round_locations=false)
    # get number of directions 
    ndirections = length(winddirections)

    # initialize excel file
    outfile = "round_38_turbines_cartesian.xlsx"
    XLSX.openxlsx(outfile, mode="w") do xf
    
        sheet = xf[1]
        sheet["A1"] = "wind direction"
        sheet["B1"] = "turbine x"
        sheet["C1"] = "turbine y"

        # loop through directions 
        for i in 1:ndirections

            # rotate turbinex and turbiney to current direction 
            rotx, roty = rotate_to_wind_direction(turbinex, turbiney, winddirections[i], center=center)

            # round if desired
            if round_locations
                rotx = round.(rotx)
                roty = round.(roty)
            end

            # plot layouts if desired, pause for 5 seconds on each
            if plot_layouts
                p = plot()
                rotordiameter = zeros(length(turbinex)) .+ 126.4
                ff.plotlayout!(p, rotx, roty, rotordiameter)
                ff.plotlayout!(p, [rotx[1], rotx[20]], [roty[1], roty[20]], rotordiameter, fillcolor=:red, title="$(winddirections[i]*180.0/pi)")
                display(p)

                sleep(5)

            end

            XLSX.rename!(sheet, "rotated layouts")

            # print rotated coordinates to specified output
            sheet["A$(i+1)"] = winddirections[i]*180.0/pi
            sheet["B$(i+1)"] = join(rotx, ",")
            sheet["C$(i+1)"] = join(roty, ",")

        end
    end
end

# the following is from floris for ccalcculating wake counts
# def number_of_wakes_iec(self, wd, return_turbines=True):
# """
# Finds the number of turbines waking each turbine for the given
# wind direction. Waked directions are determined using the formula
# in Figure A.1 in Annex A of the IEC 61400-12-1:2017 standard.
# # TODO: Add the IEC standard as a reference.
# Args:
#     wd (float): Wind direction for determining waked turbines.
#     return_turbines (bool, optional): Switch to return turbines.
#         Defaults to True.
# Returns:
#     list(int) or list( (:py:class:`~.turbine.Turbine`, int ) ):
#     Number of turbines waking each turbine and, optionally,
#     the list of Turbine objects in the map.
# TODO:
# - This could be reworked so that the return type is more consistent.
# - Describe the method used to find upstream turbines.
# """
# wake_list = []
# for coord0, turbine0 in self.items:
#     other_turbines = [
#         (coord, turbine) for coord, turbine in self.items if turbine != turbine0
#     ]
#     dists = np.array(
#         [
#             np.hypot(coord.x1 - coord0.x1, coord.x2 - coord0.x2)
#             / turbine.rotor_diameter
#             for coord, turbine in other_turbines
#         ]
#     )
#     angles = np.array(
#         [
#             np.degrees(np.arctan2(coord.x1 - coord0.x1, coord.x2 - coord0.x2))
#             for coord, turbine in self.items
#             if turbine != turbine0
#         ]
#     )
#     # angles = (-angles - 90) % 360
#     waked = dists <= 2.0
#     waked = waked | (
#         (dists <= 20.0)
#         & (
#             np.abs(wrap_180(wd - angles))
#             <= 0.5 * (1.3 * np.degrees(np.arctan(2.5 / dists + 0.15)) + 10)
#         )
#     )
#     if return_turbines:
#         wake_list.append((turbine0, waked.sum()))
#     else:
#         wake_list.append(waked.sum())
# return wake_list