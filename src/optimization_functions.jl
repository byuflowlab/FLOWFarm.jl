"""file with functions used in wind farm optimization, including constraints
created May 4, 2020
author: PJ Stanley and Jared Thomas
contributors: Nicholas F. Baker and Wesley Holt
"""


"""
    turbine_spacing(turbine_x,turbine_y)

Calculate the distance between turbines in a wind farm. There is an infinite gradient
of this function if two points are exactly the same. This can be avoided by returning the
square of the turbine spacing rather than the actual distance, but it makes the gradients
scale much more poorly. Because it is very vanishinly rare to have turbines exactly in the
same location, this function leaves the square root in the calculations.

# Arguments
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
"""
function turbine_spacing(turbine_x, turbine_y)
    nturbines = length(turbine_x)
    spacing_vec =
        zeros(typeof(turbine_x[1]), Int((nturbines) * (nturbines - 1) / 2))
    k = 1
    for i in 1:nturbines
        for j in i+1:nturbines
            spacing_vec[k] = sqrt(
                (turbine_x[j] - turbine_x[i])^2 +
                (turbine_y[j] - turbine_y[i])^2,
            )
            k += 1
        end
    end
    return spacing_vec
end

"""
    circle_boundary(center,radius,turbine_x,turbine_y)

calculate the distance from each turbine to a circular boundary. Negative means the
turbine is inside the boundary

# Arguments
- `center::Float`: circular boundary center [x,y]
- `radius::Float`: circulat boundary radius
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
"""
function circle_boundary(center, radius, turbine_x, turbine_y)
    nturbines = length(turbine_x)
    boundary_vec = zeros(typeof(turbine_x[1]), nturbines)
    for i in 1:nturbines
        boundary_vec[i] =
            sqrt((center[1] - turbine_x[i])^2 + (center[2] - turbine_y[i])^2) -
            radius
    end
    return boundary_vec
end


"""
    convex_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)

calculate the distance from each turbine to a possibly non-circular, but convex boundary. Negative means the
turbine is inside the boundary

# Arguments
- `boundary_vertices::Array{Float,2}`: vertices of the convex hull CCW in order s.t.
        boundaryVertices[i] -> first point of face for unit_normals[i]
- `boundary_normals::Array{Float,2}`: unit normal vector for each boundary face
        CCW where boundaryVertices[i] is the first point of the corresponding face
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
"""
function convex_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y)
    nturbines = length(turbine_x)
        nVertices = size(boundary_vertices)[1]
        # initialize array to hold distances from each point to each face
        face_distance = zeros(typeof(turbine_x[1]),(nturbines, nVertices))

        # loop through points and find distance to each face
        for i = 1:nturbines
                # determine if point is inside or outside of each face, and distance from each face
                for j = 1:nVertices
                        # define the vector from the point of interest to the first point of the face
                        pa = [boundary_vertices[j, 1]-turbine_x[i], boundary_vertices[j, 2]-turbine_y[i]]
                        # find perpendicular distance from point to current surface (vector projection)
                        d_vec = sum(pa .* boundary_normals[j,:]) .* boundary_normals[j,:]
                        # calculate the sign of perpendicular distance from point to current face (- is inside, + is outside)
                        face_distance[i, j] = -sum(d_vec .* boundary_normals[j,:])
                end
        end

        return vcat(face_distance...)
end


"""
    splined_boundary(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)

calculate the distance from each turbine to a closed boundary made up of zero or
more reflex angles (concavities). Boundary will have three or four user-selected
"corners", such that the "sides" between corners (that will be splined) are
injective functions (meaning that for every x-coord,there exists only one
corresponding y-coord). Returns four values for every turbine, corresponding to
the distance from the turb to the upper, lower, left, and right splined "sides".
A negative return value means the turb is outside the boundary for that "side".
Returns a single array of {Float64} of length {length(turbine_x) * 4}. Note that
all boundary coordinates must be in the first quadrant of the Cartesian
coordinate system (+x and +y values only)

# Arguments
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
- `bndry_x_clsd::Array{Float}`: x locations of boundary vertices, CCW in
        order s.t. boundaryVertices[0] is the NE corner of boundary, and with
        the first vertex duplicated at the end for completeness of calcs.
- `bndry_y_clsd::Array{Float}`: y locations of boundary vertices, CCW in
        order s.t. boundaryVertices[0] is the NE corner of boundary, and with
        the first vertex duplicated at the end for completeness of calcs.
- `bndry_corner_indcies::Array{Float}`: An array of 3 or 4 indicies in the
        bndry_x/y_clsd arrays that correspond to the three four "corners" used
        between splined "sides"
"""
function splined_boundary(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
    """ Returns if the passed turbines are within the passed closed boundary """
    num_turbs = Int8(length(turbine_x))
    num_sides = Int8(length(bndry_corner_indcies)-1)
    x_min_indx = 3                      # Default to work w/ squared boundaries
    if num_sides == 3                   # If we only have 3 corners
        x_min_indx = 2                  # denote that it's a triangle boundary
    end

    # Check to make sure our points are in
    bndry_cons = zeros(typeof(turbine_x[1]),(num_turbs, 4))   # 4 values (2 x and 2 y) for each turb

    x_max = bndry_x_clsd[bndry_corner_indcies[1]]           # Our maximum x-value
    x_min = bndry_x_clsd[bndry_corner_indcies[x_min_indx]]  # Our min x-value

    # For every turbine
    for cntr in 1:num_turbs
        #- Calc x-vals
        bndry_cons[cntr,1] = (x_max - turbine_x[cntr])   # Positive good, neg bad
        bndry_cons[cntr,2] = (turbine_x[cntr] - x_min)   # pos good, neg bad

        #- Calc y-vals
        y_max,y_min = getUpDwnYvals(turbine_x[cntr], bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
        bndry_cons[cntr,3] = (y_max - turbine_y[cntr])
        bndry_cons[cntr,4] = (turbine_y[cntr] - y_min)
    end

    # Debug code to print the constraint values
    # for i in 1:sum(num_turbs)
    #     @printf("%e %e %e %e\n", bndry_cons[i,1], bndry_cons[i,2], bndry_cons[i,3] ,bndry_cons[i,4])
    # end

    return collect(bndry_cons)
end

"""
    splined_boundary_discreet_regions(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies, turbs_per_region)

Uses FlowFarm's splined_boundary() function to calculate the turbine-boundary
constraints for one or more discreet regions, with pre-allocated turbines
for each region. Returns four values for every turbine, corresponding to the
distance from each turb to the upper, lower, left, and right splined "sides" for
the region to which it was allocated. A negative return value means the turb is
outside the "side" of boundary for which it has been allocated. Returns a single
array of {Float64} of length {length(turbine_x) * 4}. Note that all boundary
coordinates must be in the first quadrant of the Cartesian coordinate system
(+x and +y values only)

# Arguments
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
- `bndry_x_clsd::Array{Float}`: x locations of boundary vertices, CCW in
        order s.t. boundaryVertices[0] is the NE corner of boundary, and with
        the first vertex duplicated at the end for completeness of calcs.
- `bndry_y_clsd::Array{Float}`: y locations of boundary vertices, CCW in
        order s.t. boundaryVertices[0] is the NE corner of boundary, and with
        the first vertex duplicated at the end for completeness of calcs.
- `bndry_corner_indcies::Array{Float}`: An array of 3 or 4 indicies in the
        bndry_x/y_clsd arrays that correspond to the three four "corners" used
        between splined "sides"
- 'turbs_per_region::Array{Int}`: An array of length equivalent to the number of
        discrete boundary regions, with each element denoting howmany turbines
        are apportioned to the corresponding region. sum(turbs_per_region) must
        be equivalent to the total number of turbines in the windfarm
"""
function splined_boundary_discreet_regions(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies, turbs_per_region)
    """ Goes through numerous discrete splined-boundary regions and returns if the apportioned turbines are within their region """
    num_regions = length(turbs_per_region)
    bndry_constraints = [ Float64[] for i in 1:num_regions ]  # To hold cnstrnts

    #-- Loop through and do all regions --#
    prev_turb_index = 1
    for cntr in 1:num_regions
        next_turb_index = ((turbs_per_region[cntr]-1) + prev_turb_index)  # Next index for our Turbines
        region_turbine_x = turbine_x[prev_turb_index:next_turb_index]   # Simplified list of turbines preallocated to this region
        region_turbine_y = turbine_y[prev_turb_index:next_turb_index]
        bndry_constraints[cntr] = append!(bndry_constraints[cntr],splined_boundary(region_turbine_x, region_turbine_y, bndry_x_clsd[cntr], bndry_y_clsd[cntr], bndry_corner_indcies[cntr]))
        prev_turb_index = (turbs_per_region[cntr] + prev_turb_index)
    end

    # Make it a long 1D array for SNOPT
    bndry_constraints = collect(Iterators.Flatten(bndry_constraints))
    return bndry_constraints
end


"""
    getUpDwnYvals(turbine_x, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)

Supplements FlowFarm's splined_boundary() function by calculating (for a given x
location) the maximum and minimum y-value permitted to remain "inside" the
boundary. If turbine_x is located left of the boundary's leftmost vertex or
right of the boundary's rightmost vertex, it return's that corresponding
vertex's y-value as the max and min, as default. Returns two values, the minimum
and maximum interior y-values withing a boundary for the given turbine_x value.
Note that all boundary coordinates must be in the first quadrant of the
Cartesian coordinate system (+x and +y values only)

# Arguments
- `turbine_x::Array{Float}`: x-value of the turbine being examined
- `bndry_x_clsd::Array{Float}`: x locations of boundary vertices, CCW in
        order s.t. boundaryVertices[0] is the NE corner of boundary, and with
        the first vertex duplicated at the end for completeness of calcs.
- `bndry_y_clsd::Array{Float}`: y locations of boundary vertices, CCW in
        order s.t. boundaryVertices[0] is the NE corner of boundary, and with
        the first vertex duplicated at the end for completeness of calcs.
- `bndry_corner_indcies::Array{Float}`: An array of 3 or 4 indicies in the
        bndry_x/y_clsd arrays that correspond to the three four "corners" used
        between splined "sides"
"""
function getUpDwnYvals(turbine_x, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
    # Given that there are 4 "sides" (with 3&4 below, 1&2 above),
    # returns the splined y-vals the given x-coord falls between
    bool_triangle = false                # Default to work w/ squared boundaries
    if length(bndry_corner_indcies) == 4 # If we only have 3 corners (4 closed)
        bool_triangle = true             # denote that it's a triangle boundary
    end

    #-- Upper edge calculation
    if (bool_triangle) # If we're dealing with a 3-sided boundary
        # If it's out of bounds to the left of the upper boundary
        if (turbine_x <= bndry_x_clsd[bndry_corner_indcies[2]])
            # Give it the y-value of our leftmost point
            y_max = bndry_y_clsd[bndry_corner_indcies[2]]
        # If it's right of the leftmost corner, but left of the rightmost corner
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[1]])
            y_max = linear(bndry_x_clsd[bndry_corner_indcies[2]:-1:bndry_corner_indcies[1]],
                            bndry_y_clsd[bndry_corner_indcies[2]:-1:bndry_corner_indcies[1]], turbine_x)  # Make it the left upper spline #############
        else # Otherwise it's right of our top right corner
            # Give it the y-value of our rightmost corner
            y_max = bndry_y_clsd[bndry_corner_indcies[1]]
        end
    else # If we're dealing with a 4-sided boundary
        # If left of the bottom left corner
        if (turbine_x <= bndry_x_clsd[bndry_corner_indcies[3]])
            # Give it the y-value of our bottom left corner
            y_max = bndry_y_clsd[bndry_corner_indcies[3]]
        # If it's to the left of the right upper spline
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[2]])
            y_max = linear(bndry_x_clsd[bndry_corner_indcies[3]:-1:bndry_corner_indcies[2]],
                            bndry_y_clsd[bndry_corner_indcies[3]:-1:bndry_corner_indcies[2]], turbine_x)
        # If it's to the left of the rightmost point
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[1]])
            y_max = linear(bndry_x_clsd[bndry_corner_indcies[2]:-1:bndry_corner_indcies[1]],
                            bndry_y_clsd[bndry_corner_indcies[2]:-1:bndry_corner_indcies[1]], turbine_x)
        else # Otherwise, if it's out of bounds to the right
            # Give it the y-value of our rightmost point
            y_max = bndry_y_clsd[bndry_corner_indcies[1]]
        end
    end
    #-- Lower edge calculation.
    if (bool_triangle)
        # If it's to the left of the leftmost lower spline
        if (turbine_x <= bndry_x_clsd[bndry_corner_indcies[2]])
            y_min = bndry_y_clsd[bndry_corner_indcies[2]] # leftmost y-value is corner #2
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[3]])
            # Spline from bottom left to bottom middle
            y_min = linear(bndry_x_clsd[bndry_corner_indcies[2]:bndry_corner_indcies[3]],
                            bndry_y_clsd[bndry_corner_indcies[2]:bndry_corner_indcies[3]], turbine_x)
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[1]])
            # Spline from bottom middle to bottom right
            y_min = linear(bndry_x_clsd[bndry_corner_indcies[3]:bndry_corner_indcies[4]],
                            bndry_y_clsd[bndry_corner_indcies[3]:bndry_corner_indcies[4]], turbine_x)
        else
            # Give it the y-value of our rightmost point
            y_min = bndry_y_clsd[bndry_corner_indcies[1]]
        end
    else # if it's a square
        # If it's to the left of the leftmost lower spline
        if (turbine_x <= bndry_x_clsd[bndry_corner_indcies[3]])
            y_min = bndry_y_clsd[bndry_corner_indcies[3]] # leftmost y-value is corner #3
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[4]])
            # Spline from bottom left to bottom middle
            y_min = linear(bndry_x_clsd[bndry_corner_indcies[3]:bndry_corner_indcies[4]],
                            bndry_y_clsd[bndry_corner_indcies[3]:bndry_corner_indcies[4]], turbine_x)
        elseif (turbine_x < bndry_x_clsd[bndry_corner_indcies[5]])
            # Spline from bottom middle to bottom right
            y_min = linear(bndry_x_clsd[bndry_corner_indcies[4]:bndry_corner_indcies[5]],
                            bndry_y_clsd[bndry_corner_indcies[4]:bndry_corner_indcies[5]], turbine_x)
        else
            # Give it the y-value of our rightmost point
            y_min = bndry_y_clsd[bndry_corner_indcies[1]]
        end
    end

    return y_max, y_min
end

"""
    ray_trace_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)

Calculate the distance from each turbine to the nearest point on the boundary. 
Negative means the turbine is inside the boundary.

# Arguments
- `boundary_vertices::Array{Float,2}`: vertices of the convex hull CCW in order s.t.
        boundaryVertices[i] -> first point of face for unit_normals[i]
- `boundary_normals::Array{Float,2}`: unit normal vector for each boundary face
        CCW where boundaryVertices[i] is the first point of the corresponding face
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
"""
function ray_trace_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)

    # number of turbines and boundary vertices
    nturbines = length(turbine_x)
    nvertices = size(boundary_vertices)[1]

    # initialize constraint output values
    c = zeros(typeof(turbine_x[1]),(nturbines))

    # initialize array to hold distances from each turbine to closest boundary face
    turbine_to_face_distance = zeros(typeof(turbine_x[1]),(nvertices))

    # add the first boundary vertex again to the end of the boundary vertices vector (to form a closed loop)
    boundary_vertices = [boundary_vertices; boundary_vertices[1,1] boundary_vertices[1,2]]

    # iterate through each turbine location
    for i = 1:nturbines

        # initialize intersection counter
        intersection_counter = 0

        # get vector from turbine to the first vertex in first face
        turbine_to_first_facepoint = boundary_vertices[1, :] - [turbine_x[i]; turbine_y[i]]

        # iterate through each boundary
        for j = 1:nvertices

            # check if y-coordinate of turbine is less than at least one y-coordinate of the two boundary vertices
            if !(boundary_vertices[j, 2] < turbine_y[i] && boundary_vertices[j+1, 2] < turbine_y[i])        # (this might not be necessary)
                
                # check if x-coordinate of turbine is between the x-coordinates of the two boundary vertices
                if boundary_vertices[j, 1] < turbine_x[i] < boundary_vertices[j+1, 1] || boundary_vertices[j, 1] > turbine_x[i] > boundary_vertices[j+1, 1]

                    # check to see if the turbine is below the boundary
                    if turbine_y[i] < (boundary_vertices[j+1, 2] - boundary_vertices[j, 2]) / (boundary_vertices[j+1, 1] - boundary_vertices[j, 1]) * (turbine_x[i] - boundary_vertices[j, 1]) + boundary_vertices[j, 2]
                    
                        # the vertical ray intersects the boundary
                        intersection_counter += 1

                    end

                end

            end

            # define the vector from the turbine to the second point of the face
            turbine_to_second_facepoint = boundary_vertices[j+1, :] - [turbine_x[i]; turbine_y[i]]

            # find perpendicular distance from turbine to current face (vector projection)
            boundary_vector = boundary_vertices[j+1, :] - boundary_vertices[j, :]
            
            # check if perpendicular distance is the shortest
            if sum(boundary_vector .* -turbine_to_first_facepoint) > 0 && sum(boundary_vector .* turbine_to_second_facepoint) > 0
                
                # perpendicular distance from turbine to face
                turbine_to_face_distance[j] = abs(sum(turbine_to_first_facepoint .* boundary_normals[j,:]))
            
            # check if distance to first facepoint is shortest
            elseif sum(boundary_vector .* -turbine_to_first_facepoint) <= 0

                # distance from turbine to first facepoint
                turbine_to_face_distance[j] = sqrt(sum(turbine_to_first_facepoint.^2))

            # distance to second facepoint is shortest
            else

                # distance from turbine to second facepoint
                turbine_to_face_distance[j] = sqrt(sum(turbine_to_second_facepoint.^2))

            end
            
            # reset for next face iteration
            turbine_to_first_facepoint = turbine_to_second_facepoint        # (for efficiency, so we don't have to recalculate for the same vertex twice)
                
        end

        # magnitude of the constraint value
        c[i] = -ff.smooth_max(-turbine_to_face_distance, s=100.0)

        # sign of the constraint value (- is inside, + is outside)
        if mod(intersection_counter, 2) == 1
            c[i] = -c[i]
        end

    end

    return c

end