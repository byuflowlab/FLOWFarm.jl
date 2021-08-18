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
scale much more poorly. Because it is very vanishingly rare to have turbines exactly in the
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
            spacing_vec[k] = norm([turbine_x[j] - turbine_x[i], turbine_y[j] - turbine_y[i]])
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
        boundary_vec[i] = (center[1] - turbine_x[i])^2 + (center[2] - turbine_y[i])^2 - radius^2
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
injective functions (meaning that for every x-coord, there exists only one
corresponding y-coord). Returns four values for every turbine, corresponding to
the distance from the turb to the upper, lower, left, and right splined "sides".
A negative return value means the turb is inside the boundary for that "side".
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
    num_cons = Int(num_turbs * 4)
    bndry_cons = zeros(typeof(turbine_x[1]), num_cons)   # 4 values (2 x and 2 y) for each turb

    x_max = bndry_x_clsd[bndry_corner_indcies[1]]           # Our maximum x-value
    x_min = bndry_x_clsd[bndry_corner_indcies[x_min_indx]]  # Our min x-value

    # For every turbine
    for cntr in 1:num_turbs
        place = (cntr-1)*4
        #- Calc x-vals
        bndry_cons[place+1] = (x_max - turbine_x[cntr])   # Positive good, neg bad
        bndry_cons[place+2] = (turbine_x[cntr] - x_min)   # pos good, neg bad

        #- Calc y-vals
        y_max,y_min = getUpDwnYvals(turbine_x[cntr], bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
        bndry_cons[place+3] = (y_max - turbine_y[cntr])
        bndry_cons[place+4] = (turbine_y[cntr] - y_min)
    end

    return -bndry_cons # Invert values so Negative is inside boundary
end

"""
    splined_boundary_discreet_regions(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies, turbs_per_region)

Uses FLOWFarm's splined_boundary() function to calculate the turbine-boundary
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
function splined_boundary_discreet_regions(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, num_bndry_verts, bndry_corner_indcies, turbs_per_region)
    """ Goes through numerous discrete splined-boundary regions and returns if the apportioned turbines are within their region """
    num_regions = length(turbs_per_region)
    bndry_constraints = zeros(typeof(turbine_x[1]), sum(turbs_per_region)*4)#[ Float64[] for i in 1:num_regions ]  # To hold cnstrnts
    #-- Loop through and do all regions --#
    prev_turb_index = 1
    bndry_vert_index = 1
    for cntr in 1:num_regions
        next_turb_index = ((turbs_per_region[cntr]-1) + prev_turb_index)  # Next index for our Turbines
        region_turbine_x = turbine_x[prev_turb_index:next_turb_index]   # Simplified list of turbines preallocated to this region
        region_turbine_y = turbine_y[prev_turb_index:next_turb_index]
        cnstrnts_index_strt = ((prev_turb_index-1)*4)+1
        cnstrnts_index_end = ((next_turb_index-1)*4)+4
        bndry_constraints[cnstrnts_index_strt:cnstrnts_index_end] = splined_boundary(region_turbine_x, region_turbine_y, bndry_x_clsd[cntr], bndry_y_clsd[cntr], bndry_corner_indcies[bndry_vert_index:(bndry_vert_index+(num_bndry_verts[cntr]-1))])
        bndry_vert_index += num_bndry_verts[cntr]
        prev_turb_index += turbs_per_region[cntr]
    end

    # Make it a long 1D array for SNOPT
    return bndry_constraints
end


"""
    getUpDwnYvals(turbine_x, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)

Supplements FLOWFarm's splined_boundary() function by calculating (for a given x
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
    ray_casting_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)

Calculate the distance from each turbine to the nearest point on the boundary using 
the ray-casting algorithm. Negative means the turbine is inside the boundary.

# Arguments
- `boundary_vertices::Array{Float,2}`: vertices of the boundary CCW in order s.t.
        boundaryVertices[i] -> first point of face for unit_normals[i]
- `boundary_normals::Array{Float,2}`: unit normal vector for each boundary face
        CCW where boundaryVertices[i] is the first point of the corresponding face
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
"""
function ray_casting_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y; discrete=false, s=100)
    # discrete=boundary.discrete

    # single region
    if discrete == false

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
                    turbine_to_face_distance[j] = norm(turbine_to_first_facepoint)

                # distance to second facepoint is shortest
                else

                    # distance from turbine to second facepoint
                    turbine_to_face_distance[j] = norm(turbine_to_second_facepoint)

                end
                
                # reset for next face iteration
                turbine_to_first_facepoint = turbine_to_second_facepoint        # (for efficiency, so we don't have to recalculate for the same vertex twice)
                    
            end

            # magnitude of the constraint value
            c[i] = -ff.smooth_max(-turbine_to_face_distance, s=s)

            # sign of the constraint value (- is inside, + is outside)
            if mod(intersection_counter, 2) == 1
                c[i] = -c[i]
            end

        end

        return c

    # multiple discrete regions
    else

        # number of turbines
        nturbines = length(turbine_x)

        # number of regions
        nregions = length(boundary_vertices)

        # initialize constraint output values
        c = zeros(typeof(turbine_x[1]),(nturbines))

        # initialize turbine status vector
        status = zeros(Int64, nturbines)

        # iterate through each turbine location
        for i = 1:nturbines

            # iterate through each region
            for k = 1:nregions

                # number of boundary vertices
                nvertices = size(boundary_vertices[k])[1]

                # add the first boundary vertex again to the end of the boundary vertices vector (to form a closed loop)
                boundary_vertices_closed = [boundary_vertices[k]; boundary_vertices[k][1,1] boundary_vertices[k][1,2]]

                # initialize intersection counter
                intersection_counter = 0

                # get vector from turbine to the first vertex in first face
                turbine_to_first_facepoint = boundary_vertices_closed[1, :] - [turbine_x[i]; turbine_y[i]]

                # iterate through each boundary
                for j = 1:nvertices

                    # check if y-coordinate of turbine is less than at least one y-coordinate of the two boundary vertices
                    if !(boundary_vertices_closed[j, 2] < turbine_y[i] && boundary_vertices_closed[j+1, 2] < turbine_y[i])        # (this might not be necessary)
                        
                        # check if x-coordinate of turbine is between the x-coordinates of the two boundary vertices
                        if boundary_vertices_closed[j, 1] < turbine_x[i] < boundary_vertices_closed[j+1, 1] || boundary_vertices_closed[j, 1] > turbine_x[i] > boundary_vertices_closed[j+1, 1]

                            # check to see if the turbine is below the boundary
                            if turbine_y[i] < (boundary_vertices_closed[j+1, 2] - boundary_vertices_closed[j, 2]) / (boundary_vertices_closed[j+1, 1] - boundary_vertices_closed[j, 1]) * (turbine_x[i] - boundary_vertices_closed[j, 1]) + boundary_vertices_closed[j, 2]
                            
                                # the vertical ray intersects the boundary
                                intersection_counter += 1

                            end

                        end

                    end

                end

                # check if the turbine is in the current region
                if mod(intersection_counter, 2) == 1

                    # initialize array to hold distances from each turbine to closest boundary face
                    turbine_to_face_distance = zeros(typeof(turbine_x[1]),(nvertices))

                    # iterate through each boundary
                    for j = 1:nvertices

                        # define the vector from the turbine to the second point of the face
                        turbine_to_second_facepoint = boundary_vertices_closed[j+1, :] - [turbine_x[i]; turbine_y[i]]

                        # find perpendicular distance from turbine to current face (vector projection)
                        boundary_vector = boundary_vertices_closed[j+1, :] - boundary_vertices_closed[j, :]
                        
                        # check if perpendicular distance is the shortest
                        if sum(boundary_vector .* -turbine_to_first_facepoint) > 0 && sum(boundary_vector .* turbine_to_second_facepoint) > 0
                                
                            # perpendicular distance from turbine to face
                            turbine_to_face_distance[j] = abs(sum(turbine_to_first_facepoint .* boundary_normals[k][j,:]))
                        
                        # check if distance to first facepoint is shortest
                        elseif sum(boundary_vector .* -turbine_to_first_facepoint) <= 0

                            # distance from turbine to first facepoint
                            turbine_to_face_distance[j] = norm(turbine_to_first_facepoint)

                        # distance to second facepoint is shortest
                        else

                            # distance from turbine to second facepoint
                            turbine_to_face_distance[j] = norm(turbine_to_second_facepoint)

                        end
                        
                        # reset for next face iteration
                        turbine_to_first_facepoint = turbine_to_second_facepoint        # (for efficiency, so we don't have to recalculate for the same vertex twice)

                    end

                    # magnitude of the constraint value
                    c[i] = ff.smooth_max(-turbine_to_face_distance, s=100.0)
                    status[i] = 1

                end

            end

            # check if the turbine is in none of the regions
            if status[i] == 0

                # initialize array to hold distances from each turbine to closest boundary face
                turbine_to_face_distance = zeros(typeof(turbine_x[1]),0)

                # iterate through each region
                for k = 1:nregions

                    # number of boundary vertices
                    nvertices = size(boundary_vertices[k])[1]
        
                    # add the first boundary vertex again to the end of the boundary vertices vector (to form a closed loop)
                    boundary_vertices_closed = [boundary_vertices[k]; boundary_vertices[k][1,1] boundary_vertices[k][1,2]]

                    # get vector from turbine to the first vertex in first face
                    turbine_to_first_facepoint = boundary_vertices_closed[1, :] - [turbine_x[i]; turbine_y[i]]

                    # iterate through each boundary
                    for j = 1:nvertices

                        # define the vector from the turbine to the second point of the face
                        turbine_to_second_facepoint = boundary_vertices_closed[j+1, :] - [turbine_x[i]; turbine_y[i]]

                        # find perpendicular distance from turbine to current face (vector projection)
                        boundary_vector = boundary_vertices_closed[j+1, :] - boundary_vertices_closed[j, :]
                        
                        # check if perpendicular distance is the shortest
                        if sum(boundary_vector .* -turbine_to_first_facepoint) > 0 && sum(boundary_vector .* turbine_to_second_facepoint) > 0
                            
                            # perpendicular distance from turbine to face
                            push!(turbine_to_face_distance, abs(sum(turbine_to_first_facepoint .* boundary_normals[k][j,:])))
                        
                        # check if distance to first facepoint is shortest
                        elseif sum(boundary_vector .* -turbine_to_first_facepoint) <= 0

                            # distance from turbine to first facepoint
                            push!(turbine_to_face_distance, norm(turbine_to_first_facepoint))

                        # distance to second facepoint is shortest
                        else

                            # distance from turbine to second facepoint
                            push!(turbine_to_face_distance, norm(turbine_to_second_facepoint))

                        end
                        
                        # reset for next face iteration
                        turbine_to_first_facepoint = turbine_to_second_facepoint        # (for efficiency, so we don't have to recalculate for the same vertex twice)

                    end

                end

                # magnitude of the constraint value
                c[i] = -ff.smooth_max(-turbine_to_face_distance, s=s)
                status[i] = 1

            end

        end

        return c

    end

end

"""
    VR_bounary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)

Determines if the requested number of turbines can be placed along the closed
boundary with spacing and corner constraints. If the requested <num_turbs> is
too many, places as many turbines as possible along the boundary, and returns
the number of turbines not placed. NOTE: A shortcoming is that the
smallest-angled corner limits the spacing of all turbines. in the worst case,
a very thin boundary area would prevent any more than one turbine being placed
on the boundary, though more would be optimal. Future work would check to make
sure this corner (and the length of its adjacent sides) don't actually require
limiting the minimum distance between turbines.

# Arguments
- `bndry_x::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndry_y::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
- `start_dist::Float64`: the distance (positive or negative) along the boundary
        from the first boundary point where the turbines will begin to be placed
- `turb_min_spacing::Float64`: the fixed distance along the boundary's edge between
        adjacent turbines
- 'num_turbs::Float64`: the number of turbines to be placed around the boundary.
        Note that this function assumes VR_bounary_startup() has already been
        run so that the user won't attempt to place too many turbines.
- 'bndry_seg_length::Array{Int}`: an array of the lengths between adjacent
        boundary verticies, corresponding to how they appear in bndry_x and _y
"""
function VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
    #- Get arc-length -#
    bndry_seg_length = getPerimeterLength(bndry_x_clsd,bndry_y_clsd)
    bndry_tot_length = sum(bndry_seg_length)

    #- Determine how much distance between turbs for worst case (shallow corner) -#
    limiting_angle = calcSmallestAngle(bndry_x_clsd, bndry_y_clsd)  # The smallest angle in our boundary
    turb_corner_spacing = turb_min_spacing/sind(limiting_angle/2.0) # Min spacing neeed for the smallest angle
    if (turb_corner_spacing > turb_min_spacing)   # If the corners are tighter than our minimum distance
        turb_min_spacing = turb_corner_spacing    # That dictates our new minimum distance
    end

    #- Determine how many turbines we can actually place with corner constraint -#
    max_num_turbs = floor(Int, bndry_tot_length/turb_min_spacing) # Get max number of turbines that can fit the perimeter
    num_leftover_turbs = 0
    if (max_num_turbs < num_turbs)      # If the max is less than the number allocated
        num_leftover_turbs = floor(Int, num_turbs - max_num_turbs) # Note how many turbines we aren't placing
        num_turbs = floor(Int, max_num_turbs)       # Place as many as will fit
    end

    #- Evenly space the turbines we'll place -#
    turb_min_spacing = bndry_tot_length / num_turbs
    #- Make sure our start point is within the boundary length bounds
    start_dist = mod(start_dist, bndry_tot_length)
    #- Place the correct number of turbs along the boundary <start_dist> away from the first vertex -#
    turbine_x, turbine_y = VR_boundary(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)

    # Return the x- and y- coordinates of every turbine, and how many weren't placed
    return turbine_x,turbine_y, num_leftover_turbs
end

"""
    VR_boundary(bndry_x_clsd, bndry_y_clsd, start_dist, turb_spacing, num_turbs, bndry_seg_length)

Uses the Boundary portion of Boundary-Grid variable reduction method
place turbines along a closed wind farm boundary and perturb their location with
one (1) variable <start_dist>.  NOTE: Use of this function assumes prior use of
VR_bounary_startup(), which ensures the number of turbines placed on the
boundary doesn't violate any minimum spacing rules eiter along the boundary or
around corners.

# Arguments
- `bndry_x::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndry_y::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
- `start_dist::Float64`: the distance (positive or negative) along the boundary
        from the first boundary point where the turbines will begin to be placed
- `turb_spacing::Float64`: the fixed distance along the boundary's edge between
        adjacent turbines
- 'num_turbs::Float64`: the number of turbines to be placed around the boundary.
        Note that this function assumes VR_bounary_startup() has already been
        run so that the user won't attempt to place too many turbines.
"""
function VR_boundary(bndry_x_clsd, bndry_y_clsd, start_dist, turb_spacing, num_turbs)
    # Initialize necessary variables
    bndry_seg_length = getPerimeterLength(bndry_x_clsd,bndry_y_clsd)
    num_segs = length(bndry_x_clsd)-1
    # Initialize turbine locations
    turbine_x = zeros(num_turbs)
    turbine_y = zeros(num_turbs)
    
    #- Figure out where the starting point should be -#
    # "leg" is distance until next placed turbine
    # "seg" is distance between boundary verticies
    leg_remaining = start_dist
    curr_seg = 1
    for i in 1:num_segs   # Looping through the segments till we get there
        curr_seg = i
        if (bndry_seg_length[i] < leg_remaining)  # If this segment length is less than the start length
            leg_remaining -= bndry_seg_length[i]  # Clip the start length and move to the next one
        else                                        # Otherwise the start point is on this segment
            percent_to_start = leg_remaining / bndry_seg_length[i] # How far along our segment we are
            # Translate how far along the segment we are to actual x- y-coords
            turbine_x[1] = bndry_x_clsd[i] + ((bndry_x_clsd[i+1] - bndry_x_clsd[i]) * percent_to_start)
            turbine_y[1] = bndry_y_clsd[i] + ((bndry_y_clsd[i+1] - bndry_y_clsd[i]) * percent_to_start)
            break
        end
    end
    # Get how much distance is left on this leg after the starting point
    percent_left_of_segment = 1 - abs((turbine_x[1] - bndry_x_clsd[curr_seg]) / (bndry_x_clsd[curr_seg+1] - bndry_x_clsd[curr_seg]))
    seg_remaining = bndry_seg_length[curr_seg] * percent_left_of_segment
    leg_remaining = turb_spacing
    
    #- Place the rest of the turbines -#
    for i in 2:num_turbs    # For every turbine we have to place
        if(seg_remaining > leg_remaining)   #- If there's space on this leg to place the next turbine
            percent_to_place = (bndry_seg_length[curr_seg] - seg_remaining + turb_spacing) / bndry_seg_length[curr_seg]
            seg_remaining -= turb_spacing # Take out the distance we used
        else                                #- If there isn't enough space on this boundary segment to place the next turbine
            while(seg_remaining < leg_remaining)           # Keep checking boundary segments until we reach the next placement
                leg_remaining -= seg_remaining             # Subtract that much till our next placement
                curr_seg = mod(curr_seg, (num_segs))+1   # Increment what segment we're onspot
                seg_remaining = bndry_seg_length[curr_seg] # Reset how much segment length we have left
            end
            percent_to_place = leg_remaining / bndry_seg_length[curr_seg]
            seg_remaining = bndry_seg_length[curr_seg] - leg_remaining
        end
        #- Place the turbines the appropriate distance from the segment start point -#
        turbine_x[i] = bndry_x_clsd[curr_seg] + ((bndry_x_clsd[curr_seg+1] - bndry_x_clsd[curr_seg]) * percent_to_place)
        turbine_y[i] = bndry_y_clsd[curr_seg] + ((bndry_y_clsd[curr_seg+1] - bndry_y_clsd[curr_seg]) * percent_to_place)
        leg_remaining = turb_spacing  # Reset how much length till the next turbine is placed
    end
    
    return turbine_x, turbine_y
end

"""
    iea37cs4BndryVRIntPM(bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies, turbine_x, turbine_y, turb_diam, turb_min_space, num_turbs_to_place)

Uses the Variable reduction method for placing boundary turbines, and the
Partition Method (from splined_boundary()) for random interior points,
maintaining proper spacing from all previously placed turbines.

# Arguments
- `bndry_x_clsd::Array{Float,1}` : 1-D array of x-coordinates for the vertices
        around a singlar closed boundary
- `bndry_y_clsd::Array{Float,1}` : 1-D array of y-coordinates for the vertices
        around a singlar closed boundary
- `bndry_corner_indcies::Float64`: The indicies within <bndry_x_clsd> and
        <bndry_y_clsd> which denote the "corners" adjacent turbines
- 'turb_min_space::Float64`: For proximity knowledge, the minimum spacing
        required between any two turbines
- 'num_bndry_turbs::Float64`: The number of turbines desired to be placed along
        the boundary. If too many are selected (due to spacing condtraints), the
        remaining will be placed in the interior
- 'num_tot_turbs::Float64`: The number of total turbines to be placed both on
        the boundary and in the interior
"""
function iea37cs4BndryVRIntPM(bndry_x_clsd, bndry_y_clsd, bndry_corner_indicies, turb_min_space, num_bndry_turbs, num_tot_turbs)
    #-- Place all the boundary turbines we can --#
    bndry_tot_len = sum(getPerimeterLength(bndry_x_clsd,bndry_y_clsd))
    #- Make a random starting point along the boundary -#
    start_dist = rand(Float64) * bndry_tot_len
    #- Place the boundary turbines -# 
    turbine_x_bndry, turbine_y_bndry, num_leftover_turbs = VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_space, num_bndry_turbs)
    #- Determine how many will be placed in the interior -#
    num_bndry_turbs = num_bndry_turbs - num_leftover_turbs
    num_interior_turbs = num_tot_turbs - num_bndry_turbs
    
    # Initialize full list of turbine locations
    turbine_x = zeros(num_tot_turbs)
    turbine_y = zeros(num_tot_turbs)
    #- Fill in the ones ew've already placed along the boundary
    turbine_x[1:num_bndry_turbs] = turbine_x_bndry
    turbine_y[1:num_bndry_turbs] = turbine_y_bndry
    
    #-- Initialize interior space --#w
    num_sides = length(bndry_corner_indicies)-1
    #- Get the x-values -#
    x_min_indx = 3                      # Default to work w/ squared boundaries
    if num_sides == 3                   # If we only have 3 corners
        x_min_indx = 2                  # denote that it's a triangle boundary
    end
    x_max = bndry_x_clsd[bndry_corner_indicies[1]]           # Our maximum x-value
    x_min = bndry_x_clsd[bndry_corner_indicies[x_min_indx]]  # Our min x-value
    turbine_x[num_bndry_turbs+1:end] = (x_max - x_min)*rand(Float64, num_interior_turbs) .+ x_min # Get random x-values for the interior turbines

    #-- Get the y-values --#
    #- Determine the upper and lower splines to use for the given x -#
    # Fake for-loop here for proximity checking.
    i = num_bndry_turbs + 1   # Start after our last boundary turbine
    while (i <= num_tot_turbs)
        # Get the max and min y-value for this x-value
        y_min, y_max = getUpDwnYvals(turbine_x[i], bndry_x_clsd, bndry_y_clsd, bndry_corner_indicies)
        # Generate a random y-value within these limits
        turbine_y[i] = (y_max - y_min)*rand(Float64) + y_min # Get a random number in our bounds
        #- Check it doesn't conflict with aleady placed turbines -#
        for j in 1:(i-1) # Check the new ones we've place so far
            # If this turbine has a proximity conflict
            if (coordDist(turbine_x[i], turbine_y[i], turbine_x[j], turbine_y[j]) < turb_min_space)
                turbine_x[i] = (x_max - x_min)*rand(Float64) + x_min # Give it a new x-val
                i = i-1 # Redo the y-val too
                break # Stop checking for conflicts and redo the y-value
            end
        end
        i += 1
    end

    return turbine_x, turbine_y, num_interior_turbs
end