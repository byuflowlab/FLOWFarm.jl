"""file with functions used in wind farm optimization, including constraints
created May 4, 2020
author: PJ Stanley and Jared Thomas
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
function turbine_spacing(turbine_x,turbine_y)
        nturbines = length(turbine_x)
        spacing_vec = zeros(typeof(turbine_x[1]),Int((nturbines)*(nturbines-1)/2))
        k = 1
        for i = 1:nturbines
                for j = i+1:nturbines
                        spacing_vec[k] = sqrt((turbine_x[j]-turbine_x[i])^2+(turbine_y[j]-turbine_y[i])^2)
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
function circle_boundary(center,radius,turbine_x,turbine_y)
        nturbines = length(turbine_x)
        boundary_vec = zeros(typeof(turbine_x[1]),nturbines)
        for i = 1:nturbines
                boundary_vec[i] = sqrt((center[1]-turbine_x[i])^2 + (center[2]-turbine_y[i])^2) - radius
        end
        return boundary_vec
end


"""
    windfarm_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)

calculate the distance from each turbine to a circular boundary. Negative means the
turbine is inside the boundary

# Arguments
- `boundary_vertices::Array{Float,2}`: vertices of the convex hull CCW in order s.t.
        boundaryVertices[i] -> first point of face for unit_normals[i]
- `boundary_normals::Array{Float,2}`: unit normal vector for each boundary face
        CCW where boundaryVertices[i] is the first point of the corresponding face
- `turbine_x::Array{Float}`: turbine x locations
- `turbine_y::Array{Float}`: turbine y locations
"""
function windfarm_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)
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
        c[i] = -ff.smooth_max_ndim(-turbine_to_face_distance)

        # sign of the constraint value (- is inside, + is outside)
        if mod(intersection_counter, 2) == 1
            c[i] = -c[i]
        end

    end

    return c

end


"""
    smooth_max_ndim(x; s=10.0)

Calculate the smoothmax (a.k.a. softmax or LogSumExponential) of x and y.

Based on John D. Cook's writings at 
(1) https://www.johndcook.com/blog/2010/01/13/soft-maximum/
and
(2) https://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/

And based on article in FeedlyBlog
(3) https://blog.feedly.com/tricks-of-the-trade-logsumexp/

# Arguments
- `x::Array{Float64,1}`` - vector with all the input values
- `s::Float64` - controls the level of smoothing used in the smooth max
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