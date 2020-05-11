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
        face_distance = zeros(nturbines, nVertices)

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
        return face_distance
end
