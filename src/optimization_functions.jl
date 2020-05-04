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
        spacing_vec = zeros(Int((nturbines)*(nturbines-1)/2))
        k = 1
        for i = 1:nturbines
                for j = i+1:nturbines
                        spacing_vec[k] = sqrt((turbine_x[j]-turbine_x[i])^2+(turbine_y[j]-turbine_y[i])^2)
                        k += 1
                end
        end
        return spacing_vec
end


# function farm_boundary(turbine_x,turbine_y)
# end
