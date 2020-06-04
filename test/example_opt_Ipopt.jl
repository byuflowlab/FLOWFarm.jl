using JuMP, Ipopt
using DelimitedFiles 
using PyPlot
import ForwardDiff


#---FUNCTIONS FOR CONSTRAINTS---

# set up boundary constraint wrapper function at specific index
# it takes in x as separate elements instead of as a vector (needed for Ipopt)
function boundary_wrapper_atindex(ind, x...)

    # combine the inputs into a single vector
    x = collect(x)

    # include globals for optimization with Ipopt
    global xlast
    global boundary_constraint_values

    # check if the constraint calculation is necessary
    if !isequal(x, xlast)
        
        global boundary_constraint_gradients
        global spacing_constraint_values
        global spacing_constraint_gradients

        # perform all constraint calculations for this set of design variables
        boundary_constraint_values, boundary_constraint_gradients, spacing_constraint_values, spacing_constraint_gradients = con(x)

        # reset xlast
        xlast = x

    end

    # return the desired value
    return boundary_constraint_values[Int(ind)]

end

# set up boundary constraint gradient wrapper function at specific index 
# it takes in x as separate elements instead of as a vector (needed for Ipopt)
function boundary_gradient_wrapper_atindex(g, ind, x...)

    # combine the inputs into a single vector
    x = collect(x)

    # include globals for optimization with Ipopt
    global xlast
    global boundary_constraint_gradients

    # check if the constraint calculation is necessary
    if !isequal(x, xlast)
        
        global boundary_constraint_values
        global spacing_constraint_values
        global spacing_constraint_gradients

        # perform all constraint calculations for this set of design variables
        boundary_constraint_values, boundary_constraint_gradients, spacing_constraint_values, spacing_constraint_gradients = con(x)

        # reset xlast
        xlast = x

    end

    # fill in the provided gradient vector
    g[1] = 0.0 # the first spot is for the index argument, which isn't an actual design variable, so we set the partial with respect to it as zero
    for i = 2:length(g)
        g[i] = boundary_constraint_gradients[Int(ind), i-1]
    end

    # return the gradient vector
    return g

end

# set up spacing constraint wrapper function at specific index
# it takes in x as separate elements instead of as a vector (needed for Ipopt)
function spacing_wrapper_atindex(ind, x...)

    # combine the inputs into a single vector
    x = collect(x)

    # include globals for optimization with Ipopt
    global xlast
    global spacing_constraint_values

    # check if the constraint calculation is necessary
    if !isequal(x, xlast)
        
        global boundary_constraint_values
        global boundary_constraint_gradients
        global spacing_constraint_gradients

        # perform all constraint calculations for this set of design variables
        boundary_constraint_values, boundary_constraint_gradients, spacing_constraint_values, spacing_constraint_gradients = con(x)

        # reset xlast
        xlast = x

    end

    # return the desired value
    return spacing_constraint_values[Int(ind)]

end

# set up spacing constraint gradient wrapper function at specific index 
# it takes in x as separate elements instead of as a vector (needed for Ipopt)
function spacing_gradient_wrapper_atindex(g, ind, x...)

    # combine the inputs into a single vector
    x = collect(x)

    # include globals for optimization with Ipopt
    global xlast
    global spacing_constraint_gradients

    # check if the constraint calculation is necessary
    if !isequal(x, xlast)
        
        global boundary_constraint_values
        global boundary_constraint_gradients
        global spacing_constraint_values

        # perform all constraint calculations for this set of design variables
        boundary_constraint_values, boundary_constraint_gradients, spacing_constraint_values, spacing_constraint_gradients = con(x)

        # reset xlast
        xlast = x

    end

    # fill in the provided gradient vector
    g[1] = 0.0 # the first spot is for the index argument, which isn't an actual design variable, so we set the partial with respect to it as zero
    for i = 2:length(g)
        g[i] = spacing_constraint_gradients[Int(ind), i-1]
    end

    # return the gradient vector
    return g

end

# set up constraint function wrapper
# when this function is called, all constaints values and gradients are calculated
function con(x)

    # get boundary constraint values and gradients
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)
    
    # get spacing constraint values and gradients
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    return boundary_con, db_dx, spacing_con, ds_dx
    
end

# set up boundary constraint wrapper function
# this function actually calculates the constraint value
function boundary_wrapper(x)

    # include relevant globals
    global boundary_vertices
    global boundary_normals

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get boundary constraint values
    return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y)

end

# set up spacing constraint wrapper function
# this function actually calculates the constraint value
function spacing_wrapper(x)
    
    # include relevant globals
    global rotor_diameter

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]
    
    # get spacing constraint values
    return 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)

end

#---FUNCTIONS FOR OBJECTIVE---

# set up objective wrapper function
# it takes in x as separate elements instead of as a vector (needed for Ipopt)
function aep_wrapper_atindex(x...)

    # combine inputs into a single vector
    x = collect(x)

    # get the objective 
    AEP = aep_wrapper(x)[1]
    
    # return the objective
    return AEP
end

# set up objective gradient wrapper function
# it takes in x as separate elements instead of as a vector (needed for Ipopt)
function aep_gradient_wrapper(g, x...)

    x = collect(x)
    gradient = ForwardDiff.jacobian(aep_wrapper, x)
    for i = 1:size(g)[1]
        g[i] = gradient[i]
    end

    return g
end

# set up objective wrapper function
# this functino actually calculates the objective value
function aep_wrapper(x)

    # include relevant globals
    global turbine_z
    global rotor_diameter
    global hub_height
    global turbine_yaw
    global ct_model
    global generator_efficiency
    global cut_in_speed
    global cut_out_speed
    global rated_speed
    global rated_power
    global windresource
    global power_model
    global model_set
    global rotor_points_y
    global rotor_points_z
    global obj_scale    

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_model, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # return the objective as a vector
    return [AEP]
end


# import model set with wind farm and related details
include("./model_sets/model_set_6.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_vertices = ([0 0; 1 0; 1 .75; .75 .75; .75 1; 0 1] .- .5).*500 # Utah-shape boundary
boundary_normals = [0 1.0; -1 0; 0 -1; -1 0; 0 -1; 1 0]

# get the number of vertices of the boundary
n_vertices = length(boundary_vertices)

# initialize variables for constraint wrapper functions
xlast = 0
boundary_constraint_values = 0
boundary_constraint_gradients = 0
spacing_constraint_values = 0
spacing_constraint_gradients = 0
global xlast
global boundary_constraint_values
global boundary_constraint_gradients
global spacing_constraint_values
global spacing_constraint_gradients

# set globals for use in wrapper functions
global model_set
global rotor_points_y
global rotor_points_z
global turbine_z
global ambient_ti
global rotor_diameter
global boundary_vertices
global boundary_normals
global obj_scale

# initialize design variable array
x_initial = [copy(turbine_x);copy(turbine_y)]
n_designvariables = length(x_initial)

# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# get number of spacing constraints
function numberofspacingconstraints(x)
    n = 0
    for i = 1:x; n += i; end
    return n
end
n_spacingconstraints = numberofspacingconstraints(length(turbine_x)-1)

# set general lower and upper bounds for design variables
# lb = minimum(boundary_vertices)
# ub = maximum(boundary_vertices)

# initialize the model for optimization
model = Model(Ipopt.Optimizer)

# create the design variables
@variable(model, x[i = 1:n_designvariables], start = x_initial[i])

# register the required functions
register(model, :boundary, n_designvariables + 1, boundary_wrapper_atindex, boundary_gradient_wrapper_atindex)
register(model, :spacing, n_designvariables + 1, spacing_wrapper_atindex, spacing_gradient_wrapper_atindex)
register(model, :AEP, n_designvariables, aep_wrapper_atindex, aep_gradient_wrapper)

# set the objective
@NLobjective(model, Max, AEP(x...))

# set the constraints
@NLconstraint(model, boundary_con[i = 1:nturbines], boundary(i, x...) <= 0)
@NLconstraint(model, spacing_con[i = 1:n_spacingconstraints], spacing(i, x...) <= 0)

# run and time optimization
t1 = time()
optimize!(model)
t2 = time()
clkt = t2-t1
xopt = value.(x)
fopt = objective_value(model)
info = (termination_status(model), primal_status(model), dual_status(model))

# print optimization results
println()
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("start objective value: ", aep_wrapper(x_initial)[1])
println("end objective value: ", aep_wrapper(xopt)[1])

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# add wind farm boundary to plot
plt.gcf().gca().plot([boundary_vertices[:,1];boundary_vertices[1,1]],[boundary_vertices[:,2];boundary_vertices[1,2]], color="C2")

# set up and show plot
axis("square")
xlim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
ylim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
plt.show()
