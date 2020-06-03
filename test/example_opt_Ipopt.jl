using JuMP, Ipopt
using DelimitedFiles 
using PyPlot
import ForwardDiff


#---FUNCTIONS FOR CONSTRAINTS---

# set up boundary constraint wrapper function at specific index (needed for Ipopt)
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
    return boundary_constraint_values[ind]

end

# set up boundary constraint gradient wrapper function at specific index (needed for Ipopt)
function boundary_gradient_wrapper_atindex(x...)

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

    # return the desired value
    return boundary_constraint_gradients[ind, :]

end

# set up spacing constraint wrapper function at specific index (needed for Ipopt)
function spacing_wrapper_atindex(x...)

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
    return spacing_constraint_values[ind]

end

# set up spacing constraint gradient wrapper function at specific index (needed for Ipopt)
function spacing_gradient_wrapper_atindex(x...)

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

    # return the desired value
    return spacing_constraint_gradients[ind, :]

end

# set up boundary constraint wrapper function
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

# set up constraint function wrapper
function con(x)

    # get boundary constraint values and gradients
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)
    
    # get spacing constraint values and gradients
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    return boundary_con, db_dx, spacing_con, ds_dx
    
end


#---FUNCTIONS FOR OBJECTIVE---

# set up objective wrapper function
function aep_wrapper(x...)

    # put the inputs into a vector (if not already a vector)
    if isa(x, Tuple)
        x = collect(x)
    end

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
    
    # return the objective as an array
    return [AEP]
end

# set up objective gradient wrapper function
function aep_gradient_wrapper(x...)

    x = collect(x)
    dAEP_dx = ForwardDiff.jacobian(aep_wrapper,x)

end


# set up optimization problem wrapper function
function wind_farm_opt(x)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # combine constraint values and jacobians into overall constraint value and jacobian arrays
    c = [spacing_con; boundary_con]
    dcdx = [ds_dx; db_dx]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end


# import model set with wind farm and related details
include("./model_sets/model_set_6.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_vertices = [0 0; 1 0; 1 .75; .75 .75; .75 1; 0 1].*400
boundary_normals = [0 1.0; -1 0; 0 -1; -1 0; 0 -1; 1 0]
n_vertices = length(boundary_vertices)

# initialize variables for wrapper functions
xlast = 0
boundary_constraint_values = 0
boundary_constraint_gradients = 9
spacing_constraint_values = 0
spacing_constraint_gradients = 0

# set globals for use in wrapper functions
global model_set
global rotor_points_y
global rotor_points_z
global turbine_z
global ambient_ti
global rotor_diameter
global boundary_center
global boundary_radius
global obj_scale

global xlast
global boundary_constraint_values
global boundary_constraint_gradients
global spacing_constraint_values
global spacing_constraint_gradients

# initialize design variable array
x_initial = [copy(turbine_x);copy(turbine_y)]
n_designvariables = length(x_initial)

# initialize the model for optimization
model = Model(Ipopt.Optimizer)

# create the design variables
@variable(model, x[i = 1:n_designvariables], start = x_initial[i])

# register the required functions
register(model, :boundary, n_designvariables, boundary_wrapper, boundary_gradient_wrapper)
register(model, :spacing, n_designvariables, spacing_wrapper, spacing_gradient_wrapper)
register(model, :AEP, n_designvariables, aep_wrapper, aep_gradient_wrapper)

# set the objective
@NLobjective(model, Max, AEP(x...))

# set auxiliary variables for constraints
# @variable(model, )

# set the constraints
# @NLexpression(model, my_expr[i = 1:n_vertices], boundary(x...)[i])
@NLconstraint(model, boundary_con[i = 1:3], boundary(x...) <= 0)
@NLconstraint(model, spacing_con, spacing(x...) .<= 0)

# report initial objective value
println("starting objective value: ", aep_wrapper(x)[1])

# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# set general lower and upper bounds for design variables
lb = zeros(length(x)) .- boundary_radius
ub = zeros(length(x)) .+ boundary_radius

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 3
options["Major optimality tolerance"] = 1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary.out"
options["Print file"] = "print.out"

# run and time optimization
t1 = time()
xopt, fopt, info = snopt(wind_farm_opt, x, lb, ub, options)
t2 = time()
clkt = t2-t2

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end objective value: ", aep_wrapper(xopt))

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# add wind farm boundary to plot
plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))

# set up and show plot
axis("square")
xlim(-boundary_radius-200,boundary_radius+200)
ylim(-boundary_radius-200,boundary_radius+200)
plt.show()
