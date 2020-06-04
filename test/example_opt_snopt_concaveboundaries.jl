using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff

# set up boundary constraint wrapper function
function boundary_wrapper(x)
    # include relevant globals
    # global boundary_center
    # global boundary_radius
    global boundary_vertices
    global boundary_normals

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.ray_trace_boundary(boundary_vertices,boundary_normals,turbine_x,turbine_y)
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

    # get and return spacing distances
    return 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
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
    
    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt(x)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
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
# boundary_center = [0.0,0.0]
# boundary_radius = 300.0
boundary_vertices = ([0 0; 1 0; 1 .75; .75 .75; .75 1; 0 1] .- .5).*500 # Utah-shape boundary
boundary_normals = [0 1.0; -1 0; 0 -1; -1 0; 0 -1; 1 0]

# set globals for use in wrapper functions
global model_set
global rotor_points_y
global rotor_points_z
global turbine_z
global ambient_ti
global rotor_diameter
# global boundary_center
# global boundary_radius
global boundary_vertices
global boundary_normals
global obj_scale

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]

# report initial objective value
println("starting objective value: ", aep_wrapper(x)[1])

# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# set general lower and upper bounds for design variables\
lb = zeros(length(x)) .+ minimum(boundary_vertices) # create a square that completely encloses the wind farm boundary
ub = zeros(length(x)) .+ maximum(boundary_vertices)

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 1
options["Major optimality tolerance"] = 1e-12
options["Major iteration limit"] = 1e8
options["Summary file"] = "snopt_summary.out"
options["Print file"] = "snopt_print.out"

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
plt.gcf().gca().plot([boundary_vertices[:,1];boundary_vertices[1,1]],[boundary_vertices[:,2];boundary_vertices[1,2]], color="C2")

# set up and show plot
axis("square")
xlim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
ylim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
plt.show()
