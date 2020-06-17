using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_center
    params.boundary_radius

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(boundary_center, boundary_radius, turbine_x, turbine_y)
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # include relevant globals
    params.rotor_diameter

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return spacing distances
    return 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
function aep_wrapper(x, params)
    # include relevant globals
    params.turbine_z
    params.rotor_diameter
    params.hub_height
    params.turbine_yaw
    params.ct_model
    params.generator_efficiency
    params.cut_in_speed
    params.cut_out_speed
    params.rated_speed
    params.rated_power
    params.windresource
    params.power_models
    params.model_set
    params.rotor_points_y
    params.rotor_points_z
    params.obj_scale

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
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
include("./model_sets/model_set_9_38turb_round_farm.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_center = [0.0,0.0]
boundary_radius = 1225.8227848101264

# set globals for use in wrapper functions
struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_center
    boundary_radius
    obj_scale
    hub_height
    turbine_yaw
    ct_model
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    windresource
    power_models
end

params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_center, boundary_radius, obj_scale, hub_height, turbine_yaw, 
    ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]
global x

# report initial objective value
println("starting objective value: ", aep_wrapper(x, params)[1])

plot(0,0)
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
options["Major optimality tolerance"] = 1e-6
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary.out"
options["Print file"] = "print.out"

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)
obj_func(x) = wind_farm_opt(x)

# set up for WEC optimization
wec_steps = 6
wec_max = 3.0
wec_end = 1.0
wec_values = collect(LinRange(wec_max, wec_end, wec_steps))
println(wec_values)
info = fill("",wec_steps)
# run and time optimization
t1t = time()
for i in 1:length(wec_values)
    global x
    global xopt
    println("Running with WEC = ", wec_values[i])
    params.model_set.wake_deficit_model.wec_factor[1] = wec_values[i]
    
    t1 = time()
    xopt, fopt, info[i] = snopt(obj_func, x, lb, ub, options)
    t2 = time()
    clk = t2-t1
    # print optimization results
    println("Finished in : ", clk, " (s)")
    println("info: ", info)
    println("end objective value: ", -fopt)
    println("locations ", x[1:5])
    println("locations opt ", xopt[1:5])
    x = deepcopy(xopt)
end
t2t = time()
clkt = t2t - t1t
# print optimization results
println("Overall Stats: ")
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
