using SNOW
using Snopt
using DelimitedFiles 
using PyPlot
# using Distrib@uted

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_center
    params.boundary_radius

    # get number of turbines
    # println("x in bw: ", size(x))
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
    # println("x in sw: ", size(x))
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
    turbine_z = params.turbine_z
    rotor_diameter = params.rotor_diameter
    hub_height = params.hub_height
    turbine_yaw =params.turbine_yaw
    ct_models = params.ct_models
    generator_efficiency = params.generator_efficiency
    cut_in_speed = params.cut_in_speed
    cut_out_speed = params.cut_out_speed
    rated_speed = params.rated_speed
    rated_power = params.rated_power
    windresource = params.windresource
    power_models = params.power_models
    model_set = params.model_set
    rotor_points_y = params.rotor_points_y
    rotor_points_z = params.rotor_points_z
    obj_scale = params.obj_scale
    # println("x in aw: ", size(x))
    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # return the objective as an array
    return AEP
end

# set up optimization problem wrapper function
function wind_farm_opt!(g, x)
    # println("x in wfo: ", size(x))
    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    # ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x)
    # db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[1:(end-nturbines)] = spacing_con[:]
    g[end-nturbines+1:end] = boundary_con[:]
    # println(g)
    # quit()
    # dcdx = [ds_dx; db_dx]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    # dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)

    # set fail flag to false
    # fail = false
    # println("constraints: ", minimum(g[length(g)-nturbines:end]), " ", maximum(g[length(g)-nturbines:end]))
    # println("obj: ", AEP)
    # return objective, constraint, and jacobian values
    # quit()
    return AEP #, dAEP_dx, dcdx, fail
end

# import model set with wind farm and related details
include("./model_sets/model_set_9_38turb_round_farm.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-0#1E-11
xyscale = 1 #1E4

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
    xyscale
    hub_height
    turbine_yaw
    ct_models
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    windresource
    power_models
end

params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_center, boundary_radius, obj_scale, xyscale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x0 = [copy(turbine_x);copy(turbine_y)]
xinit = deepcopy(x0)
println(size(x0))
# report initial objective value
aep_init = aep_wrapper(x0, params)[1]
println("starting objective value: ", aep_init)

plot(0,0)
# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)
obj_func!(g, x) = wind_farm_opt!(g, x)

# set SNOPT options
snopt_opt = Dict(
    "Derivative option" => 1,
    "Verify level" => 0,
    "Major optimality tolerance" => 1e-4,
    "Major iterations limit" => 1E5,
    "Summary file" => "snopt_summary1.out",
    "Print file" => "snopt_print1.out"
)
solver = SNOPT(options=snopt_opt)
# ip_options = Dict(
#     "max_iter" => 3,
#     "tol" => 1e-6
# )
# solver = IPOPT(ip_options)
options = Options(;solver, derivatives=ForwardAD())

# set general lower and upper bounds for design variables
lx = zeros(length(x0)) .- boundary_radius
ux = zeros(length(x0)) .+ boundary_radius

# set general lower and upper bounds for constraints
ng = Int(nturbines + (nturbines)*(nturbines - 1)/2)
# println("ng: ", ng)
lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines)]
println(minimum(lg), maximum(lg))
# quit()
ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines)]
# println("ug: ", minimum(lg), maximum(ug))
# run and time optimization
# quit()
t1 = time()

# println(length(x0))
# println(length(lg))
xopt, fopt, info, out = minimize(obj_func!, x0, ng, lx, ux, lg, ug, options)
t2 = time()
clk = t2-t1

println("xopt ", xopt)
aep_final = aep_wrapper(xopt, params)

# print optimization results
println("Finished in : ", clk, " (s)")
println("info: ", info)
println("end objective value: ", -fopt)
# println("major iter = ", out.major_iter)
# println("iterations = ", out.iterations)
# println("solve time = ", out.run_time)
println("AEP improvement (%) = ", 100*(aep_final - aep_init)/aep_init) 
println("opt locs: ", xopt)
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
