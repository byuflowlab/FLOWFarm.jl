using Ipopt
using DelimitedFiles 
using PyPlot
import ForwardDiff

using CSV
using DataFrames
# using YAML

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices
    params.boundary_normals
    
    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y)
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
    params.power_model
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
                cut_out_speed, rated_speed, rated_power, windresource, power_model, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # return the objective as an array
    return [AEP]
end

# objective function
function obj(x) 
    global funcalls_AEP
    AEP = -aep_wrapper(x)[1]
    append!(funcalls_AEP, -AEP)
    return AEP
end

# constraint function
function con(x, g)
    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)

    # calculate boundary constraint
    boundary_con = boundary_wrapper(x)

    # combine spacing and boundary constraint values
    g[:] = [spacing_con; boundary_con]
end

# objective gradient function
function obj_grad(x, grad_f)
    grad_f[:] = -ForwardDiff.jacobian(aep_wrapper,x)
end

# constraint gradients function
function con_grad(x, mode, rows, cols, values)
    if mode == :Structure
        # report the sparcity structure of the jacobian
        for i = 1:prob.m
            rows[(i-1)*prob.n+1:i*prob.n] = i*ones(Int,prob.n)
            cols[(i-1)*prob.n+1:i*prob.n] = 1:prob.n
        end
        return rows, cols
    else
        # calculate spacing constraint jacobian
        ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

        # calculate boundary constraint jacobian
        db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

        # combine constraint jacobians into an overall array of constraint derivatives
        for i = 1:prob.m
            for j = 1:prob.n
                values[(i-1)*prob.n+j] = [ds_dx; db_dx][i, j]
            end
        end
    end
    return values
end

# this function will be called at each iteration of the optimization
function intermediate(
    alg_mod::Int,
    iter_count::Int,
    obj_value::Float64,
    inf_pr::Float64, inf_du::Float64,
    mu::Float64, d_norm::Float64,
    regularization_size::Float64,
    alpha_du::Float64, alpha_pr::Float64,
    ls_trials::Int)
    global iter_AEP
    append!(iter_AEP, -obj_value)
    
    return true
end

# set globals for iteration history
iter_AEP = zeros(Float64, 0)
funcalls_AEP = zeros(Float64, 0)
global iter_AEP
global funcalls_AEP

# import model set with wind farm and related details
include("./model_sets/model_set_7_ieacs3.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-12

# set wind farm boundary parameters
boundary_vertices = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5;
    8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
boundary_normals = [0.9829601758936983 -0.1838186405319916; 0.9934614633172962 -0.11416795042154541; 0.9987121579438882 -0.050734855622757584; 
    0.9998686751666075 -0.01620593781838486; 0.9999954987444023 0.0030004151269687495; -0.9998078216567232 -0.019604074934516894; -0.6957179389375846 -0.718315076718037; 
    -0.6957275377423737 -0.7183057797532565; -0.8019887481131871 0.5973391397688945; 0.5138086803485797 0.8579047965820281; 0.4252760929807897 0.905063668886888; 
    0.2645057513093967 0.9643841078762402; -0.0684295708121141 0.9976559496331737; -0.39636379138742883 0.9180935381958544; -0.6828023205475376 0.7306031693435896; 
    -0.7996740386176392 0.6004343694034798; -0.8578802011411015 0.5138497450520954; 0.42552559023380465 0.9049463918134445]

# set globals for use in wrapper functions
struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_vertices
    boundary_normals
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
    power_model
end

params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, obj_scale, hub_height, turbine_yaw, 
    ct_model, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_model)

# initialize design variable array
x_initial = [copy(turbine_x);copy(turbine_y)]

# calculate initial state AEPs
state_aeps = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_model, model_set;
                rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.0*24.0)

dir_aep = zeros(20)
for i in 1:20
    for j in 1:20
        dir_aep[i] += state_aeps[(i-1)*20 + j]
    end
end

# add initial turbine location to plot
plt.clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# get number of design variables
n_designvariables = length(x_initial)

# get number of constraints
function numberofspacingconstraints(nturb)
    # calculates number of spacing constraints needed for given number of turbines
    ncon = 0
    for i = 1:nturb-1; ncon += i; end
    return ncon
end
n_spacingconstraints = numberofspacingconstraints(nturbines)
n_boundaryconstraints = length(boundary_wrapper(x_initial, params))
n_constraints = n_spacingconstraints + n_boundaryconstraints

# set general lower and upper bounds for design variables
lb = ones(n_designvariables) * -Inf
ub = ones(n_designvariables) * Inf

# set lower and upper bounds for constraints
lb_g = ones(n_constraints) * -Inf
ub_g = zeros(n_constraints)

# create the problem
prob = createProblem(n_designvariables, lb, ub, n_constraints, lb_g, ub_g, n_designvariables*n_constraints, 0,
    obj, con, obj_grad, con_grad)
addOption(prob, "hessian_approximation", "limited-memory")
prob.x = x_initial

# set function to be ran at every iteration
setIntermediateCallback(prob, intermediate)

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)

# run and time optimization
t1 = time()
status = solveProblem(prob)
t2 = time()
clkt = t2-t1
xopt = prob.x
fopt = prob.obj_val
info = Ipopt.ApplicationReturnStatus[status]

# print initial objective value
println("nturbines: ", nturbines)
println("rotor diameter: ", rotor_diameter[1])
println("starting AEP value (MWh): ", aep_wrapper(x_initial, params)[1]*1E5)

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end AEP value (MWh): ", aep_wrapper(xopt)[1]*1E5)

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# add wind farm boundary to plot
plt.gcf().gca().plot([boundary_vertices[:,1];boundary_vertices[1,1]],[boundary_vertices[:,2];boundary_vertices[1,2]], color="C2")

# set up and show turbine layout plot
axis("equal")
xlim(5500, 11000)
ylim(-500, 7000)
plt.show()
savefig("ipoptresults_turbinelayout.png")

# create obj value history plots
plt.clf()
plt.gcf().gca().plot(1:length(funcalls_AEP), funcalls_AEP)
title("Objective Value History")
xlabel("Function Call")
ylabel("AEP (TWh)")
savefig("ipoptresults_plot_AEP_history_by_funcall.png")

plt.clf()
plt.gcf().gca().plot(1:length(iter_AEP), iter_AEP)
title("Objective Value at Each Iteration")
xlabel("Iteration")
ylabel("AEP (TWh)")
savefig("ipoptresults_plot_AEP_history_by_iter.png")

# write text file with optimization summary info
io = open("ipoptresults_summary.txt", "a")
write(io, "number of turbines: ")
writedlm(io, nturbines)
write(io, "rotor diameter: ")
writedlm(io, rotor_diameter[1])
write(io, "starting AEP value (MWh): ")
writedlm(io, aep_wrapper(x_initial)[1]*1E5)
write(io, "Finished in : ")
writedlm(io, clkt)
write(io, "end AEP value (MWh): ")
writedlm(io, aep_wrapper(xopt)[1]*1E5)
write(io, "\noptimal turbine coordinates:\n")
writedlm(io, [turbine_x turbine_y])
write(io, "\neach function call value\n")
writedlm(io, funcalls_AEP)
close(io)

# write results to csv files
dataforcsv_xopt = DataFrame(xopt = xopt)
dataforcsv_funceval = DataFrame(function_value = funcalls_AEP)
dataforcsv_funceval_at_iteration = DataFrame(function_value = iter_AEP)
CSV.write("cs3_xopt.csv", dataforcsv_xopt)
CSV.write("cs3_functionvalue_log.csv", dataforcsv_funceval)
CSV.write("cs3_functionvalue_at_iteration_log.csv", dataforcsv_funceval_at_iteration)

# # write results to a YAML file (not complete)
# dataforyaml = [
#     "optimal_turbine_coordinates" => [turbine_x turbine_y]
#     "function_evaluation_log" => funcalls_AEP
#     "function_evaluation_at_iterations_log" => iter_AEP
#     ]
# YAML.write_file("cs3_results.yml", dataforyaml)