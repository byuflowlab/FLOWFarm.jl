
using Distributed
using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff
using BenchmarkTools


# set up boundary constraint wrapper function
function boundary_wrapper(x, params)

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(params.boundary_center, params.boundary_radius, turbine_x, turbine_y)
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return spacing distances
    return 2.0*params.rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
function aep_wrapper(x, params)

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, params.turbine_z, params.rotor_diameter,
    params.hub_height, params.turbine_yaw, params.ct_models, params.generator_efficiency, params.cut_in_speed,
    params.cut_out_speed, params.rated_speed, params.rated_power, params.windresource, params.power_models, params.model_set,
                rotor_sample_points_y=params.rotor_points_y,rotor_sample_points_z=params.rotor_points_z, hours_per_year=365.0*24.0)
    
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
include("./model_sets/model_set_7_ieacs3.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-12

# set wind farm boundary parameters
boundary_center = [6000.0,6000.0]
boundary_radius = 5500.0

# set globals for use in wrapper functions
struct params_struct{MS, AF, F, I, ACTM, WR, APM}
    model_set::MS
    rotor_points_y::AF
    rotor_points_z::AF
    turbine_z::AF
    rotor_diameter::AF
    boundary_center::AF
    boundary_radius::F
    obj_scale::I
    hub_height::AF
    turbine_yaw::AF
    ct_models::ACTM
    generator_efficiency::AF
    cut_in_speed::AF
    cut_out_speed::AF
    rated_speed::AF
    rated_power::AF
    windresource::WR
    power_models::APM
end

params1 = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z,
    rotor_diameter, boundary_center, boundary_radius, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]
xinit = deepcopy(x)
state_aeps = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set;
                rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.0*24.0)

dir_aep = zeros(20)
for i in 1:20
    for j in 1:20
        dir_aep[i] += state_aeps[(i-1)*20 + j]
    end
end

# report initial objective value
println("nturbines: ", nturbines)
println("rotor diameter: ", rotor_diameter[1])
println("starting AEP value (MWh): ", aep_wrapper(x, params1)[1])
# println("frequencies ", windprobabilities)
println("Directional AEP at start: ", dir_aep.*1E-6)

# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# set general lower and upper bounds for design variables
lb = zeros(length(x)) .- boundary_radius .+ boundary_center[1]
ub = zeros(length(x)) .+ boundary_radius .+ boundary_center[2]

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 3
options["Major optimality tolerance"] = 1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary.out"
options["Print file"] = "print.out"

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params1)
@everywhere aep_wrapper(x) = aep_wrapper(x, params1)
boundary_wrapper(x) = boundary_wrapper(x, params1)

# # run and time optimization
# t1 = time()
# xopt, fopt, info = snopt(wind_farm_opt, x, lb, ub, options)
# t2 = time()
# clkt = t2-t2

# # print optimization results
# println("Finished in : ", clkt, " (s)")
# println("info: ", info)
# println("end objective value: ", aep_wrapper(xopt))

# # extract final turbine locations
# turbine_x = copy(xopt[1:nturbines])
# turbine_y = copy(xopt[nturbines+1:end])

# # add final turbine locations to plot
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
# end

# # add wind farm boundary to plot
# plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))

# # set up and show plot
# axis("square")
# xlim(-boundary_radius+boundary_center[1],boundary_radius+boundary_center[1])
# ylim(-boundary_radius+boundary_center[2],boundary_radius+boundary_center[2])
# plt.show()
println("CPUs: ", nprocs())
@benchmark ForwardDiff.jacobian(aep_wrapper, x) setup=(x=xinit)
# @benchmark aep_wrapper(xinit)