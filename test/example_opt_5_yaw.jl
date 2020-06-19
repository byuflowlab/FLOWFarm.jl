using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff

# set up objective wrapper function
function aep_wrapper(x, params)
    # include relevant globals
    turbine_z = params.turbine_z
    turbine_x = params.turbine_x
    turbine_y = params.turbine_y
    rotor_diameter = params.rotor_diameter
    hub_height = params.hub_height
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

    # get number of turbines
    nturbines = length(x)

    # extract yaw of turbines from design variables vector
    turbine_yaw = x

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z, hours_per_year=365.0*24.0)
    
    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt(x)

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)

    # set fail flag to false
    fail = false

    c = []
    dcdx = []
    
    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# import model set with wind farm and related details
include("./model_sets/model_set_8_shiloh.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set globals for use in wrapper functions
struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_x
    turbine_y
    turbine_z
    rotor_diameter
    obj_scale
    hub_height
    ct_models
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    windresource
    power_models
end

params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_x, turbine_y, turbine_z,
    rotor_diameter, obj_scale, hub_height, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x = zeros(nturbines)

# report initial objective value
println("Nturbines: ", nturbines)
println("Rotor diameter: ", rotor_diameter[1])
t1 = time()
println("Starting AEP value (GWh): ", aep_wrapper(x, params)*1e-9/obj_scale)
t2 = time()
println("AEP calc took ", t2-t1, "sec.")
# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 0
options["Major optimality tolerance"] = 1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "snopt_yaw_summary.out"
options["Print file"] = "snopt_yaw_print.out"
lb = zeros(nturbines) .- 30.0*pi/180.0
ub = zeros(nturbines) .+ 30.0*pi/180.0

# generate wrapper function surrogate
aep_wrapper(x) = aep_wrapper(x, params)

# set up holders for optimized yaw values for each state
optyaw = zeros((nstates,nturbines))
diraep = zeros(nstates)
# run and time optimization
t1 = time()
for i in 1:nstates
    println("Optimizing for state: ", i)
    println("Direction: ", winddirections[i])
    println("Wind speed: ", windspeeds[i])
    println("Probability: ", windprobabilities[i])
    params.windresource.wind_directions[1] = winddirections[i]
    params.windresource.wind_speeds[1] = windspeeds[i]
    params.windresource.wind_probabilities[1] = windprobabilities[i]
    params.windresource.measurement_heights[1] = measurementheight[i]
    params.windresource.ambient_tis[1] = ambient_tis[i]
    x = zeros(nturbines)
    xopt, fopt, info = snopt(wind_farm_opt, x, lb, ub, options)

    diraep[i] = fopt
    for j in 1:nturbines
        optyaw[i,j] = xopt[j]
    end
end

t2 = time()
clkt = t2-t2

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end AEP value: ", sum(diraep))

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# set up and show plot
axis("square")
plt.show()
