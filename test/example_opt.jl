using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff
import ReverseDiff



function boundary_wrapper(x)
    global boundary_center
    global boundary_radius
    global scale
    nturbines = Int(length(x)/2)
    turbine_x = x[1:nturbines].*scale
    turbine_y = x[nturbines+1:end].*scale
    return ff.circle_boundary(boundary_center,boundary_radius,turbine_x,turbine_y)
end

function spacing_wrapper(x)
    global rotor_diameter
    global scale
    nturbines = Int(length(x)/2)
    turbine_x = x[1:nturbines].*scale
    turbine_y = x[nturbines+1:end].*scale
    return 2.0*rotor_diameter .- ff.turbine_spacing(turbine_x,turbine_y)
end

function aep_wrapper(x)
    global ms6
    global rotor_points_y
    global rotor_points_z
    global turbine_z
    global turbine_definition_ids
    global turbine_definitions
    global ambient_ti
    global scale

    nturbines = Int(length(x)/2)
    turbine_x = x[1:nturbines].*scale
    turbine_y = x[nturbines+1:end].*scale
    turbine_inflow_velcities = zeros(typeof(turbine_x[1]),nturbines)
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities,
                            zeros(typeof(turbine_x[1]),nturbines), zeros(typeof(turbine_x[1]),nturbines).+ambient_ti,sorted_turbine_index)

    windfarmstate_array = Vector{typeof(windfarmstate)}(undef, length(winddirections))

    for i = 1:length(winddirections)
                  tx, ty = ff.rotate_to_wind_direction(turbine_x, turbine_y, windresource.wind_directions[i])
                  windfarmstate_array[i] = ff.SingleWindFarmState(i, tx, ty, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                            turbine_inflow_velcities, zeros(typeof(turbine_x[1]),nturbines), zeros(typeof(turbine_x[1]),nturbines).+ambient_ti,sorted_turbine_index)
    end



    pd3 = ff.WindFarmProblemDescription(windfarm, windresource, windfarmstate_array)
    AEP = ff.calculate_aep(ms6, pd6, rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)/1e11
    return [AEP]
end

function wind_farm_opt(x)

    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)
    # println("objective gradient: ", dAEP_dx)
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper,x)

    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper,x)
    #
    c = [spacing_con;boundary_con]
    dcdx = [ds_dx;db_dx]

    # c = boundary_con
    # dcdx = db_dx

    fail = false
    return AEP, c, dAEP_dx, dcdx, fail
end

global ms6
global rotor_points_y
global rotor_points_z
global turbine_z
global turbine_definition_ids
global turbine_definitions
global ambient_ti
global rotor_diameter
global boundary_center
global boundary_radius
global scale

scale = 1.

boundary_center = [0.0,0.0]
boundary_radius = 300.0

include("./model_sets/model_set_6.jl")
x = [copy(turbine_x);copy(turbine_y)]./scale

println("start: ", aep_wrapper(x))

turbine_x = copy(x[1:nturbines]).*scale
turbine_y = copy(x[nturbines+1:end]).*scale
area = zeros(length(turbine_x)) .+ pi*(rotor_diameter/2.0)^2
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter/2.0, fill=false, color="C0"))
end

lb = zeros(length(x)) .-boundary_radius
ub = zeros(length(x)) .+boundary_radius
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 3
options["Major optimality tolerance"] = 1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary.out"
options["Print file"] = "print.out"

t1 = time()
xopt, fopt, info = snopt(wind_farm_opt, x, lb, ub, options)
println("Finished in : ", time()-t1, "s")
# println("xopt: ", xopt)
# println("fopt: ", fopt)
println("info: ", info)

println("finish: ", aep_wrapper(xopt))

turbine_x = copy(xopt[1:nturbines]).*scale
turbine_y = copy(xopt[nturbines+1:end]).*scale

for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter/2.0, fill=false,color="C1"))
end
plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))

axis("square")
xlim(-boundary_radius-200,boundary_radius+200)
ylim(-boundary_radius-200,boundary_radius+200)
