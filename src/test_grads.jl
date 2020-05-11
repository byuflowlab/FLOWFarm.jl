using Snopt
import ForwardDiff
import ReverseDiff
include("fatigue_model.jl")
include("model.jl")

# function rosenbrock(x)
#     f = (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2
#
#     g1 = -2*(1 - x[1]) + 200*(x[2] - x[1]^2)*-2*x[1]
#     g2 = 200*(x[2] - x[1]^2)
#     println("g1: ", g1[1,1])
#     println("g2: ", g2[1])
#     #
#     # fail = false
#     #
#     # # c = []
#     # # dcdx = []
#     c = x[1]^3+10
#     # dcdx = 3*x[1]^2
#     # println(dcdx[1])
#     # # dcdx[2,2] = 2*x[2]
#
#     # return f, c, g, dcdx, fail
#     return [f;c]
#
# end
#
#
# x0 = [4.0; 4.0]
# # println(rosenbrock(x0))
# J = ForwardDiff.jacobian(rosenbrock,x0)
# println("J: ", J)
# # lb = [-50.0; -50.0]
# # ub = [50.0; 50.0]
# # options = Dict{String, Any}()
# # options["Derivative option"] = 1
# # options["Verify level"] = 1
# # options["Major optimality tolerance"] = 1e-6
# # options["Major iteration limit"] = 1e6
# #
# # xopt, fopt, info = snopt(rosenbrock, x0, lb, ub, options)
# #
# # println("xopt: ", xopt)
# # println("fopt: ", fopt)
# # println("info: ", info)


# function turbulence_function(loc)
#     r = sqrt(loc[2]^2+(loc[3]-90.0)^2)
#     if loc[1] > 10.0
#         if r < 70.0
#             return 0.17
#         else
#             return 0.05
#         end
#     else
#         return 0.05
#     end
# end
#
#


function aep_wrapper(x)
    global ms3
    # global pd3
    global rotor_points_y
    global rotor_points_z
    global turbine_z
    global turbine_definition_ids
    global turbine_definitions
    global ambient_ti
    global windfarmstate_array

    nturbines = Int(length(x)/2)
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]
    turbine_inflow_velcities = zeros(typeof(turbine_x[1]),nturbines)
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcities,
                            zeros(typeof(turbine_x[1]),nturbines), zeros(typeof(turbine_x[1]),nturbines).+ambient_ti,sorted_turbine_index)
    # wfsa = [windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
    #                             windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
    #                             windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,
    #                             windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate,windfarmstate]
    # windfarmstate_array = []
    # windfarmstate_array = Array(typeof(windfarmstate),length(winddirections))
    windfarmstate_array = Vector{typeof(windfarmstate)}(undef, length(winddirections))

    # windfarmstate_array = similar(wfsa)

    for i = 1:length(winddirections)
                  tx, ty = ff.rotate_to_wind_direction(turbine_x, turbine_y, windresource.wind_directions[i])
                  windfarmstate_array[i] = ff.SingleWindFarmState(i, tx, ty, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                            turbine_inflow_velcities, zeros(typeof(turbine_x[1]),nturbines), zeros(typeof(turbine_x[1]),nturbines).+ambient_ti,sorted_turbine_index)
                  # append!(windfarmstate_array,[ff.SingleWindFarmState(i, tx, ty, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                  #           turbine_inflow_velcities, zeros(typeof(turbine_x[1]),nturbines), zeros(typeof(turbine_x[1]),nturbines).+ambient_ti,sorted_turbine_index)])
    end


    # problem_description.wind_farm_states[i].turbine_x[:],problem_description.wind_farm_states[i].turbine_y[:] =
    #         rotate_to_wind_direction(wind_farm.turbine_x, wind_farm.turbine_y, wind_resource.wind_directions[i])
    #
    # problem_description.wind_farm_states[i].sorted_turbine_index[:] = sortperm(problem_description.wind_farm_states[i].turbine_x)
    #
    # turbine_velocities_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
    #     model_set, problem_description, wind_farm_state_id=i)
    #
    # turbine_powers_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
    #     problem_description, wind_farm_state_id=i)
    #
    # state_power = sum(problem_description.wind_farm_states[i].turbine_generators_powers)
    # state_energy[i] = state_power*hours_per_year*wind_probabilities[i]




    # println(windfarmstate[1].turbine_inflow_velcities[1])
    pd3 = ff.WindFarmProblemDescription(windfarm, windresource, windfarmstate_array)
    # pd3.wind_farm.turbine_x[:] = turbine_x
    # pd3.wind_farm.turbine_y[:] = turbine_y
    AEP = ff.calculate_aep(ms3, pd3,rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)/1e9
    # AEP = 1.0
    # println(AEP)
    return [AEP]
    # return [x[1]]
end

global ms3
# global pd3
global rotor_points_y
global rotor_points_z
global turbine_z
global turbine_definition_ids
global turbine_definitions
global ambient_ti
global windfarmstate_array




# include("model_set_3.jl")
# modelAEP = ff.calculate_aep(ms3, pd3,rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)/1e9

include("model_set_3.jl")
# windfarmstate_array = Vector{Any}(undef, length(winddirections))

x = [copy(turbine_x);copy(turbine_y)]
# x = [0.0;0.0]

# scatter(pd3.wind_farm.turbine_x,pd3.wind_farm.turbine_y)
t1 = time()
# config = ForwardDiff.GradientConfig(aep_wrapper, x, ForwardDiff.Chunk{4}())
# Jfor = ForwardDiff.jacobian(aep_wrapper,x,config)
Jfor = ForwardDiff.jacobian(aep_wrapper,x)

# out = similar(x)
# cfg4 = ForwardDiff.GradientConfig(aep_wrapper, x, ForwardDiff.Chunk{1}())
# Jfor = ForwardDiff.gradient!(out,aep_wrapper,x,cfg4)
println("Forward time: ", time()-t1)
# println("forward diff: ", J)


# include("model_set_3.jl")
# x = [copy(turbine_x);copy(turbine_y)]
# t1 = time()
# Jrev = ReverseDiff.jacobian(aep_wrapper,x)
# println("Reverse time: ", time()-t1)
# println("reverse diff: ", J)


#FD
# include("model_set_3.jl")
x = [copy(turbine_x);copy(turbine_y)]
step = 1E-4
grad = zeros(length(x))
t1 = time()
AEP_orig = aep_wrapper(x)
println("call time: ", time()-t1)
println("AEP: ", AEP_orig)
for i = 1:length(x)
    xtemp = copy(x)
    xtemp[i] += step
    AEP = aep_wrapper(xtemp)
    grad[i] = (AEP[1]-AEP_orig[1])/step
end
#
# println(Jrev)
println("FD: ", maximum((grad .- Jfor[:])./grad))
# println("RD: ", maximum((grad .- Jrev[:])./grad))

# println(grad)
# println(Jfor)
# println(Jrev)
