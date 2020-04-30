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
#     dcdx = 3*x[1]^2
#     println(dcdx[1])
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


function turbulence_function(loc)
    r = sqrt(loc[2]^2+(loc[3]-90.0)^2)
    if loc[1] > 10.0
        if r < 70.0
            return 0.17
        else
            return 0.05
        end
    else
        return 0.05
    end
end


function damage_wrapper(x)
    global turbine_x
    global turbine_z
    global turbine_definition_ids
    global turbine_definitions
    global turbine_yaw
    global turbine_ct
    global turbine_ai
    global sorted_turbine_index
    global init_inflow_velcities
    global nturbines
    global windresource
    global ms
    global nCycles
    global az_arr
    global turb_samples
    global points_x
    global points_y
    global omega_func
    global pitch_func
    global turbulence_function
    global r
    global rotor
    global sections
    global Rhub
    global Rtip
    global rotor_diameter

    turbine_y = [0.0,rotor_diameter.*x[1]]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, init_inflow_velcities, zeros(nturbines))
    pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    # println("turbine_y: ,", turbine_y)
    # println(turbine_y[1])
    # println(turbine_y[2])
    total_damage = get_total_farm_damage(ms,pd,nCycles,az_arr,
        turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,r,rotor,sections,Rhub,Rtip)
    # fail = false
    # return [total_damage[2]]
    return [7.0]

end

global turbine_x
global turbine_z
global turbine_definition_ids
global turbine_definitions
global turbine_yaw
global turbine_ct
global turbine_ai
global sorted_turbine_index
global init_inflow_velcities
global nturbines
global windresource
global ms
global nCycles
global az_arr
global turb_samples
global points_x
global points_y
global omega_func
global pitch_func
global turbulence_function
global r
global rotor
global sections
global Rhub
global Rtip
global rotor_diameter

turb = "low"
ws = 11.0
windspeeds = [ws]
TI_free = 0.046

windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])

sep = 4.0

ka = 0.38
kb = 0.004

initial = 0.313
constant = 1.931
ai = 0.435
downstream = -0.855
alpha_star = 2.32
beta_star = 0.154

points_x = [0.69,0,-0.69,0]
points_y = [0,0.69,0,-0.69]

fos = 1.15

Nlocs = 20
xlocs = zeros(Nlocs)
ylocs = zeros(Nlocs)
angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
root_rad=3.542/2.
for i in 1:Nlocs
    xlocs[i] = cos(angles[i])*root_rad
    ylocs[i] = sin(angles[i])*root_rad
end

rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
sections = CCBlade.Section.(r,chord,theta,airfoils)

off = -0.6

turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
nturbines = length(turbine_x)
init_inflow_velcities = zeros(length(turbine_x)).+ws

omega_func = Akima(speeds, omegas)
pitch_func = Akima(speeds, pitches)

tilt = deg2rad(5.0)
rho = 1.225

nCycles = 400
naz = 2
az_arr = [pi/2,3*pi/2]

turb_samples = randn(naz*nCycles)
x0 = [off;8.0]

J = ForwardDiff.jacobian(damage_wrapper,x0)
# println("J: ", J)
#
#
# println(damage_wrapper(x0))
# println("done")
