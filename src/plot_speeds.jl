using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics


function calc_TI(constant,ai,TI_free,initial,sep,downstream)

    ti_calculation = constant * (1.0/3.0)^ai * TI_free^initial * sep^downstream

    TI = sqrt(ti_calculation^2 + TI_free^2)

    return TI
end


include("data_12_low.jl")

horizontal = range(-200.0,stop=200.0,length=161)
vertical = range(0.0,stop=300.0,length=121)

# subplot(131)
# plot(horizontal,Uinf4)
# ylim([8.0,13.0])
#
# subplot(132)
# plot(horizontal,Uinf7)
# ylim([8.0,13.0])
#
# subplot(133)
# plot(horizontal,Uinf10)
# ylim([8.0,13.0])

# subplot(131)
# plot(horizontal,Vinf4)
# ylim([0,300])
#
# subplot(132)
# plot(horizontal,Vinf7)
# ylim([0,300])
#
# subplot(133)
# plot(horizontal,Vinf10)
# ylim([0,300])

subplot(131)
plot(horizontal,Winf4)
# ylim([0,300])

subplot(132)
plot(horizontal,Winf7)
# ylim([0,300])

subplot(133)
plot(horizontal,Winf10)
# ylim([0,300])


# subplot(131)
# plot(Uinf4_V,vertical)
# ylim([0,300])
# xlim([0,15])
#
# subplot(132)
# plot(Uinf7_V,vertical)
# ylim([0,300])
# xlim([0,15])
#
# subplot(133)
# plot(Uinf10_V,vertical)
# ylim([0,300])
# xlim([0,15])

# subplot(131)
# plot(Vinf4_V,vertical)
# ylim([0,300])
# # xlim([-1,1])
#
# subplot(132)
# plot(Vinf7_V,vertical)
# ylim([0,300])
# # xlim([-1,1])
#
# subplot(133)
# plot(Vinf10_V,vertical)
# ylim([0,300])
# # xlim([-1,1])


# subplot(131)
# plot(Winf4_V,vertical)
# ylim([0,300])
# xlim([-1,1])
#
# subplot(132)
# plot(Winf7_V,vertical)
# ylim([0,300])
# xlim([-1,1])
#
# subplot(133)
# plot(Winf10_V,vertical)
# ylim([0,300])
# xlim([-1,1])



include("model.jl")


TI_free = 0.046
wind_speed = 10.
ka = 0.38
kb = 0.004
initial = 0.313
constant = 1.931
ai = 0.435
downstream = -0.855
# alpha_star = 2.
# beta_star = 0.1
alpha_star = 2.32
beta_star = 0.154


windspeeds = [wind_speed]
turbine_inflow_velcities = zeros(nturbines) .+ wind_speed
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])



"""U speeds"""



#
# sep = 4.0
# turbine_x = [0.0,sep*rotor_diameter]
# turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# ky = ka*turbulence_intensity + kb
# kz = ka*turbulence_intensity + kb
# horizontal_spread_rate = ky
# vertical_spread_rate = kz
# wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
# wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
# ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)
#
# turbine_y = [0.0,0.0]
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
# pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
#
#
# npts = length(horizontal)
# point_velocities = zeros(npts)
# for i in 1:npts
#       loc = [sep*rotor_diameter, horizontal[i], hub_height]
#       point_velocities[i] = FlowFarm.point_velocity(loc, ms2, pd2)
# end
#
# subplot(131)
# plot(horizontal,point_velocities)
#
#
# sep = 7.0
# turbine_x = [0.0,sep*rotor_diameter]
# turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# ky = ka*turbulence_intensity + kb
# kz = ka*turbulence_intensity + kb
# horizontal_spread_rate = ky
# vertical_spread_rate = kz
# wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
# wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
# ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)
#
# turbine_y = [0.0,0.0]
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
# pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
#
#
# npts = length(horizontal)
# point_velocities = zeros(npts)
# for i in 1:npts
#       loc = [sep*rotor_diameter, horizontal[i], hub_height]
#       point_velocities[i] = FlowFarm.point_velocity(loc, ms2, pd2)
# end
#
# subplot(132)
# plot(horizontal,point_velocities)
#
#
# sep = 10.0
# turbine_x = [0.0,sep*rotor_diameter]
# turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# ky = ka*turbulence_intensity + kb
# kz = ka*turbulence_intensity + kb
# horizontal_spread_rate = ky
# vertical_spread_rate = kz
# wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
# wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
# ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)
#
# turbine_y = [0.0,0.0]
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
# pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
#
#
# npts = length(horizontal)
# point_velocities = zeros(npts)
# for i in 1:npts
#       loc = [sep*rotor_diameter, horizontal[i], hub_height]
#       point_velocities[i] = FlowFarm.point_velocity(loc, ms2, pd2)
# end
#
#
#
# subplot(133)
# plot(horizontal,point_velocities)


"""V or W speeds"""

omega_func = Akima(speeds, omegas)
Omega_rpm = omega_func(wind_speed) #in rpm
Omega = Omega_rpm*0.10471975512 #convert to rad/s
println("Omega: ", Omega)
println("Omega: ", Omega_rpm)

pitch_func = Akima(speeds, pitches)
pitch_deg = pitch_func(wind_speed) #in degrees
pitch = pitch_deg*pi/180.0
println("pitch: ", pitch)


"""vertical"""
# # 90 deg
# azimuth = deg2rad(0.)
#
# z_heights = 90.0 .+ r
# velocities = wind_speed.*(z_heights./90.0).^shearexponent
# rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
# sections = CCBlade.Section.(r,chord,theta,airfoils)
# op = ff.distributed_velocity_op.(velocities, Omega, r, precone, yaw, deg2rad(5.0), azimuth, 1.225)
# out = CCBlade.solve.(Ref(rotor), sections, op)
# V = 2.0.*out.v #farfield swirl
# W = zeros(length(r))
#
# subplot(131)
# scatter(V,z_heights)
#
# subplot(132)
# scatter(V,z_heights)
#
# subplot(133)
# scatter(V,z_heights)
#
# azimuth = deg2rad(180.)
#
# z_heights = 90.0 .- r
# velocities = wind_speed.*(z_heights./90.0).^shearexponent
# rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
# sections = CCBlade.Section.(r,chord,theta,airfoils)
# op = ff.distributed_velocity_op.(velocities, Omega, r, precone, yaw, deg2rad(5.0), azimuth, 1.225)
# out = CCBlade.solve.(Ref(rotor), sections, op)
# V = -2.0.*out.v #farfield swirl
# W = zeros(length(r))
#
# subplot(131)
# scatter(V,z_heights)
#
# subplot(132)
# scatter(V,z_heights)
#
# subplot(133)
# scatter(V,z_heights)

"""horizontal"""
# 90 deg
azimuth = deg2rad(90.)

z_heights = 90.0 .+ zeros(length(r))
velocities = wind_speed.+ zeros(length(r))
rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
sections = CCBlade.Section.(r,chord,theta,airfoils)
op = ff.distributed_velocity_op.(velocities, Omega, r, precone, yaw, deg2rad(5.0), azimuth, 1.225)
out = CCBlade.solve.(Ref(rotor), sections, op)
W = -2.0.*out.v #farfield swirl
V = zeros(length(r))

subplot(131)
scatter(-r,W)

subplot(132)
scatter(-r,W)

subplot(133)
scatter(-r,W)

azimuth = deg2rad(270.)

z_heights = 90.0 .+ zeros(length(r))
velocities = wind_speed.+ zeros(length(r))
rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
sections = CCBlade.Section.(r,chord,theta,airfoils)
op = ff.distributed_velocity_op.(velocities, Omega, r, precone, yaw, deg2rad(5.0), azimuth, 1.225)
out = CCBlade.solve.(Ref(rotor), sections, op)
W = 2.0.*out.v #farfield swirl
V = zeros(length(r))

subplot(131)
scatter(r,W)

subplot(132)
scatter(r,W)

subplot(133)
scatter(r,W)

show()
