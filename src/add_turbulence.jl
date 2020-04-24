using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics
using NPZ

const ff=FlowFarm

function get_peaks(A)
    peaks = [A[1]]
    if A[2] > A[1]
        slope = "p"
    elseif A[2] < A[1]
        slope = "m"
    else
        slope = "z"
    end
    for i = 1:length(A)-2
        ind = i+1
        if slope == "p"
            if A[ind+1] < A[ind]
                peaks = append!(peaks,A[ind])
                slope = "m"
    end

end


function delete_repeats(A,t)
    B = [A[1]]
    T = [t[1]]
    for i = 1:length(A)-1
        if A[i+1] != A[i]
            B = append!(B,A[i])
            T = append!(T,t[i])
        end
    end
    return B,T
end


function calc_TI(constant,ai,TI_free,initial,sep,downstream)

    ti_calculation = constant * (1.0/3.0)^ai * TI_free^initial * sep^downstream

    TI = sqrt(ti_calculation^2 + TI_free^2)

    return TI
end


function get_moments(out,Rhub,Rtip,r,az)
    loads_flap = out.Np/1000.
    loads_edge = out.Tp/1000.
    #
    #approximate loads at r = Rhub
    dL_flap = loads_flap[2]-loads_flap[1]
    dL_edge = loads_edge[2]-loads_edge[1]
    dr = r[2]-r[1]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_flap[1] + m_flap*(Rhub-r[1])
    Lhub_edge = loads_edge[1] + m_edge*(Rhub-r[1])

    loads_flap = append!([Lhub_flap],loads_flap)
    loads_edge = append!([Lhub_edge],loads_edge)
    r = append!([Rhub],r)

    #approximate loads at r = Rtip
    dL_flap = loads_flap[end]-loads_flap[end-1]
    dL_edge = loads_edge[end]-loads_edge[end-1]
    dr = r[end]-r[end-1]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_flap[end] + m_flap*(Rtip-r[end])
    Lhub_edge = loads_edge[end] + m_edge*(Rtip-r[end])

    loads_flap = append!(loads_flap,[Lhub_flap])
    loads_edge = append!(loads_edge,[Lhub_edge])
    r = append!(r,[Rtip])

    M_flap = trapz(r,loads_flap.*(r.-Rhub))
    M_edge = trapz(r,loads_edge.*(r.-Rhub))

    # add gravity loads
    blade_mass=17536.617
    blade_cm=20.650
    grav=9.81
    M_edge += sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.
    return M_flap, M_edge
end


function multiple_components_op(U, V, W, Omega, r, precone, yaw, tilt, azimuth, rho, mu=1.81206e-05, asound=1.0)

    sy = sin(yaw)
    cy = cos(yaw)
    st = sin(tilt)
    ct = cos(tilt)
    sa = sin(azimuth)
    ca = cos(azimuth)
    sc = sin(precone)
    cc = cos(precone)

    magnitude = sqrt.(U.^2 + V.^2 + W.^2)
    x0 = U./magnitude
    y0 = V./magnitude
    z0 = W./magnitude

    #first
    x1 = x0.*cy + y0.*sy
    y1 = -x0.*sy + y0.*cy
    z1 = z0

    #second
    x2 = -z1.*st + x1.*ct
    y2 = y1
    z2 = z1.*ct + x1.*st

    #third
    x3 = x2
    y3 = y2.*ca + z2.*sa
    z3 = -y2.*sa + z2.*ca

    #fourth
    x4 = z3.*sc + x3.*cc
    y4 = y3
    z4 = z3.*cc - x3.*sc

    # transform wind to blade c.s.
    Vwind_x = magnitude.*x4
    Vwind_y = magnitude.*y4
    Vwind_z = magnitude.*z4

    # coordinate in azimuthal coordinate system
    z_az = r.*cos(precone)

    # wind from rotation to blade c.s.
    Vrot_x = -Omega*sc
    Vrot_y = Omega.*z_az

    # total velocity
    Vx = Vwind_x .+ Vrot_x
    Vy = Vwind_y + Vrot_y

    # operating point
    return OperatingPoint(Vx, Vy, rho, mu, asound)

end

include("FAST_data.jl")

turb = "low"
ws = 11.0

u_turb = npzread("flowfields_lowTI/speeds11.npy").+ws
time = range(0.0,stop=600.0,length=length(u_turb))
u_turb,time = delete_repeats(u_turb,time)

# include("model.jl")
#
# TI_free = 0.046
# TI = "low"
#
#
# ka = 0.38
# kb = 0.004
#
# initial = 0.313
# constant = 1.931
# ai = 0.435
# downstream = -0.855
# alpha_star = 2.32
# beta_star = 0.154
# turb_index = 2
#
# points_x = [0.69,0,-0.69,0]
# points_y = [0,0.69,0,-0.69]
#
# fos = 1.15
#
# Nlocs = 20
# xlocs = zeros(Nlocs)
# ylocs = zeros(Nlocs)
# angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
# root_rad=3.542/2.
# for i in 1:Nlocs
#     xlocs[i] = cos(angles[i])*root_rad
#     ylocs[i] = sin(angles[i])*root_rad
# end
#
# rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
# sections = CCBlade.Section.(r,chord,theta,airfoils)
#
# # off = [-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
# # off = range(-1.5,stop=1.5,length=50)
# off = -100.0
# dams = zeros(length(off))
#
#
# zero = true
#
# sep = 4.0
# turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# ky = ka*turbulence_intensity + kb
# kz = ka*turbulence_intensity + kb
# horizontal_spread_rate = ky
# vertical_spread_rate = kz
# wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
# wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
# ms = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)
#
#
# turbine_x = [0.0,sep*rotor_diameter]
# init_inflow_velcities = zeros(length(turbine_x)).+ws
#
# omega_func = Akima(speeds, omegas)
# pitch_func = Akima(speeds, pitches)
#
# tilt = deg2rad(5.0)
# rho = 1.225
#
# offset = off*rotor_diameter
# turbine_y = [0.0,offset]
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, init_inflow_velcities, zeros(nturbines))
#
# windspeeds = [ws]
# turbine_inflow_velcities = zeros(nturbines) .+ ws
# windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
#
# pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
# ff.turbine_velocities_one_direction!(points_x, points_y, ms, pd)
# turbine_inflow_velcities = pd.wind_farm_states[1].turbine_inflow_velcities
# Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
# Omega = Omega_rpm*0.10471975512 #convert to rad/s
# pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
# pitch = pitch_deg*pi/180.0
#
# flap = zeros(length(u_turb))
# edge = zeros(length(u_turb))
#
# for i = 1:length(u_turb)
#     println(i)
#     windspeeds = [u_turb[i]]
#     turbine_inflow_velcities = zeros(nturbines) .+ u_turb[i]
#     windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
#     pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
#
#     """find the azimuth"""
#     az = (Omega*time[i])%(2.0*pi)
#     U = ff.get_speeds(turbine_x,turbine_y,turb_index,hubHt,r,yaw,az,ms,pd)
#     op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)
#     out = CCBlade.solve.(Ref(rotor), sections, op)
#     flap[i],edge[i] = get_moments(out,Rhub,Rtip,r,az)
# end
#
# Nlocs = 20
# xlocs = zeros(Nlocs)
# ylocs = zeros(Nlocs)
# angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
# root_rad=3.542/2.
# for i in 1:Nlocs
#     xlocs[i] = cos(angles[i])*root_rad
#     ylocs[i] = sin(angles[i])*root_rad
# end
#
# damage = 0.0
# for i  = 1:Nlocs
#     sigma = ff.calc_moment_stress.(edge,flap,xlocs[i],ylocs[i])
#
# end
#
# #
# #     damage = 0.
# #
# #     su = 70000. # hard coded in now, maybe worth it to let it be passed in
# #     m = 10.
# #     years = 25.
# #
# #     for j in 1:Nlocs
# #
# #         stress0 = ff.calc_moment_stress(M_edge0,M_flap0,xlocs[j],ylocs[j])
# #         stress90 = ff.calc_moment_stress(M_edge90,M_flap90,xlocs[j],ylocs[j])
# #         stress180 = ff.calc_moment_stress(M_edge180,M_flap180,xlocs[j],ylocs[j])
# #         stress270 = ff.calc_moment_stress(M_edge270,M_flap270,xlocs[j],ylocs[j])
# #
# #         stresses = [stress0,stress90,stress180,stress270]
# #         cycle = [maximum(stresses),minimum(stresses)]
# #         mean = sum(cycle)/2.
# #         alternate = (maximum(cycle)-minimum(cycle))/2.
# #         effective = alternate/(1.0-mean/su)
# #         effective = effective*omega_mult
# #
# #         Nfail = (su/(effective*fos))^m #mLife
# #         omega_rpm = Omega/0.10471975512
# #         nCycles = omega_rpm*60.0*24.0*365.25*years
# #         d = nCycles/Nfail
# #         if d > damage
# #                 damage = d
# #         end
# #
# #     end
# #
# #     dams[i] = damage
# #     end
# # end
#
# figure(1)
# plot(time,flap)
# figure(2)
# plot(time,edge)
# #
# # # plot(off,dams)
# # # title("SOWFA")
# # # # ylim([0,1.5])
# # # println("damage: ", dams)
# # #
# # # FS,FD = fastdata(turb,ws,sep)
# # # scatter(FS,FD)
# # #
# # # legend(["mean SOWFA velocity","SOWFA+FAST"])
# show()
# # #
# # #
# # # println("finished")
