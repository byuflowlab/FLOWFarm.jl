using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics


function distributed_velocity_op(V, Omega, r, precone, yaw, tilt, azimuth, rho, mu=1.81206e-05, asound=1.0)

    sy = sin(yaw)
    cy = cos(yaw)
    st = sin(tilt)
    ct = cos(tilt)
    sa = sin(azimuth)
    ca = cos(azimuth)
    sc = sin(precone)
    cc = cos(precone)

    # coordinate in azimuthal coordinate system
    x_az = -r.*sin(precone)
    z_az = r.*cos(precone)

    # get section heights in wind-aligned coordinate system
    heightFromHub = (z_az*ca)*ct - x_az*st

    # transform wind to blade c.s.
    Vwind_x = V * ((cy*st*ca + sy*sa)*sc + cy*ct*cc)
    Vwind_y = V * (cy*st*sa - sy*ca)

    # wind from rotation to blade c.s.
    Vrot_x = -Omega*sc
    Vrot_y = Omega*z_az

    # total velocity
    Vx = Vwind_x + Vrot_x
    Vy = Vwind_y + Vrot_y

    # operating point
    return OperatingPoint(Vx, Vy, rho, mu, asound)

end

function calculate_root_moment(r,velocity,chord,theta,af,Rhub,Rtip,hubHt,Omega,pitch,azimuth,
                B=3,rho=1.225,precone=0.,yaw=0.,tilt=0.,blade_mass=17536.617,blade_cm=20.650,grav=9.81)

        rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
        sections = CCBlade.Section.(r,chord,theta,airfoils)
        op = distributed_velocity_op.(velocity, Omega, r, precone, yaw, tilt, azimuth, rho)
        out = CCBlade.solve.(Ref(rotor), sections, op)

        loads_flap = out.Np/1000.
        loads_edge = out.Tp/1000.
        # println("flap: ", loads_flap)
        # println("edge: ", loads_edge)

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
        M_edge += sin(azimuth)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.

        return M_flap,M_edge

end

function find_xyz_simple(x_hub,y_hub,z_hub,r,yaw,azimuth)

        sy = sin(yaw)
        cy = cos(yaw)
        sa = sin(azimuth)
        ca = cos(azimuth)

        x_locs = x_hub .- r.*(sy*sa)
        y_locs = y_hub .- r.*(cy*sa)
        z_locs = z_hub .+ r.*(ca)

        return x_locs, y_locs, z_locs

end


"""
    calc_moment_stress(mx,my,dx,dy,Rcyl=1.771,tcyl=0.06)

Calculate the stress on a cylinder caused by bending moments as a specific location.
"""
function calc_moment_stress(mx,my,dx,dy,Rcyl=1.771,tcyl=0.06)

        I = 0.25*pi*(Rcyl^4-(Rcyl-tcyl)^4)
        stress = mx*dy/I + my*dx/I

        return stress

end


function get_speeds(turbineX,turbineY,turb_index,hubHt,r,yaw,azimuth,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
        x_locs, y_locs, z_locs = find_xyz_simple(turbineX[turb_index],turbineY[turb_index],hubHt,r,yaw,azimuth)
        npts = length(x_locs)
        point_velocities = zeros(npts)
            for i in 1:npts
                loc = [x_locs[i], y_locs[i], z_locs[i]]
                point_velocities[i] = FlowFarm.point_velocity(loc, windfarm, windfarmstate, windresource, wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, 0)
            end
        return point_velocities
end


function get_damage(turbineX,turbineY,turb_index,TSR,pitch,free_speed,Rhub,Rtip,r,chord,theta,af,hubHt,
                        windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel,
                        rho=1.225,mu=1.81206e-05,B=3,precone=0,yaw=0,TI=0.11,Nlocs=10,fos=1.25,shearExp=0.15,root_rad=3.542/2.)

        """get the TSR"""

        """get the root moments at 0, 90, 180, 270 degrees (need more?)"""
        # 90 deg
        azimuth = deg2rad(0.)
        point_velocities = get_speeds(turbineX,turbineY,turb_index,hubHt,r,yaw,azimuth,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
        # println("0")
        Omega = TSR*turbine_inflow_velcities[turb_index]/Rtip
        m0 = calculate_root_moment.(Ref(r),Ref(point_velocities),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,azimuth)

        azimuth = deg2rad(90.)
        point_velocities = get_speeds(turbineX,turbineY,turb_index,hubHt,r,yaw,azimuth,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
        # println("90")
        # println(point_velocities)
        Omega = TSR*turbine_inflow_velcities[turb_index]/Rtip
        m90 = calculate_root_moment.(Ref(r),Ref(point_velocities),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,azimuth)

        azimuth = deg2rad(180.)
        point_velocities = get_speeds(turbineX,turbineY,turb_index,hubHt,r,yaw,azimuth,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
        # println("180")
        Omega = TSR*turbine_inflow_velcities[turb_index]/Rtip
        m180 = calculate_root_moment.(Ref(r),Ref(point_velocities),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,azimuth)

        azimuth = deg2rad(270.)
        point_velocities = get_speeds(turbineX,turbineY,turb_index,hubHt,r,yaw,azimuth,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
        # println("270")
        # println(point_velocities)
        Omega = TSR*turbine_inflow_velcities[turb_index]/Rtip
        m270 = calculate_root_moment.(Ref(r),Ref(point_velocities),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,azimuth)

        # println("0: ", m0)
        # println("90: ", m90)
        # println("180: ", m180)
        # println("270: ", m270)

        angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
        xlocs = zeros(Nlocs)
        ylocs = zeros(Nlocs)
        for i in 1:Nlocs
            xlocs[i] = cos(angles[i])*root_rad
            ylocs[i] = sin(angles[i])*root_rad
        end

        damage = 0.

        su = 70000. # hard coded in now, maybe worth it to let it be passed in
        m = 10.
        years = 25.

        for i in 1:Nlocs

            stress0 = calc_moment_stress(m0[2],m0[1],xlocs[i],ylocs[i])
            stress90 = calc_moment_stress(m90[2],m90[1],xlocs[i],ylocs[i])
            stress180 = calc_moment_stress(m180[2],m180[1],xlocs[i],ylocs[i])
            stress270 = calc_moment_stress(m270[2],m270[1],xlocs[i],ylocs[i])

            stresses = [stress0,stress90,stress180,stress270]
            cycle = [maximum(stresses),minimum(stresses)]
            mean = sum(cycle)/2.
            # println("i: ", i)
            # println("mean: ", mean)
            alternate = (maximum(cycle)-minimum(cycle))/2.
            # println("alternate: ", alternate)
            effective = alternate/(1.0-mean/su)
            # println("effective: ", effective)

            Nfail = (su/(effective*fos))^m #mLife
            # println("Nfail: ", Nfail)
            nCycles = Omega*60.0*24.0*365.25*years
            d = nCycles/Nfail
            if d > damage
                    damage = d
                    # println(damage)
            end

        end

        return damage

end
#         actual_speed = get_eff_turbine_speeds(turbineX, turbineY, free_speed,TI=TI)[turb_index]
#         if isinstance(TSR_TSR,float):
#                 TSR = TSR_TSR
#         else:
#                 f_TSR = interp1d(TSR_speeds, TSR_TSR, kind='cubic')
#                 TSR = f_TSR(actual_speed)
#
#
#         Omega = TSR*actual_speed/Rtip
#         Omega = Omega*9.5492965964254 #rad/s to RPM
#
#         f_pitch = interp1d(pitch_speeds, pitch_angle)
#         # print 'speed: ', actual_speed
#         if actual_speed < pitch_speeds[0]:
#                 pitch = 0.
#         else:
#                 pitch = f_pitch(actual_speed)
#         pitch = 0.
#
#         az = 90.
#         x_locs,y_locs,z_locs = findXYZ(turbineX[turb_index],turbineY[turb_index],hubHt,r,yaw_deg,az)
#         speeds, _ = get_speeds(turbineX, turbineY, x_locs, y_locs, z_locs, free_speed, TI=TI,shearExp=shearExp)
#         flap90,edge90 = calc_moment(speeds,Rhub,r,chord,theta,af,Rhub,Rtip,B,rho,mu,precone,hubHt,nSector,Omega,pitch,azimuth=az)
#         az = 270.
#         x_locs,y_locs,z_locs = findXYZ(turbineX[turb_index],turbineY[turb_index],hubHt,r,yaw_deg,az)
#         speeds, _ = get_speeds(turbineX, turbineY, x_locs, y_locs, z_locs, free_speed, TI=TI,shearExp=shearExp)
#         flap270,edge270 = calc_moment(speeds,Rhub,r,chord,theta,af,Rhub,Rtip,B,rho,mu,precone,hubHt,nSector,Omega,pitch,azimuth=az)
#         az = 0.
#         x_locs,y_locs,z_locs = findXYZ(turbineX[turb_index],turbineY[turb_index],hubHt,r,yaw_deg,az)
#         speeds, _ = get_speeds(turbineX, turbineY, x_locs, y_locs, z_locs, free_speed, TI=TI,shearExp=shearExp)
#         flap0,edge0 = calc_moment(speeds,Rhub,r,chord,theta,af,Rhub,Rtip,B,rho,mu,precone,hubHt,nSector,Omega,pitch,azimuth=az)
#         az = 180.
#         x_locs,y_locs,z_locs = findXYZ(turbineX[turb_index],turbineY[turb_index],hubHt,r,yaw_deg,az)
#         speeds, _ = get_speeds(turbineX, turbineY, x_locs, y_locs, z_locs, free_speed, TI=TI,shearExp=shearExp)
#         flap180,edge180 = calc_moment(speeds,Rhub,r,chord,theta,af,Rhub,Rtip,B,rho,mu,precone,hubHt,nSector,Omega,pitch,azimuth=az)
#
#         angles = np.linspace(-np.pi/2,np.pi/2.,Nlocs)
#         xlocs = np.zeros(Nlocs)
#         ylocs = np.zeros(Nlocs)
#         # R = 0.5
#         R = 3.542/2.
#         for i in range(Nlocs):
#             xlocs[i] = np.cos(angles[i])*R
#             ylocs[i] = np.sin(angles[i])*R
#
#         damage = 0.
#
#         su = 70000.
#         # su = 535000.
#         # su = 596000.
#         # su = 800000.
#         m = 10.
#
#         years = 25.
#
#         for i in range(Nlocs):
#
#             stress0 = calc_moment_stress(edge0,flap0,xlocs[i],ylocs[i])/1000.
#             stress90 = calc_moment_stress(edge90,flap90,xlocs[i],ylocs[i])/1000.
#             stress180 = calc_moment_stress(edge180,flap180,xlocs[i],ylocs[i])/1000.
#             stress270 = calc_moment_stress(edge270,flap270,xlocs[i],ylocs[i])/1000.
#
#             # stress90 = calc_moment_stress(edge90,0.,xlocs[i],ylocs[i])/1000.
#             # stress270 = calc_moment_stress(edge270,0.,xlocs[i],ylocs[i])/1000.
#
#             # stress90 = calc_moment_stress(0.,flap90,xlocs[i],ylocs[i])/1000.
#             # stress270 = calc_moment_stress(0.,flap270,xlocs[i],ylocs[i])/1000.
#
#             stress_max = np.max(np.array([stress0,stress90,stress180,stress270]))
#             stress_min = np.min(np.array([stress0,stress90,stress180,stress270]))
#
#             # cycle = np.array([stress90,stress270])
#             cycle = np.array([stress_max,stress_min])
#             mean = np.sum(cycle)/2.
#             alternate = (np.max(cycle)-np.min(cycle))/2.
#
#             effective_stress = alternate/(1.-mean/su)
#
#             Nfail = (su/(effective_stress*fos))**m #mLife
#             nCycles = Omega*60.*24.*365.25*years
#             d = nCycles/Nfail
#             if d > damage:
#                     damage = d
#
#         return damage
#
# end









#
#
#
#
# const ff=FlowFarm
# # # #testing during development
# r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
#               28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
#               56.1667, 58.9000, 61.6333]
# chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
#                   3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
# theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
#                   6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
# # #
# #
# #
# # x_hub = 0.
# # y_hub = 0.
# # z_hub = 90.
# # yaw = pi/8
# # azimuth = range(0.,stop=8*pi,length=100)
# # precone = 0.
# # tilt = 0*pi/6
# #
# # fig = figure()
# # for n = 1:length(azimuth)
# #
# #         # x_locs, y_locs, z_locs = find_xyz(x_hub,y_hub,z_hub,r,yaw,azimuth[n],precone,tilt)
# #         x_locs, y_locs, z_locs = find_xyz_simple(x_hub,y_hub,z_hub,r,yaw,azimuth[n])
# #         clf()
# #
# #         subplot(131)
# #         plot(x_locs,y_locs)
# #         # axis("equal")
# #         xlim(-70.,70.)
# #         ylim(-70.,70.)
# #
# #         subplot(132)
# #         plot(y_locs,z_locs)
# #         # axis("equal")
# #         ylim(20.,160.)
# #         xlim(-70.,70.)
# #
# #         subplot(133)
# #         plot(x_locs,z_locs)
# #         # axis("equal")
# #         ylim(20.,160.)
# #         xlim(-70.,70.)
# #
# #         sleep(0.05)
# #         display(fig)
# #
# # end
# #
# #
# #
# Rhub = 1.5
# Rtip = 63.0
# B = 3
# # turbine = true
# pitch = 0.0
# precone = 2.5*pi/180
#
# af_path = "/Users/ningrsrch/Dropbox/Projects/waked-loads/5MW_AFFiles_julia"
# path1 = af_path*"/Cylinder1.dat"
# path2 = af_path*"/Cylinder2.dat"
# path3 = af_path*"/DU40_A17.dat"
# path4 = af_path*"/DU35_A17.dat"
# path5 = af_path*"/DU30_A17.dat"
# path6 = af_path*"/DU25_A17.dat"
# path7 = af_path*"/DU21_A17.dat"
# path8 = af_path*"/NACA64_A17.dat"
#
# af1 = af_from_files(path1)
# af2 = af_from_files(path2)
# af3 = af_from_files(path3)
# af4 = af_from_files(path4)
# af5 = af_from_files(path5)
# af6 = af_from_files(path6)
# af7 = af_from_files(path7)
# af8 = af_from_files(path8)
#
# af = [af1,af2,af3,af4,af5,af6,af7,af8]
#
# af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]
# airfoils = af[af_idx]
#
# TSR = 5.
# pitch = 0.
# free_speed = 14.
# Rhub = 1.5
# hubHt = 90.
# #
# # turbineX = [0.,500.]
# # turbineY = [0.,100.]
#
#
# #
# # #
#
#
#
# # # rotor = Rotor(Rhub, Rtip, B, turbine, pitch, precone)
# # #
# # # sections = Section.(r,chord,theta,airfoils)
# # #
# # # # Vy = zeros(17)
# # # # Vx = ones(17).*8.
# # # # rho = 1.225
# # # # op = CCBlade.OperatingPoint(Vx, Vy, rho)
# # #
# # # yaw = 0.0*pi/180
# # # tilt = 5.0*pi/180
# # # hubHt = 90.0
# # # shearExp = 0.2
# # #
# # # Vinf = 10.0
# # # tsr = 7.55
# # # rotorR = Rtip*cos(precone)
# # # Omega = Vinf*tsr/rotorR
# # # azimuth = 90.0*pi/180
# # # rho = 1.225
# # # yaw = 0.
# # # tilt = 0.
# # #
# # # # op = windturbine_op.(Vinf, Omega, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
# # # velocity_dist = ones(17).*Vinf
# # # op = distributed_velocity_op.(velocity_dist, Omega, r, precone, yaw, tilt, azimuth, rho)
# # #
# # # # out = solve.(Ref(rotor), sections, op)
# # # azimuth = range(0.,stop=pi*2,length=100)
# # # M = calculate_root_moment.(Ref(r),Ref(velocity_dist),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,B,rho,precone,hubHt,Omega,pitch,azimuth,yaw,tilt)
# # # figure()
# # # plot(azimuth*180/pi, M)
# # # xlabel("azimuth angle")
# # # ylabel("root moments (kN-m)")
#
#
# rotor_diameter = 126.4
# hub_height = 90.0
# ct = 0.689
# cp = 0.42
# generator_efficiency = 0.944
# ai = 1.0/3.0
# wind_speed = 12.0
# air_density = 1.1716  # kg/m^3
# turbine_x = [0.,500.]
# nturbines = length(turbine_x)
# turbine_y = [0.,-50.]
# turbine_z = zeros(nturbines)
# turbine_yaw = zeros(nturbines)
# turbine_ct = zeros(nturbines) .+ ct
# turbine_ai = zeros(nturbines) .+ ai
# winddirections = [270.0*pi/180.0,0.0*pi/180.0]
# windspeeds = [wind_speed,wind_speed]
# windprobabilities = [1.0,1.0]
# measurementheight = [hub_height,hub_height]
# shearexponent = 0.15
# turbine_inflow_velcities = zeros(nturbines) .+ wind_speed
#
# # rotor sample points
# rotor_points_y = [0.0]
# rotor_points_z = [0.0]
#
# ct_model = ff.ConstantCt(ct)
# power_model = ff.ConstantCp([cp], [generator_efficiency])
# wind_shear_model = ff.PowerLawWindShear(shearexponent)
#
# turbine_definitions = [] # or `Vector{Coords}(undef, x)`
# turbine = ff.Turbine(1, [rotor_diameter], [hub_height], [ct_model], [power_model])
# turbine_definitions = [turbine for i in 1:nturbines]
# sorted_turbine_index = zeros(Int64, nturbines)
# for i in 1:nturbines
#     turbine = ff.Turbine(i, [rotor_diameter], [hub_height], [ct_model], [power_model])
#     turbine_definitions[i] = turbine
#     sorted_turbine_index[i] = Int64(i)
# end
#
# windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definitions)
# windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, [0.0])
# windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
#
# alpha = 0.1
# k_star = 0.022 # [1]  p. 525
# turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
# horizontal_spread_rate = k_star
# vertical_spread_rate = k_star
# alpha_star = 2.32 #[1] p. 534
# beta_star = 0.154 #[1] p. 534
# # wakedeficitmodel = ff.JensenTopHat(alpha)
# wakedeficitmodel = ff.GaussYaw(turbulence_intensity, horizontal_spread_rate , vertical_spread_rate, alpha_star, beta_star)
# horizontal_spread_rate = alpha
# wakedeflectionmodel = ff.NoYawDeflection()
# wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
#
# # ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, windfarm, windfarmstate, windresource, wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)
# # println(windfarmstate.turbine_inflow_velcities)
#
# turb_index = 2
#
# get_damage(turbine_x,turbine_y,turb_index,TSR,pitch,free_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
