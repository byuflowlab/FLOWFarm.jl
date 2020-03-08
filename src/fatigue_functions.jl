using CCBlade
using PyPlot
using FLOWMath

function distributed_velocity_op(V, Omega, r, precone, yaw, tilt, azimuth, rho,
                mu=1.0, asound=1.0)

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

function calculate_root_moment(r,velocity,chord,theta,af,Rhub,Rtip,B,rho,precone,
                hubHt,Omega,pitch,azimuth,yaw,tilt,blade_mass=17536.617,blade_cm=20.650,grav=9.81)


        rotor = Rotor(Rhub, Rtip, B, turbine, pitch, precone)
        sections = Section.(r,chord,theta,airfoils)
        op = distributed_velocity_op.(velocity, Omega, r, precone, yaw, tilt, azimuth, rho)
        out = solve.(Ref(rotor), sections, op)

        loads_flap = out.Np/1000.
        loads_edge = out.Tp/1000.

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

function calc_moment_stress(mx,my,dx,dy,Rcyl=1.771,tcyl=0.06)

        I = 0.25*np.pi*(Rcyl^4-(Rcyl-tcyl)^4)
        stress = mx*dy/I + my*dx/I

        return stress

end


# #testing during development
r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
              28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
              56.1667, 58.9000, 61.6333]
# chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
#                   3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
# theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
#                   6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
#


x_hub = 0.
y_hub = 0.
z_hub = 90.
yaw = pi/8
azimuth = range(0.,stop=8*pi,length=100)
precone = 0.
tilt = 0*pi/6

fig = figure()
for n = 1:length(azimuth)

        # x_locs, y_locs, z_locs = find_xyz(x_hub,y_hub,z_hub,r,yaw,azimuth[n],precone,tilt)
        x_locs, y_locs, z_locs = find_xyz_simple(x_hub,y_hub,z_hub,r,yaw,azimuth[n])
        clf()

        subplot(131)
        plot(x_locs,y_locs)
        # axis("equal")
        xlim(-70.,70.)
        ylim(-70.,70.)

        subplot(132)
        plot(y_locs,z_locs)
        # axis("equal")
        ylim(20.,160.)
        xlim(-70.,70.)

        subplot(133)
        plot(x_locs,z_locs)
        # axis("equal")
        ylim(20.,160.)
        xlim(-70.,70.)

        sleep(0.05)
        display(fig)

end



# Rhub = 1.5
# Rtip = 63.0
# B = 3
# turbine = true
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
# rotor = Rotor(Rhub, Rtip, B, turbine, pitch, precone)
#
# sections = Section.(r,chord,theta,airfoils)
#
# # Vy = zeros(17)
# # Vx = ones(17).*8.
# # rho = 1.225
# # op = CCBlade.OperatingPoint(Vx, Vy, rho)
#
# yaw = 0.0*pi/180
# tilt = 5.0*pi/180
# hubHt = 90.0
# shearExp = 0.2
#
# Vinf = 10.0
# tsr = 7.55
# rotorR = Rtip*cos(precone)
# Omega = Vinf*tsr/rotorR
# azimuth = 90.0*pi/180
# rho = 1.225
# yaw = 0.
# tilt = 0.
#
# # op = windturbine_op.(Vinf, Omega, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
# velocity_dist = ones(17).*Vinf
# op = distributed_velocity_op.(velocity_dist, Omega, r, precone, yaw, tilt, azimuth, rho)
#
# # out = solve.(Ref(rotor), sections, op)
# azimuth = range(0.,stop=pi*2,length=100)
# M = calculate_root_moment.(Ref(r),Ref(velocity_dist),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,B,rho,precone,hubHt,Omega,pitch,azimuth,yaw,tilt)
# figure()
# plot(azimuth*180/pi, M)
# xlabel("azimuth angle")
# ylabel("root moments (kN-m)")
