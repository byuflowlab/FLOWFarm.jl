using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics
using NPZ

const ff=FlowFarm

include("FAST_data.jl")


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

turb = "high"
ws = 13.0

if turb == "low"
    if ws == 10
        u4 = npzread("flowfields_lowTI/u10_4.npy")
        v4 = npzread("flowfields_lowTI/v10_4.npy")
        w4 = npzread("flowfields_lowTI/w10_4.npy")

        u7 = npzread("flowfields_lowTI/u10_7.npy")
        v7 = npzread("flowfields_lowTI/v10_7.npy")
        w7 = npzread("flowfields_lowTI/w10_7.npy")

        u10 = npzread("flowfields_lowTI/u10_10.npy")
        v10 = npzread("flowfields_lowTI/v10_10.npy")
        w10 = npzread("flowfields_lowTI/w10_10.npy")
    elseif ws == 11
        u4 = npzread("flowfields_lowTI/u11_4.npy").-1.0
        v4 = npzread("flowfields_lowTI/v11_4.npy")
        w4 = npzread("flowfields_lowTI/w11_4.npy")

        u7 = npzread("flowfields_lowTI/u11_7.npy").-1.0
        v7 = npzread("flowfields_lowTI/v11_7.npy")
        w7 = npzread("flowfields_lowTI/w11_7.npy")

        u10 = npzread("flowfields_lowTI/u11_10.npy").-1.0
        v10 = npzread("flowfields_lowTI/v11_10.npy")
        w10 = npzread("flowfields_lowTI/w11_10.npy")
    elseif ws == 12
        u4 = npzread("flowfields_lowTI/u12_4.npy")
        v4 = npzread("flowfields_lowTI/v12_4.npy")
        w4 = npzread("flowfields_lowTI/w12_4.npy")

        u7 = npzread("flowfields_lowTI/u12_7.npy")
        v7 = npzread("flowfields_lowTI/v12_7.npy")
        w7 = npzread("flowfields_lowTI/w12_7.npy")

        u10 = npzread("flowfields_lowTI/u12_10.npy")
        v10 = npzread("flowfields_lowTI/v12_10.npy")
        w10 = npzread("flowfields_lowTI/w12_10.npy")
    elseif ws == 13
        u4 = npzread("flowfields_lowTI/u13_4.npy").+1.0
        v4 = npzread("flowfields_lowTI/v13_4.npy")
        w4 = npzread("flowfields_lowTI/w13_4.npy")

        u7 = npzread("flowfields_lowTI/u13_7.npy").+1.0
        v7 = npzread("flowfields_lowTI/v13_7.npy")
        w7 = npzread("flowfields_lowTI/w13_7.npy")

        u10 = npzread("flowfields_lowTI/u13_10.npy").+1.0
        v10 = npzread("flowfields_lowTI/v13_10.npy")
        w10 = npzread("flowfields_lowTI/w13_10.npy")
    end

elseif turb == "high"
    if ws == 11
        u4 = npzread("flowfields_highTI/u11_4.npy").+1.0
        v4 = npzread("flowfields_highTI/v11_4.npy")
        w4 = npzread("flowfields_highTI/w11_4.npy")

        u7 = npzread("flowfields_highTI/u11_7.npy").+1.0
        v7 = npzread("flowfields_highTI/v11_7.npy")
        w7 = npzread("flowfields_highTI/w11_7.npy")

        u10 = npzread("flowfields_highTI/u11_10.npy").+1.0
        v10 = npzread("flowfields_highTI/v11_10.npy")
        w10 = npzread("flowfields_highTI/w11_10.npy")
    elseif ws == 13
        u4 = npzread("flowfields_lowTI/u13_4.npy")
        v4 = npzread("flowfields_lowTI/v13_4.npy")
        w4 = npzread("flowfields_lowTI/w13_4.npy")

        u7 = npzread("flowfields_lowTI/u13_7.npy")
        v7 = npzread("flowfields_lowTI/v13_7.npy")
        w7 = npzread("flowfields_lowTI/w13_7.npy")

        u10 = npzread("flowfields_lowTI/u13_10.npy")
        v10 = npzread("flowfields_lowTI/v13_10.npy")
        w10 = npzread("flowfields_lowTI/w13_10.npy")
    end

end

x = range(-200.0,stop=200.0,length=161)
y = range(0.0,stop=300.0,length=121)
# XX = zeros(161,121)
# YY = zeros(161,121)
#
# for i=1:161
#    for j=1:121
#        XX[i,j] = x[i]
#        YY[i,j] = y[j]
#    end
# end

# A = range(-200.0,stop=200.0,length=100)
# B = range(0.0,stop=300.0,length=100)
# AA = zeros(length(A),length(B))
# BB = zeros(length(A),length(B))
# UU = zeros(length(A),length(B))
#
# for i=1:length(A)
#    for j=1:length(B)
#        AA[i,j] = A[i]
#        BB[i,j] = B[j]
#        # UU[i,j] = interp2d(akima,x,y,u4,A[i],B[j])[1,1]
#    end
# end

# figure(3)
# pcolormesh(XX,YY,u4,vmin=0.0,vmax=15.0)
# colorbar()
# show()

# low 4 12
# om = [12.08440107, 12.03642098, 11.84670305, 11.42030999, 10.52858835,
#        10.02159635, 10.3127903 , 11.07607054, 11.65258906, 11.90763177,
#        12.00021499]
# low 7 12
# om = [12.07616724, 12.01987292, 11.83136661, 11.50052123, 10.87809062,
#        10.50776205, 10.65150827, 11.12983651, 11.59940669, 11.87620266,
#        11.96680138]
# low 10 12
# om = [12.05448273, 11.99626349, 11.82102287, 11.59350652, 11.23757302,
#        10.98082676, 11.04007275, 11.31570935, 11.62990417, 11.84696304,
#        11.94673847]

#low 4
# om = [11.79803175, 11.53212616, 10.77807941,  9.70458289,  8.75945065,  8.47301996,
#   8.83018058,  9.81286971, 10.80966481, 11.49976043, 11.77493021]

#low 7
# om = [11.72596267, 11.45299696, 10.81810204,  9.93350831,  9.1844484,   8.97723403,
#   9.37077388, 10.18559968, 10.99909754, 11.52076122, 11.74183784]

speeds = range(3.,stop=25.,length=23)
omegas = [6.972,7.183,7.506,7.942,8.469,9.156,10.296,11.431,11.89,
                    12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,
                    12.1,12.1,12.1]
omega_func = Akima(speeds,omegas)
pitches = -1.0.*[0.,0.,0.,0.,0.,0.,0.,0.,0.,3.823,6.602,8.668,10.45,12.055,
                        13.536,14.92,16.226,17.473,18.699,19.941,21.177,22.347,
                        23.469]
pitch_func = Akima(speeds,pitches)

r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
              28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
              56.1667, 58.9000, 61.6333]
chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
                  3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
                  6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]


yaw = 0.
tilt = 5.0*pi/180
rho = 1.225

Rhub = 1.5
Rtip = 63.0
B = 3

precone = 2.5*pi/180

af_path = "/Users/ningrsrch/Dropbox/Projects/waked-loads/5MW_AFFiles_julia"
path1 = af_path*"/Cylinder1.dat"
path2 = af_path*"/Cylinder2.dat"
path3 = af_path*"/DU40_A17.dat"
path4 = af_path*"/DU35_A17.dat"
path5 = af_path*"/DU30_A17.dat"
path6 = af_path*"/DU25_A17.dat"
path7 = af_path*"/DU21_A17.dat"
path8 = af_path*"/NACA64_A17.dat"

af1 = af_from_files(path1)
af2 = af_from_files(path2)
af3 = af_from_files(path3)
af4 = af_from_files(path4)
af5 = af_from_files(path5)
af6 = af_from_files(path6)
af7 = af_from_files(path7)
af8 = af_from_files(path8)

af = [af1,af2,af3,af4,af5,af6,af7,af8]

af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]
airfoils = af[af_idx]


D = 126.4
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

sections = CCBlade.Section.(r,chord,theta,airfoils)

# off = [-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
off = range(-1.5,stop=1.5,length=50)
dams = zeros(length(off))

n = 1.0
mult = 1.


zero = true

sepa = [4.0,7.0,10.0]

#LOW TI
# omega_mult = 1.125

#HIGH TI
omega_mult = 1.12

for k=1:length(sepa)

    sep = sepa[k]

    if sep == 4.0
        u = u4
        v = v4
        w = w4
    elseif sep == 7.0
        u = u7
        v = v7
        w = w7
    elseif sep == 10.0
        u = u10
        v = v10
        w = w10
    end


    rotor_diameter = 126.4
    for i in 1:length(off)

        offset = off[i]*D

        # Omega = om[i]*0.10471975512
        # Omega = Omega*omega_mult
        #
        # println(off[i])
        # println("Omega data: ", om[i])
        # u1 = interp2d(akima,x,y,u,[offset+0.69*rotor_diameter],[90.0])[1,1]
        # u2 = interp2d(akima,x,y,u,[offset-0.69*rotor_diameter],[90.0])[1,1]
        # u3 = interp2d(akima,x,y,u,[offset],[90.0+0.69*rotor_diameter])[1,1]
        # u4 = interp2d(akima,x,y,u,[offset],[90.0-0.69*rotor_diameter])[1,1]
        # ua = (u1+u2+u3+u4)/4.0
        # println("Omega function: ", omega_func(ua))


        u1 = interp2d(akima,x,y,u,[offset+0.69*rotor_diameter],[90.0])[1,1]
        u2 = interp2d(akima,x,y,u,[offset-0.69*rotor_diameter],[90.0])[1,1]
        u3 = interp2d(akima,x,y,u,[offset],[90.0+0.69*rotor_diameter])[1,1]
        u4 = interp2d(akima,x,y,u,[offset],[90.0-0.69*rotor_diameter])[1,1]
        ua = (u1+u2+u3+u4)/4.0
        # ua = interp2d(akima,x,y,u,[offset],[90.0])[1,1]
        Omega = omega_func(ua)*0.10471975512
        Omega = Omega*omega_mult
        pitch = deg2rad(pitch_func(ua))
        rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)


        """0 degrees"""
        az = 0.0
        U = interp2d(akima,x,y,u,zeros(length(r)).+offset,90.0.+r)[1,:]
        # if i == 25
        #     println("0: ", U)
        # end
        V = interp2d(akima,x,y,v,zeros(length(r)).+offset,90.0.+r)[1,:]
        W = interp2d(akima,x,y,w,zeros(length(r)).+offset,90.0.+r)[1,:]
        if zero == true
            V = zeros(length(r))
            W = zeros(length(r))
        end
        V = V.*n
        W = W.*n

        op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)

        out = CCBlade.solve.(Ref(rotor), sections, op)
        M_flap0,M_edge0 = get_moments(out,Rhub,Rtip,r,az)


        """90 degrees"""
        az = pi/2.0
        U = interp2d(akima,x,y,u,offset.-r,90.0.+zeros(length(r)))[:,1]
        # if i == 25
        #     println("90: ", U)
        # end
        V = interp2d(akima,x,y,v,offset.-r,90.0.+zeros(length(r)))[:,1]
        W = interp2d(akima,x,y,w,offset.-r,90.0.+zeros(length(r)))[:,1]
        if zero == true
            V = zeros(length(r))
            W = zeros(length(r))
        end
        V = V.*n
        W = W.*n

        op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)

        out = CCBlade.solve.(Ref(rotor), sections, op)
        M_flap90,M_edge90 = get_moments(out,Rhub,Rtip,r,az)

        """180 degrees"""
        az = pi
        U = interp2d(akima,x,y,u,zeros(length(r)).+offset,90.0.-r)[1,:]
        # if i == 25
        #     println("180: ", U)
        # end
        V = interp2d(akima,x,y,v,zeros(length(r)).+offset,90.0.-r)[1,:]
        W = interp2d(akima,x,y,w,zeros(length(r)).+offset,90.0.-r)[1,:]
        if zero == true
            V = zeros(length(r))
            W = zeros(length(r))
        end
        V = V.*n
        W = W.*n

        op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)

        out = CCBlade.solve.(Ref(rotor), sections, op)
        M_flap180,M_edge180 = get_moments(out,Rhub,Rtip,r,az)

        """270 degrees"""
        az = 3.0*pi/2.0
        U = interp2d(akima,x,y,u,offset.+r,90.0.+zeros(length(r)))[:,1]
        # if i == 25
        #     println("270: ", U)
        # end
        V = interp2d(akima,x,y,v,offset.+r,90.0.+zeros(length(r)))[:,1]
        W = interp2d(akima,x,y,w,offset.+r,90.0.+zeros(length(r)))[:,1]
        if zero == true
            V = zeros(length(r))
            W = zeros(length(r))
        end
        V = V.*n
        W = W.*n

        op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)

        out = CCBlade.solve.(Ref(rotor), sections, op)
        M_flap270,M_edge270 = get_moments(out,Rhub,Rtip,r,az)



        damage = 0.

        su = 70000. # hard coded in now, maybe worth it to let it be passed in
        m = 10.
        years = 25.

        # println(off[i])
        # println(M_flap0)
        # println(M_flap90)
        # println(M_flap180)
        # println(M_flap270)

        for j in 1:Nlocs

            stress0 = ff.calc_moment_stress(M_edge0,M_flap0,xlocs[j],ylocs[j])
            stress90 = ff.calc_moment_stress(M_edge90,M_flap90,xlocs[j],ylocs[j])
            stress180 = ff.calc_moment_stress(M_edge180,M_flap180,xlocs[j],ylocs[j])
            stress270 = ff.calc_moment_stress(M_edge270,M_flap270,xlocs[j],ylocs[j])

            stresses = [stress0,stress90,stress180,stress270]
            cycle = [maximum(stresses),minimum(stresses)]
            mean = sum(cycle)/2.
            alternate = (maximum(cycle)-minimum(cycle))/2.
            effective = alternate/(1.0-mean/su)
            effective = effective*omega_mult

            Nfail = (su/(effective*fos))^m #mLife
            omega_rpm = Omega/0.10471975512
            nCycles = omega_rpm*60.0*24.0*365.25*years
            nCycles = nCycles*mult
            d = nCycles/Nfail
            if d > damage
                    damage = d
            end

        end

        dams[i] = damage
    end

    if k == 1
        subplot(131)
        xlabel("4D")
        ylabel("damage")
    elseif k == 2
        subplot(132)
        xlabel("7D")
        title("high turbulence: 11 m/s")
    else
        subplot(133)
        xlabel("10D")
    end
    plot(off,dams)
    # title("SOWFA")
    ylim([0,1.75])
    println("damage: ", dams)

    FS,FD = fastdata(turb,ws,sep)
    scatter(FS,FD)

    # legend(["mean SOWFA velocity","SOWFA+FAST"])
end

tight_layout()
show()
println("finished")
