using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics
using NPZ

const ff=FlowFarm


# function add_turbulence(mean,TI,)
#
# end

function rainflow(array_ext;uc_mult=0.5)

    # """ Rainflow counting of a signal's turning points with Goodman correction
    #
    #     Args:
    #         array_ext (numpy.ndarray): array of turning points
    #
    #     Keyword Args:
    #         uc_mult (float): partial-load scaling [opt, default=0.5]
    #
    #     Returns:
    #         array_out (numpy.ndarray): (3 x n_cycle) array of rainflow values:
    #                                     1) load range
    #                                     2) range mean
    #                                     3) cycle count
    # """

    tot_num = length(array_ext)             # total size of input array
    array_out = zeros(3, tot_num-1)         # initialize output array

    pr = 1                                  # index of input array
    po = 1                                  # index of output array
    j = 0                                   # index of temporary array "a"
    a  = zeros(tot_num)                     # temporary array for algorithm

    # loop through each turning point stored in input array
    for i = 1:tot_num
        j += 1                  # increment "a" counter
        a[j] = array_ext[pr]    # put turning point into temporary array
        pr += 1                 # increment input array pointer

        while j >= 3 && abs( a[j-1] - a[j-2]) <= abs(a[j] - a[j-1])

            lrange = abs( a[j-1] - a[j-2] )

            # partial range
            if j == 3
                mean      = ( a[1] + a[2] ) / 2.0
                a[1]=a[2]
                a[2]=a[3]
                j=2
                if lrange > 0
                    array_out[1,po] = lrange
                    array_out[2,po] = mean
                    array_out[3,po] = uc_mult
                    po += 1
                end

            # full range
            else
                mean      = ( a[j-1] + a[j-2] ) / 2.0
                a[j-2]=a[j]
                j=j-2
                if (lrange > 0)
                    array_out[1,po] = lrange
                    array_out[2,po] = mean
                    array_out[3,po] = 1.00
                    po += 1
                end
            end
        end
    end

    # partial range
    for i = 1:j-1
        lrange    = abs( a[i] - a[i+1] )
        mean      = ( a[i] + a[i+1] ) / 2.0
        if lrange > 0
            array_out[1,po] = lrange
            array_out[2,po] = mean
            array_out[3,po] = uc_mult
            po += 1
        end
    end

    # get rid of unused entries
    out = array_out[:,1:po-1]

    return out
end


function get_peaks(array)
    A = array[:]
    # get rid of any zero slope in the beginning
    while A[2] == A[1]
        A = A[2:length(A)]
    end
    peaks = [A[1]]
    if A[2] > A[1]
        slope = "p"
    elseif A[2] < A[1]
        slope = "m"
    end
    for i = 1:length(A)-2
        ind = i+1
        if slope == "p"
            if A[ind+1] < A[ind]
                peaks = append!(peaks,A[ind])
                slope = "m"
            end
        elseif slope == "m"
            if A[ind+1] > A[ind]
                peaks = append!(peaks,A[ind])
                slope = "p"
            end
        end
    end
    peaks = append!(peaks,A[length(A)])
    return peaks
end


function get_peaks_indices(array)
    A = array[:]
    # get rid of any zero slope in the beginning
    while A[2] == A[1]
        A = A[2:length(A)]
    end
    peaks = [1]
    if A[2] > A[1]
        slope = "p"
    elseif A[2] < A[1]
        slope = "m"
    end
    for i = 1:length(A)-2
        ind = i+1
        if slope == "p"
            if A[ind+1] < A[ind]
                peaks = append!(peaks,ind)
                slope = "m"
            end
        elseif slope == "m"
            if A[ind+1] > A[ind]
                peaks = append!(peaks,ind)
                slope = "p"
            end
        end
    end
    peaks = append!(peaks,length(A))
    return peaks
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
#
turb = "low"
ws = 11.0

include("model.jl")

TI_free = 0.046
TI = "low"


ka = 0.38
kb = 0.004

initial = 0.313
constant = 1.931
ai = 0.435
downstream = -0.855
alpha_star = 2.32
beta_star = 0.154
turb_index = 2

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

# off = [-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
# off = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0]
off = range(-1.5,stop=1.5,length=20)
# off = [-0.4]
dams = zeros(length(off))


zero = true

sep = 4.0
turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)


turbine_x = [0.0,sep*rotor_diameter]
init_inflow_velcities = zeros(length(turbine_x)).+ws

omega_func = Akima(speeds, omegas)
pitch_func = Akima(speeds, pitches)

tilt = deg2rad(5.0)
rho = 1.225
vw = 0.0

TI_check = 0.1

TIdat = [0.04895652, 0.05049786, 0.04961281, 0.04932443, 0.04677447,
       0.04837576, 0.0477134 , 0.04751765, 0.04458461, 0.04660241,
       0.0453109 , 0.04441194, 0.0406133 , 0.04219035, 0.04144826,
       0.04123373, 0.03863729, 0.04068373, 0.04062233, 0.0410512 ,
       0.03957269, 0.04177059, 0.0411834 , 0.04119085, 0.03889296,
       0.04037401, 0.04017531, 0.04073167, 0.03930292, 0.04286042,
       0.04310193, 0.04375134, 0.04129519, 0.04432437, 0.04490191,
       0.0459663 , 0.04482169, 0.0492744 , 0.05076854, 0.05278087,
       0.05167714, 0.05642463, 0.05785846, 0.05997315, 0.0602225 ,
       0.06707714, 0.07219531, 0.07850824, 0.08500409, 0.10277199,
       0.1117703 , 0.12174509, 0.12251037, 0.13973445, 0.14484071,
       0.15102847, 0.1425917 , 0.14825097, 0.14962787, 0.15270284,
       0.14898782, 0.15669054, 0.15610436, 0.15769035, 0.15041006,
       0.15285879, 0.15239825, 0.15375462, 0.15034   , 0.15399772,
       0.15354361, 0.15408179, 0.15109033, 0.15266166, 0.15162482,
       0.1512091 , 0.148946  , 0.14712916, 0.1462822 , 0.14624127,
       0.14813103, 0.15139834, 0.15438181, 0.15790578, 0.15992832,
       0.1640869 , 0.16770551, 0.17188595, 0.17625248, 0.18590122,
       0.18781925, 0.19070745, 0.18601297, 0.19079224, 0.18929965,
       0.18905714, 0.18220195, 0.18635539, 0.18231762, 0.18043071,
       0.17050113, 0.1732225 , 0.16786195, 0.16468851, 0.15057797,
       0.14780192, 0.13916778, 0.13267622, 0.1188933 , 0.11544137,
       0.10749484, 0.10158743, 0.08959862, 0.08671188, 0.07978171,
       0.07428687, 0.06475171, 0.0627464 , 0.05871523, 0.05611433,
       0.05166367, 0.05387342, 0.05308169, 0.05322432, 0.05020915,
       0.05370841, 0.05376202, 0.05451429, 0.05135909, 0.05394679,
       0.05357917, 0.05365005, 0.05086475, 0.05244815, 0.05186676,
       0.05170529, 0.04959481, 0.05203156, 0.05204569, 0.05261812,
       0.05114477, 0.05440483, 0.05455887, 0.05520471, 0.05275623,
       0.05472298, 0.05393864, 0.05377447, 0.0506379 , 0.0525567 ,
       0.05198626, 0.05212971, 0.0497818 , 0.05223933, 0.05196204,
       0.05234997, 0.05112163, 0.05376753, 0.05415094, 0.05499336,
       0.05394576]

TIlocs = range(-200.0, stop=200.0, length=length(TIdat))
TIfunc = Akima(TIlocs, TIdat)

# plot(TIlocs,TIdat)


nCycles = 200
# dist = 25.5
# dist = 15.0
dist = 30.0
# dist = 40.0
naz = 8
az_arr = range(0.0,stop=2.0*pi-2.0*pi/naz,length=naz)
for k=1:length(off)
    println("offset: ", off[k])
    offset = off[k]*rotor_diameter
    turbine_y = [0.0,offset]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, init_inflow_velcities, zeros(nturbines))
    # TI_check = TIfunc(offset)
    # TI_check = 0.09
    flap = []
    edge = []
    oms = []

    for i = 1:nCycles

        for j = 1:naz

            az = az_arr[j]
            # TI_check = TIfunc(offset-sin(az)*dist)
            # if offset-dist < -200.0
            #     TI_check = TIdat[1]
            # elseif offset+dist > 200.0
            #     TI_check = TIdat[end]
            # end

            TI_check_arr = zeros(length(r))
            for k = 1:length(r)
                Ly = offset-sin(az)*r[k]
                Lz = cos(az)*r[k]
                R = sqrt(Ly^2+Lz^2)
                if R < 70.
                    TI_check_arr[k] = 0.16
                else
                    TI_check_arr[k] = 0.05
                end
            end

            TI_check = sum(TI_check_arr)/length(TI_check_arr)
            # TI_check = sum(TI_check_arr .* r)/(length(TI_check_arr)*22.5)
            windspeeds = [ws] + randn(1).*TI_check.*ws
            turbine_inflow_velcities = zeros(nturbines) .+ windspeeds
            windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
            pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])

            ff.turbine_velocities_one_direction!(points_x, points_y, ms, pd)
            turbine_inflow_velcities = pd.wind_farm_states[1].turbine_inflow_velcities
            Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
            Omega = Omega_rpm*0.10471975512 #convert to rad/s
            pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
            pitch = pitch_deg*pi/180.0

            V = -ones(length(r)).*(cos(az)*vw)
            W = -ones(length(r)).*(sin(az)*vw)
            U = ff.get_speeds(turbine_x,turbine_y,turb_index,hubHt,r,yaw,az,ms,pd)
            op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)
            out = CCBlade.solve.(Ref(rotor), sections, op)
            flap_temp,edge_temp = get_moments(out,Rhub,Rtip,r,az)
            flap = append!(flap,flap_temp)
            edge = append!(edge,edge_temp)
            oms = append!(oms,Omega)
        end

        # az = deg2rad(0.0)
        # TI_check = TIfunc(offset)
        # windspeeds = [ws] + randn(1).*TI_check.*ws
        # turbine_inflow_velcities = zeros(nturbines) .+ windspeeds
        # windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
        # pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
        #
        # ff.turbine_velocities_one_direction!(points_x, points_y, ms, pd)
        # turbine_inflow_velcities = pd.wind_farm_states[1].turbine_inflow_velcities
        # Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
        # Omega = Omega_rpm*0.10471975512 #convert to rad/s
        # pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
        # pitch = pitch_deg*pi/180.0
        #
        # V = -ones(length(r)).*(cos(az)*vw)
        # W = -ones(length(r)).*(sin(az)*vw)
        # U = ff.get_speeds(turbine_x,turbine_y,turb_index,hubHt,r,yaw,az,ms,pd)
        # op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)
        # out = CCBlade.solve.(Ref(rotor), sections, op)
        # flap_temp,edge_temp = get_moments(out,Rhub,Rtip,r,az)
        # flap = append!(flap,flap_temp)
        # edge = append!(edge,edge_temp)
        # oms = append!(oms,Omega)
        #
        #
        # az = deg2rad(90.0)
        # TI_check = TIfunc(offset-dist)
        # if offset-dist < -200.0
        #     TI_check = TIdat[1]
        # end
        # windspeeds = [ws] + randn(1).*TI_check.*ws
        # turbine_inflow_velcities = zeros(nturbines) .+ windspeeds
        # windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
        # pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
        #
        # ff.turbine_velocities_one_direction!(points_x, points_y, ms, pd)
        # turbine_inflow_velcities = pd.wind_farm_states[1].turbine_inflow_velcities
        # Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
        # Omega = Omega_rpm*0.10471975512 #convert to rad/s
        # pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
        # pitch = pitch_deg*pi/180.0
        #
        # V = -ones(length(r)).*(cos(az)*vw)
        # W = -ones(length(r)).*(sin(az)*vw)
        # U = ff.get_speeds(turbine_x,turbine_y,turb_index,hubHt,r,yaw,az,ms,pd)
        # op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)
        # out = CCBlade.solve.(Ref(rotor), sections, op)
        # flap_temp,edge_temp = get_moments(out,Rhub,Rtip,r,az)
        # flap = append!(flap,flap_temp)
        # edge = append!(edge,edge_temp)
        # oms = append!(oms,Omega)
        #
        #
        # az = deg2rad(180.0)
        # TI_check = TIfunc(offset)
        # windspeeds = [ws] + randn(1).*TI_check.*ws
        # turbine_inflow_velcities = zeros(nturbines) .+ windspeeds
        # windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
        # pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
        #
        # ff.turbine_velocities_one_direction!(points_x, points_y, ms, pd)
        # turbine_inflow_velcities = pd.wind_farm_states[1].turbine_inflow_velcities
        # Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
        # Omega = Omega_rpm*0.10471975512 #convert to rad/s
        # pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
        # pitch = pitch_deg*pi/180.0
        #
        # V = -ones(length(r)).*(cos(az)*vw)
        # W = -ones(length(r)).*(sin(az)*vw)
        # U = ff.get_speeds(turbine_x,turbine_y,turb_index,hubHt,r,yaw,az,ms,pd)
        # op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)
        # out = CCBlade.solve.(Ref(rotor), sections, op)
        # flap_temp,edge_temp = get_moments(out,Rhub,Rtip,r,az)
        # flap = append!(flap,flap_temp)
        # edge = append!(edge,edge_temp)
        # oms = append!(oms,Omega)
        #
        #
        # az = deg2rad(270.0)
        # TI_check = TIfunc(offset+dist)
        # if offset+dist > 200.
        #     TI_check = TIdat[end]
        # end
        # windspeeds = [ws] + randn(1).*TI_check.*ws
        # turbine_inflow_velcities = zeros(nturbines) .+ windspeeds
        # windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])
        # pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
        #
        # ff.turbine_velocities_one_direction!(points_x, points_y, ms, pd)
        # turbine_inflow_velcities = pd.wind_farm_states[1].turbine_inflow_velcities
        # Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
        # Omega = Omega_rpm*0.10471975512 #convert to rad/s
        # pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
        # pitch = pitch_deg*pi/180.0
        #
        # V = -ones(length(r)).*(cos(az)*vw)
        # W = -ones(length(r)).*(sin(az)*vw)
        # U = ff.get_speeds(turbine_x,turbine_y,turb_index,hubHt,r,yaw,az,ms,pd)
        # op = multiple_components_op.(U, V, W, Omega, r, precone, yaw, tilt, az, rho)
        # out = CCBlade.solve.(Ref(rotor), sections, op)
        # flap_temp,edge_temp = get_moments(out,Rhub,Rtip,r,az)
        # # if i == 1
        # #     println(out.Np)
        # #     println(flap_temp,edge_temp)
        # # end
        # flap = append!(flap,flap_temp)
        # edge = append!(edge,edge_temp)
        # oms = append!(oms,Omega)
    end

    Nlocs = 20
    xlocs = zeros(Nlocs)
    ylocs = zeros(Nlocs)
    angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
    root_rad=3.542/2.
    for i in 1:Nlocs
        xlocs[i] = cos(angles[i])*root_rad
        ylocs[i] = sin(angles[i])*root_rad
    end

    # xlocs = [cos(pi/4.0)*root_rad]
    # ylocs = [sin(pi/4.0)*root_rad]

    avg_omega = sum(oms)/length(oms)
    total_time = nCycles/(avg_omega/(2.0*pi))

    global damage
    damage = 0.0
    for i  = 1:Nlocs
        global damage
        sigma = ff.calc_moment_stress.(edge,flap,xlocs[i],ylocs[i])
        peaks = get_peaks(sigma)
        out = rainflow(peaks)

        alternate = out[1,:]/2.
        mean = out[2,:]
        count = out[3,:]

        su = 70000.
        m = 10.
        years = 25.
        freq = 1.0

        mar = alternate./(1.0.-mean./su)

        # figure(1)
        # title("mar, julia")
        # hist(mar,bins=50,range=[6000,11000])

        npts = length(mar)
        #damage calculations
        d = 0.0
        for i = 1:npts
            Nfail = ((su)/(mar[i]*fos))^m
            mult = years*365.25*24.0*3600.0*freq/total_time
            d += count[i]*mult/Nfail
        end
        if d > damage
            damage = d
        end
    end
    dams[k] = damage
end
println("damage: ", dams)
plot(off,dams)

FS,FD = fastdata(turb,ws,sep)
scatter(FS,FD)
