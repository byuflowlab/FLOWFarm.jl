using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics

const ff=FlowFarm

function calc_TI(constant,ai,TI_free,initial,sep,downstream)

    ti_calculation = constant * (1.0/3.0)^ai * TI_free^initial * sep^downstream

    TI = sqrt(ti_calculation^2 + TI_free^2)

    return TI
end

function calc_c(TI,ws)
    alpha = 22.
    beta = 0.01
    c = 1.
    return TI*alpha*(ws/11.4)
    # return 1.1
end


function turbulence_damage(wind_speed,turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180,Nlocs=20,root_rad=3.542/2.,fos=1.15,mult=1.0,dt=0.025)

    omega_func = Akima(speeds, omegas)
    Omega_rpm = omega_func(wind_speed) #in rpm
    Omega = Omega_rpm*0.10471975512 #convert to rad/s

    pitch_func = Akima(speeds, pitches)
    pitch_deg = pitch_func(wind_speed) #in degrees
    pitch = pitch_deg*pi/180.0


    damage = 0.

    m1 = 1000.0
    m2 = 20.0
    m3 = 1.0
    m1 = 0
    m2 = 0
    m3 = 0

    tm = 1.
    high_speed = wind_speed*(1+tm*turbulence_intensity)
    low_speed = wind_speed*(1-tm*turbulence_intensity)

    moments_high = ff.calculate_root_moment.(Ref(r),Ref(ones(length(r))*high_speed),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,0.0,precone=precone,tilt=tilt)
    moments_low = ff.calculate_root_moment.(Ref(r),Ref(ones(length(r))*low_speed),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,0.0,precone=precone,tilt=tilt)

    angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
    xlocs = zeros(Nlocs)
    ylocs = zeros(Nlocs)
    for i in 1:Nlocs
        xlocs[i] = cos(angles[i])*root_rad
        ylocs[i] = sin(angles[i])*root_rad
    end

    dp = 0.

    su = 70000. # hard coded in now, maybe worth it to let it be passed in
    m = 10.
    years = 25.

    for i in 1:Nlocs

        stress_high = ff.calc_moment_stress(moments_high[2],moments_high[1],xlocs[i],ylocs[i])
        stress_low = ff.calc_moment_stress(moments_low[2],moments_low[1],xlocs[i],ylocs[i])

        cycle = [stress_high,stress_low]
        mean = (stress_high+stress_low)/2.
        alternate = (stress_high-stress_low)/2.
        effective = alternate/(1.0-mean/su)

        Nfail = (su/(effective*fos))^m #mLife
        # nCycles = Omega_rpm*60.0*24.0*365.25*years
        # nCycles = 3600.0/dt*24.0*365.25*years
        nCycles = m1*24.0*365.25*years

        nCycles = nCycles*mult
        d = nCycles/Nfail
        if d > dp
                dp = d
                # println(damage)
        end

    end

    damage += dp


    tm = 2.
    high_speed = wind_speed*(1+tm*turbulence_intensity)
    low_speed = wind_speed*(1-tm*turbulence_intensity)

    moments_high = ff.calculate_root_moment.(Ref(r),Ref(ones(length(r))*high_speed),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,0.0,precone=precone,tilt=tilt)
    moments_low = ff.calculate_root_moment.(Ref(r),Ref(ones(length(r))*low_speed),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,0.0,precone=precone,tilt=tilt)

    angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
    xlocs = zeros(Nlocs)
    ylocs = zeros(Nlocs)
    for i in 1:Nlocs
        xlocs[i] = cos(angles[i])*root_rad
        ylocs[i] = sin(angles[i])*root_rad
    end

    dp = 0.

    su = 70000. # hard coded in now, maybe worth it to let it be passed in
    m = 10.
    years = 25.

    for i in 1:Nlocs

        stress_high = ff.calc_moment_stress(moments_high[2],moments_high[1],xlocs[i],ylocs[i])
        stress_low = ff.calc_moment_stress(moments_low[2],moments_low[1],xlocs[i],ylocs[i])

        cycle = [stress_high,stress_low]
        mean = (stress_high+stress_low)/2.
        alternate = (stress_high-stress_low)/2.
        effective = alternate/(1.0-mean/su)

        Nfail = (su/(effective*fos))^m #mLife
        # nCycles = Omega_rpm*60.0*24.0*365.25*years
        # nCycles = 3600.0/dt*24.0*365.25*years
        nCycles = m2*24.0*365.25*years

        nCycles = nCycles*mult
        d = nCycles/Nfail
        if d > dp
                dp = d
                # println(damage)
        end

    end

    damage += dp


    tm = 3.
    high_speed = wind_speed*(1+tm*turbulence_intensity)
    low_speed = wind_speed*(1-tm*turbulence_intensity)

    moments_high = ff.calculate_root_moment.(Ref(r),Ref(ones(length(r))*high_speed),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,0.0,precone=precone,tilt=tilt)
    moments_low = ff.calculate_root_moment.(Ref(r),Ref(ones(length(r))*low_speed),Ref(chord),Ref(theta),Ref(airfoils),Rhub,Rtip,hubHt,Omega,pitch,0.0,precone=precone,tilt=tilt)

    angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
    xlocs = zeros(Nlocs)
    ylocs = zeros(Nlocs)
    for i in 1:Nlocs
        xlocs[i] = cos(angles[i])*root_rad
        ylocs[i] = sin(angles[i])*root_rad
    end

    dp = 0.

    su = 70000. # hard coded in now, maybe worth it to let it be passed in
    m = 10.
    years = 25.

    for i in 1:Nlocs

        stress_high = ff.calc_moment_stress(moments_high[2],moments_high[1],xlocs[i],ylocs[i])
        stress_low = ff.calc_moment_stress(moments_low[2],moments_low[1],xlocs[i],ylocs[i])

        cycle = [stress_high,stress_low]
        mean = (stress_high+stress_low)/2.
        alternate = (stress_high-stress_low)/2.
        effective = alternate/(1.0-mean/su)

        Nfail = (su/(effective*fos))^m #mLife
        # nCycles = Omega_rpm*60.0*24.0*365.25*years
        # nCycles = 3600.0/dt*24.0*365.25*years
        nCycles = m3*24.0*365.25*years

        nCycles = nCycles*mult
        d = nCycles/Nfail
        if d > dp
                dp = d
                # println(damage)
        end

    end

    damage += dp

    return damage
end

include("model.jl")

TI_free = 0.046
TI = "low"
wind_speed = 12.
windspeeds = [wind_speed]
turbine_inflow_velcities = zeros(nturbines) .+ wind_speed
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])


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

# sweep = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
sweep = range(-5,stop=5,length=50)
damage = zeros(length(sweep))
turb_index = 2

# mult = 1.12 # (for 12)
c = 1.15
b = 0.07
a = (TI_free/0.046)^b
rated_power = 11.4

points_x = [0.69,0,-0.69,0]
points_y = [0,0.69,0,-0.69]

sep = 4.0
# turbulence_intensity = 0.12
# turbulence_intensity = 0.14699300659238304
turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# mult = calc_c(TI_free,wind_speed)
mult = 1.0
m2 = 2.0


ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
for i = 1:length(sweep)
    turbine_y = [0.0,sweep[i]*rotor_diameter]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
    pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    waked_damage = ff.get_damage(turbine_x,turbine_y,turb_index,speeds,omegas,pitches,wind_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,ms2,pd2,wind_speed,Nlocs=20,shearExp=shearexponent,mult=mult,m2=m2)
    ff.turbine_velocities_one_direction!(points_x, points_y, ms2, pd2)
    in_speeds = pd2.wind_farm_states[1].turbine_inflow_velcities
    turbulent_damage = turbulence_damage(in_speeds[turb_index],turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180)
    damage[i] = waked_damage + turbulent_damage
end
# subplot(131)
# plot(sweep,damage)
# ylim([0,3])
# title("4")

sep = 7.0

# turbulence_intensity = 0.07
# turbulence_intensity = 0.0979887388404962
turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# mult = calc_c(TI_free,wind_speed)
mult = 1.0

ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
for i = 1:length(sweep)
    turbine_y = [0.0,sweep[i]*rotor_diameter]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
    pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    waked_damage = ff.get_damage(turbine_x,turbine_y,turb_index,speeds,omegas,pitches,wind_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,ms2,pd2,wind_speed,Nlocs=20,shearExp=shearexponent,mult=mult,m2=m2)
    ff.turbine_velocities_one_direction!(points_x, points_y, ms2, pd2)
    in_speeds = pd2.wind_farm_states[1].turbine_inflow_velcities
    turbulent_damage = turbulence_damage(in_speeds[turb_index],turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180)
    damage[i] = waked_damage + turbulent_damage
end
# subplot(132)
# plot(sweep,damage)
# ylim([0,3])
# title("7")

sep = 10.0
# turbulence_intensity = 0.05
# turbulence_intensity = 0.07863689518866869
turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# mult = calc_c(TI_free,wind_speed)
mult = 1.0

ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
for i = 1:length(sweep)
    turbine_y = [0.0,sweep[i]*rotor_diameter]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
    pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    waked_damage = ff.get_damage(turbine_x,turbine_y,turb_index,speeds,omegas,pitches,wind_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,ms2,pd2,wind_speed,Nlocs=20,shearExp=shearexponent,mult=mult,m2=m2)
    ff.turbine_velocities_one_direction!(points_x, points_y, ms2, pd2)
    in_speeds = pd2.wind_farm_states[1].turbine_inflow_velcities
    turbulent_damage = turbulence_damage(in_speeds[turb_index],turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180)
    damage[i] = waked_damage + turbulent_damage
end
# subplot(133)
# plot(sweep,damage)
# ylim([0,3])
# title("10")
# println(damage)


FS = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]

if TI == "high"
    if wind_speed == 11.0

        FD4 = [0.2728365 , 0.65977472, 0.77260459, 0.4778191 , 0.24800923,
       0.10146671, 0.04171074, 0.06843597, 0.06012162, 0.05531389,
       0.15369841]

        FD7 = [0.30243331, 0.37563574, 0.38496825, 0.30413719, 0.17994517,
       0.08185298, 0.04101291, 0.04436714, 0.04838322, 0.09869199,
       0.17469088]

        FD10 = [0.25667372, 0.28045112, 0.31380656, 0.26435168, 0.13541739,
       0.08130023, 0.06937866, 0.05838691, 0.06783909, 0.09664072,
       0.19967812]

    elseif wind_speed == 12.0

        """high TI"""
        FD4 = [0.59436966, 0.86531902, 1.47720059, 1.29034929, 0.39374063,
               0.0675834 , 0.09811606, 0.0944617 , 0.0546642 , 0.1133815 ,
               0.30375038]
        FD7 = [0.74362544, 0.94617913, 1.12020256, 0.82845987, 0.41520515,
               0.09437026, 0.08533575, 0.09938972, 0.13268212, 0.24422962,
               0.28269637]
        FD10 = [0.8467823 , 0.76667692, 0.7433349 , 0.69224899, 0.38166549,
               0.17593557, 0.12199661, 0.12614291, 0.12928943, 0.19476268,
               0.26522156]

    elseif wind_speed == 13.0

        FD4 = [0.48608728, 0.90662317, 1.39549172, 2.04401001, 1.46901554,
       0.54139495, 0.14399162, 0.1125801 , 0.14678671, 0.1724072 ,
       0.33570147]

        FD7 = [0.84783926, 0.90382795, 1.66442715, 2.02387409, 1.79021644,
       0.47391596, 0.27280782, 0.17559805, 0.20302434, 0.21569954,
       0.44578776]

        FD10 = [0.91365132, 1.29885624, 1.38659419, 1.02947128, 0.91010369,
       0.41046279, 0.24486986, 0.20180934, 0.24375777, 0.26666614,
       0.6751095 ]

   elseif wind_speed == 8.0
       FD4 = [0.0413936 , 0.06184613, 0.0739998 , 0.0547903 , 0.03016046,
       0.01753587, 0.01669262, 0.01942399, 0.01987035, 0.02022905,
       0.02742349]
       FD7 = [0.03588626, 0.04584446, 0.04958232, 0.04257099, 0.02926373,
       0.01875865, 0.0168872 , 0.01832914, 0.01869439, 0.02128756,
       0.02407287]
       FD10 = [0.03309448, 0.03835529, 0.03819291, 0.03135918, 0.02562621,
       0.02056727, 0.01861202, 0.01865587, 0.02037396, 0.02085615,
       0.02464827]

    end

elseif TI == "low"
    if wind_speed == 11.0

        FD4 = [0.12474905, 0.44955552, 0.98879786, 0.55243793, 0.13710362,
                0.02830596, 0.02459855, 0.04351297, 0.0610912 , 0.04793649,
                0.04335935]
        FD7 = [0.15679473, 0.33452526, 0.5630296 , 0.43909217, 0.1547891 ,
                0.04392345, 0.02885432, 0.04226767, 0.04382467, 0.0332834 ,
                0.0453312 ]
        FD10 = [0.13796973, 0.24582433, 0.3031502 , 0.26217174, 0.12378311,
                0.04605147, 0.02505028, 0.02983001, 0.03059458, 0.03172991,
                0.04559384]

    elseif wind_speed == 12.0
        """low TI"""
        FD4 = [0.14841399, 0.49595138, 1.34803516, 1.25430712, 0.29605346,
                       0.05068178, 0.02653278, 0.04368748, 0.05165779, 0.03799159,
                       0.07522253]
        FD7 = [0.18327844, 0.41799378, 0.76867414, 0.86818366, 0.35801523,
                       0.0841998 , 0.02661469, 0.03694555, 0.04142776, 0.03735182,
                       0.07156566]
        FD10 = [0.2084088 , 0.3576736 , 0.48538222, 0.50291998, 0.29240833,
                       0.1060008 , 0.04258391, 0.02906886, 0.03303259, 0.04580801,
                       0.07174415]

    elseif wind_speed == 13.0
        FD4 = [0.120385  , 0.43966496, 1.49963083, 1.283816  , 0.34256675,
                   0.05646229, 0.03843661, 0.05306928, 0.05401545, 0.08259596,
                   0.08577113]
        FD7 = [0.17079266, 0.52412727, 1.10421169, 0.91336942, 0.32727372,
                    0.07652119, 0.04013585, 0.04518639, 0.06678548, 0.06778572,
                    0.09812841]
        FD10 = [0.20879636, 0.4840231 , 0.88389687, 0.74301718, 0.3106862 ,
                    0.09396246, 0.04306845, 0.0497348 , 0.07675588, 0.07744574,
                    0.10162294]
    elseif wind_speed == 10.0
        FD4 = [0.09182389, 0.27367474, 0.43737724, 0.25577379, 0.07315111,
       0.02212047, 0.0219701 , 0.03314924, 0.04282983, 0.03492318,
       0.03608086]
        FD7 = [0.11118407, 0.20795058, 0.27237578, 0.19238532, 0.07942963,
       0.03005412, 0.02408143, 0.02873857, 0.03003227, 0.02703879,
       0.03532939]
        FD10 = [0.09962295, 0.14386566, 0.1529574 , 0.10837712, 0.06500144,
       0.03931054, 0.02685782, 0.02341648, 0.02369127, 0.02745877,
       0.03709964]
    end
end
# subplot(131)
# scatter(FS,FD4)
# xlabel("offset (D)")
# ylabel("damage")
#
# subplot(132)
# scatter(FS,FD7)
# xlabel("offset (D)")
#
# subplot(133)
# scatter(FS,FD10)
# xlabel("offset (D)")
#
# suptitle("10 m/s, low turbulence")
# tight_layout()
#
# show()
