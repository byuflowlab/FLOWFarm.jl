using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics

const ff=FlowFarm

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


# # #testing during development
r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
              28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
              56.1667, 58.9000, 61.6333]
chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
                  3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
                  6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

Rhub = 1.5
Rtip = 63.0
B = 3
pitch = 0.0
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

pitch = 0.
Rhub = 1.5
hubHt = 90.

speeds = range(3.,stop=25.,length=23)
omegas = [6.972,7.183,7.506,7.942,8.469,9.156,10.296,11.431,11.89,
                    12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,
                    12.1,12.1,12.1]
pitches = -1.0.*[0.,0.,0.,0.,0.,0.,0.,0.,0.,3.823,6.602,8.668,10.45,12.055,
                        13.536,14.92,16.226,17.473,18.699,19.941,21.177,22.347,
                        23.469]


rotor_diameter = 126.4
hub_height = 90.0
yaw = 0.0
ct = 8.0/9.0
cp = 0.42

generator_efficiency = 0.944
ai = 1.0/3.0
wind_speed = 12.
air_density = 1.225  # kg/m^3
nturbines = 2
# turbine_y = zeros(nturbines)
turbine_y = [0.,-rotor_diameter/2.]
turbine_z = zeros(nturbines)
turbine_yaw = zeros(nturbines)
turbine_ct = zeros(nturbines) .+ ct
turbine_ai = zeros(nturbines) .+ ai
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheight = [hub_height]
shearexponent = 0.1
turbine_inflow_velcities = zeros(nturbines) .+ wind_speed

# rotor sample points
rotor_points_y = [0.0]
rotor_points_z = [0.0]

ct_model = ff.ConstantCt(ct)
power_model = ff.ConstantCp([cp], [generator_efficiency])
wind_shear_model = ff.PowerLawWindShear(shearexponent)

turbine = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [ct_model], [power_model])
# turbine_definitions = [turbine for i in 1:nturbines]
turbine_definitions = [turbine]
sorted_turbine_index = [i for i  in 1:nturbines]
turbine_definition_ids = ones(Int, nturbines)

windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])

# turbulence_intensity = 0.046
# turbulence_intensity = 0.147
# turbulence_intensity = 0.157
ka = 0.38
kb = 0.004
# ka = 0.38424141992410843
# kb = 0.01850438686814998
# ka = 0.1
# kb = 0.005


# ky = ka*turbulence_intensity + kb
# kz = ka*turbulence_intensity + kb
# ky = 0.0284
# kz = 0.0284
# ky = 0.041
# kz = 0.041
# horizontal_spread_rate = ky
# vertical_spread_rate = kz
# alpha_star = 2.32
# beta_star = 0.154

alpha_star = 2.
beta_star = 0.1
# wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
#
# wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()

# ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

# sweep = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
sweep = range(-2.0,stop=2.0,length=100)
# sweep = [-0.5]
damage = zeros(length(sweep))
turb_index = 2

# figure(figsize=(9,4))

mult = 1.0
points_x = [0.69,0,-0.69,0]
points_y = [0,0.69,0,-0.69]

sep = 4.0
# turbulence_intensity = 0.14699300659238304
# turbulence_intensity = 0.2138695361889833
turbulence_intensity = 0.12
# turbulence_intensity = 0.1
# turbulence_intensity = 0.07
# turbulence_intensity = 0.1338627218981348
# turbulence_intensity = 0.1

ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
for i = 1:length(sweep)
    # println("sweep[i]: ", sweep[i])
    turbine_y = [0.0,sweep[i]*rotor_diameter]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
    pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    waked_damage = ff.get_damage(turbine_x,turbine_y,turb_index,speeds,omegas,pitches,wind_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,ms2,pd2,wind_speed,Nlocs=20,shearExp=shearexponent,mult=mult)
    # println("damage loads: ", damage[i])

    ff.turbine_velocities_one_direction!(points_x, points_y, ms2, pd2)
    in_speeds = pd2.wind_farm_states[1].turbine_inflow_velcities
    turbulent_damage = turbulence_damage(in_speeds[turb_index],turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180)
    # println("damage turbulence: ", turbulent_damage)
    damage[i] = waked_damage + turbulent_damage
end
subplot(131)
plot(sweep,damage)
ylim([0,1.5])
title("4")

sep = 7.0
# turbulence_intensity = 0.0979887388404962
# turbulence_intensity = 0.15817718313411153
turbulence_intensity = 0.07
# turbulence_intensity = 0.055
# turbulence_intensity = 0.04
# turbulence_intensity = 0.07563260011618758

ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
for i = 1:length(sweep)
    # println("sweep[i]: ", sweep[i])
    turbine_y = [0.0,sweep[i]*rotor_diameter]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
    pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    waked_damage = ff.get_damage(turbine_x,turbine_y,turb_index,speeds,omegas,pitches,wind_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,ms2,pd2,wind_speed,Nlocs=20,shearExp=shearexponent,mult=mult)
    # println("damage loads: ", damage[i])
    ff.turbine_velocities_one_direction!(points_x, points_y, ms2, pd2)
    in_speeds = pd2.wind_farm_states[1].turbine_inflow_velcities
    turbulent_damage = turbulence_damage(in_speeds[turb_index],turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180)
    # println("damage turbulence: ", turbulent_damage)
    damage[i] = waked_damage + turbulent_damage
end
subplot(132)
plot(sweep,damage)
ylim([0,1.5])
title("7")

sep = 10.0
# turbulence_intensity = 0.07863689518866869
# turbulence_intensity = 0.13827764626005523
turbulence_intensity = 0.05
# turbulence_intensity = 0.04
# turbulence_intensity = 0.046
# turbulence_intensity = 0.05933821834067302
ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms2 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
for i = 1:length(sweep)
    # println("sweep[i]: ", sweep[i])
    turbine_y = [0.0,sweep[i]*rotor_diameter]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, zeros(nturbines))
    pd2 = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])
    waked_damage = ff.get_damage(turbine_x,turbine_y,turb_index,speeds,omegas,pitches,wind_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,ms2,pd2,wind_speed,Nlocs=20,shearExp=shearexponent,mult=mult)
    # println("damage loads: ", damage[i])
    ff.turbine_velocities_one_direction!(points_x, points_y, ms2, pd2)
    in_speeds = pd2.wind_farm_states[1].turbine_inflow_velcities
    turbulent_damage = turbulence_damage(in_speeds[turb_index],turbulence_intensity,r,chord,theta,airfoils,Rhub,Rtip,hubHt,speeds,omegas,pitch;precone=2.5*pi/180,tilt=5*pi/180)
    # println("damage turbulence: ", turbulent_damage)
    damage[i] = waked_damage + turbulent_damage
end
subplot(133)
plot(sweep,damage)
ylim([0,1.5])
title("10")


FS = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]

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

    """high TI"""
    # FD4 = [0.59436966, 0.86531902, 1.47720059, 1.29034929, 0.39374063,
    #        0.0675834 , 0.09811606, 0.0944617 , 0.0546642 , 0.1133815 ,
    #        0.30375038]
    # FD7 = [0.74362544, 0.94617913, 1.12020256, 0.82845987, 0.41520515,
    #        0.09437026, 0.08533575, 0.09938972, 0.13268212, 0.24422962,
    #        0.28269637]
    # FD10 = [0.8467823 , 0.76667692, 0.7433349 , 0.69224899, 0.38166549,
    #        0.17593557, 0.12199661, 0.12614291, 0.12928943, 0.19476268,
    #        0.26522156]

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
end
subplot(131)
scatter(FS,FD4)

subplot(132)
scatter(FS,FD7)

subplot(133)
scatter(FS,FD10)


# println(damage)
show()
# # ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, ms2, pd2)
# # println(windfarmstate.turbine_inflow_velcities)
#
# # get_damage(turbine_x,turbine_y,turb_index,TSR,pitch,free_speed,Rhub,Rtip,r,chord,theta,airfoils,hubHt,windfarm,windfarmstate,windresource,wakedeficitmodel,wakedeflectionmodel,wakecombinationmodel)
