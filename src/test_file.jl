# include("fatigue_model.jl")

# include("FAST_data.jl")
include("model.jl")
# include("model_set_2.jl")

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


turb = "low"
ws = 11.0
TI_free = 0.046

sep = 4.0

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

println("this file")
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
off = range(-1.0,stop=1.0,length=12)
# off = [-0.6]
dams = zeros(length(off))


# turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
# ky = ka*turbulence_intensity + kb
# kz = ka*turbulence_intensity + kb

# println(ky)
turbulence_intensity = 0.3
ky = 0.06
kz = 0.06
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel,local_ti_model)

turbine_x = [0.0,sep*rotor_diameter]
init_inflow_velcities = zeros(length(turbine_x)).+ws

omega_func = Akima(speeds, omegas)
pitch_func = Akima(speeds, pitches)

tilt = deg2rad(5.0)
rho = 1.225

diff_vel = 0.0
nCycles = 1000
# naz = 4
# az_arr = range(0.0,stop=2.0*pi-2.0*pi/naz,length=naz)
naz = 2
az_arr = [pi/2,3*pi/2]

turb_samples = randn(naz*nCycles)
start = time()

global t1_damage
global t2_damage

t1_damage = zeros(length(off))
t2_damage = zeros(length(off))

for k=1:length(off)

    global t1_damage
    global t2_damage

    println("offset: ", off[k])
    offset = off[k]*rotor_diameter
    turbine_y = [0.0,offset]
    windspeeds = [ws]
    windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
    windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, init_inflow_velcities, zeros(nturbines),zeros(nturbines).+ambient_ti)
    windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_ti, [wind_shear_model])
    pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])



    # dams[k] = get_single_damage(ms,pd,2,1,nCycles,az_arr,turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,
    #         r,rotor,sections,Rhub,Rtip)

    # # state_damage = get_single_state_damage(ms,pd,1,nCycles,az_arr,
    # #     turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,r,rotor,sections,Rhub,Rtip)
    t1 = time()
    total_damage = ff.get_total_farm_damage(ms,pd,nCycles,az_arr,
        turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,r,rotor,sections,Rhub,Rtip,precone,tilt,rho)
    t1_damage[k] = total_damage[1]
    t2_damage[k] = total_damage[2]
    println(time()-t1)
    # println(t2_damage)
end
# println((time()-start)/length(off))
# println("damage: ", dams)



scatter(off,t1_damage)
scatter(off,t2_damage)
# #
# FS,FD = fastdata(turb,ws,sep)
# scatter(FS,FD)
#
xlabel("offset")
ylabel("damage")
show()
