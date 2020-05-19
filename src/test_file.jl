# include("fatigue_model.jl")

# include("FAST_data.jl")
include("model.jl")
# include("model_set_2.jl")

function turbulence_function(loc,windfarm,windfarmstate,ambient_ti)
    r = sqrt(loc[2]^2+(loc[3]-90.0)^2)
    if loc[1] > 10.0
        if r < 70.0
            return 0.17
        else
            return 0.046
        end
    else
        return 0.046
    end
end


turb = "high"
ws = 11.0

sep = 7.0

points_x = [0.69,0,-0.69,0]
points_y = [0,0.69,0,-0.69]

fos = 1.15

rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
sections = CCBlade.Section.(r,chord,theta,airfoils)

# off = [-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
# off = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0]
off = range(-1.0,stop=1.0,length=25)
# off = range(-0.5,stop=-0.2,length=3)
# off = [-0.6]
dams = zeros(length(off))


# println(ky)
# ky = 0.02
# kz = 0.02
# alpha_star = 5.5
# beta_star = 0.5

#low turb
# ky = 0.015
# kz = 0.015
# alpha_star = 5.0
# beta_star = 0.5

#high turb
ky = 0.035
kz = 0.035
alpha_star = 6.5
beta_star = 0.8

# alpha_star = 1.0
# beta_star = 0.3
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star,1.0)
# wakedeficitmodel = ff.GaussYawVariableSpread(alpha_star,beta_star,1.0)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel,local_ti_model)

turbine_x = [0.0,sep*rotor_diameter]
init_inflow_velcities = zeros(length(turbine_x)).+ws

omega_func = Akima(speeds, omegas)
pitch_func = Akima(speeds, pitches)

tilt = deg2rad(5.0)
rho = 1.225

nCycles = 10000
# naz = 2
# az_arr = range(0.0,stop=2.0*pi-2.0*pi/naz,length=naz)
naz = 2
Nlocs = 200
az_arr = [pi/2,3*pi/2]
# az_arr = [0.0,pi/2,pi,3*pi/2]

# turb_samples = zeros(naz*nCycles)
# using Distributions
# alph = 2.0
# thet = 1.0
# weib = Weibull(alph,thet)
# ray = Rayleigh(alph)
# for i = 1:length(turb_samples)
#     turb_samples[i] = rand(weib)
# end
# using Statistics
# turb_samples = turb_samples./std(turb_samples)
# turb_samples = turb_samples .- mean(turb_samples)

#
turb_samples = randn(naz*nCycles)
# turb_samples = turb_samples + (rand(naz*nCycles).^2.0.-1.0)
# println(std(turb_samples))
# turb_samples .-= minimum(turb_samples)
# #0.92 good
# turb_samples = turb_samples.^0.7
# using Statistics
# turb_samples = turb_samples .- mean(turb_samples)
# turb_samples = turb_samples./std(turb_samples)
# turb_samples = turb_samples .- mean(turb_samples)
# # hist(turb_samples,bins=50)
# println(std(turb_samples))

# turb_samples = rand(naz*nCycles).*2.0 .- 1.0
#
#
run = true
if run == false
    # hist(turb_samples,bins=20)
    plot(velpoints,ctpoints)
elseif run == true
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
        windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, init_inflow_velcities, zeros(nturbines),zeros(nturbines).+ambient_ti,sorted_turbine_index)
        windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_ti, [wind_shear_model])
        pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])



        # dams[k] = get_single_damage(ms,pd,2,1,nCycles,az_arr,turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,
        #         r,rotor,sections,Rhub,Rtip)

        # # state_damage = get_single_state_damage(ms,pd,1,nCycles,az_arr,
        # #     turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,r,rotor,sections,Rhub,Rtip)
        t1 = time()
        total_damage = ff.get_total_farm_damage(ms,pd,nCycles,az_arr,
            turb_samples,points_x,points_y,omega_func,pitch_func,ff.GaussianTI,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,Nlocs=Nlocs)
        # total_damage = ff.get_total_farm_damage(ms,pd,nCycles,az_arr,
        #     turb_samples,points_x,points_y,omega_func,pitch_func,ff.GaussianTI_stanley,r,rotor,sections,Rhub,Rtip,precone,tilt,rho)
        # total_damage = ff.get_total_farm_damage(ms,pd,nCycles,az_arr,
        #     turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_function,r,rotor,sections,Rhub,Rtip,precone,tilt,rho)


        t1_damage[k] = total_damage[1]
        t2_damage[k] = total_damage[2]
        # println(time()-t1)
        # println(t2_damage)
    end
    # println((time()-start)/length(off))
    # println("damage: ", dams)



    scatter(off,t1_damage,color="C0")
    scatter(off,t2_damage,color="C3")

     include("FAST_data.jl")
    FS,FD = fastdata(turb,ws,sep)
    scatter(FS,FD,color="C2")
    #
    xlabel("offset (D)")
    ylabel("damage")


    legend(["front turbine","back turbine", "FAST"])
    title(sep)
    # ylim(-0.1,2.0)
    ylim(-0.1,1.1)
end
