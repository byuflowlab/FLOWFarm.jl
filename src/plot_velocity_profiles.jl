using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics
using NPZ

const ff=FlowFarm

function calc_TI(constant,ai,TI_free,initial,sep,downstream)

    ti_calculation = constant * (1.0/3.0)^ai * TI_free^initial * sep^downstream

    TI = sqrt(ti_calculation^2 + TI_free^2)

    return TI
end

turb = "low"
ws = 12.0

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

include("model.jl")

sep = 4.0

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

L = 100
sweep = range(-1.5,stop=1.5,length=L)
U_SOWFA = interp2d(akima,x,y,u,sweep.*rotor_diameter,zeros(L).+hubHt)[:,1]
plot(sweep,U_SOWFA)
# sweep = range(0.0,stop=300.0,length=L)
# U_SOWFA = interp2d(akima,x,y,u,zeros(L),sweep)[1,:]
# plot(U_SOWFA,sweep)




TI_free = 0.046
TI = "low"
windspeeds = [ws]
turbine_inflow_velcities = zeros(nturbines) .+ ws
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])


ka = 0.38
kb = 0.004

initial = 0.313
constant = 1.931
ai = 0.435
downstream = -0.855
alpha_star = 2.32
beta_star = 0.154
turbulence_intensity = calc_TI(constant,ai,TI_free,initial,sep,downstream)
ky = ka*turbulence_intensity + kb
kz = ka*turbulence_intensity + kb
horizontal_spread_rate = ky
vertical_spread_rate = kz
wakedeficitmodel = ff.GaussYaw(turbulence_intensity,horizontal_spread_rate,vertical_spread_rate,alpha_star,beta_star)
wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
ms = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

turbine_x = [0.0,sep*rotor_diameter]
turbine_y = [0.0,0.0]
init_inflow_velcities = zeros(nturbines).+ws
windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definition_ids, turbine_definitions)
windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, init_inflow_velcities, zeros(nturbines))
pd = ff.WindFarmProblemDescription(windfarm, windresource, [windfarmstate])

U_model = zeros(L)
for i = 1:L
        offset = sweep[i]*rotor_diameter
        loc = [sep*rotor_diameter, sweep[i]*rotor_diameter, hubHt]
        # loc = [sep*rotor_diameter, 0.0, sweep[i]]
        U_model[i] = ff.point_velocity(loc, ms, pd)
end
plot(sweep,U_model)
# plot(U_model,sweep)
show()
