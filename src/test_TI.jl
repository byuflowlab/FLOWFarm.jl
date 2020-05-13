
include("model2.jl")

sep = 10.0
x = rotor_diameter*sep
x_vec = range(turbine_x[1],stop = turbine_x[1]+10.0*rotor_diameter,length=100)
y_vec = range(-2.0*rotor_diameter,stop = 2.0*rotor_diameter,length=100)
z_vec = range(0.0,stop = 2.0*rotor_diameter,length=100)
z = hub_height
println(z)
TI = zeros(length(y_vec))

for i = 1:length(TI)
           loc = [x,y_vec[i],z]
           # loc = [x_vec[i],0.0,z]
           # loc = [4.0*rotor_diameter,0.0,z_vec[i]]
           TI[i] = FlowFarm.GaussianTI_stanley(loc,windfarm,windfarmstate,ambient_ti[1])
end

figure(1)
# plot((x_vec.-turbine_x[1])./rotor_diameter,TI.-ambient_ti)
plot(y_vec./rotor_diameter,TI)
include("TIdata.jl")
dat = get_TIdat(11.0,"low",sep)
println(dat)
x_TI = range(-200.0/rotor_diameter,stop=200.0/rotor_diameter,length=length(dat))
scatter(x_TI,dat)




# plot((z_vec.-turbine_x[1])./rotor_diameter,TI.-ambient_ti)

# figure(2)
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter/2.0, fill=false,color="C1"))
# end
#
# plot([-1000.0,1000.0],[0.0,0.0])
# # plot([0.0,0.0],[-1000.0,1000.0])
# axis("square")
# xlim(minimum(turbine_x)-200,maximum(turbine_x)+200)
# ylim(minimum(turbine_y)-200,maximum(turbine_y)+200)

# axis("square")
# xlim(-boundary_radius-200,boundary_radius+200)
# ylim(-boundary_radius-200,boundary_radius+200)
