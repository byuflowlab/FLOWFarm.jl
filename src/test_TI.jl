
include("model_set_3.jl")

x = rotor_diameter*4.0
x_vec = range(turbine_x[1],stop = turbine_x[1]+10.0*rotor_diameter,length=100)
y_vec = range(-4.0*rotor_diameter,stop = 4.0*rotor_diameter,length=100)
z_vec = range(0.0,stop = 2.0*rotor_diameter,length=100)
z = hub_height + rotor_diameter/2.0
println(z)
TI = zeros(length(x_vec))

for i = 1:length(TI)
           # loc = [x,y_vec[i],z]
           loc = [x_vec[i],0.0,z]
           # loc = [2.0*rotor_diameter,0.0,z_vec[i]]
           TI[i] = FlowFarm.GaussianTI(loc,windfarm,windfarmstate,ambient_ti)
end

figure(1)
plot((x_vec.-turbine_x[1])./rotor_diameter,TI.-ambient_ti)
# plot((y_vec.-turbine_x[1])./rotor_diameter,TI.-ambient_ti)
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
