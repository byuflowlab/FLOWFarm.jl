include("turbines.jl")
include("wake_models.jl")
using Plots

model = Jensen(0.1)

coord = Coord(0.0,0.0,0.0)
turbine = Turbine(coord,100.,90.,1.0/3.0,0.0)
println(turbine)

num = 100
xlocs = ones(num).*800.0
ylocs = range(-100,stop=100,length=num)
zlocs = ones(num).*turbine.hub_height
locs = zeros(length(xlocs),3)
locs[:,1] = xlocs[:]
locs[:,2] = ylocs[:]
locs[:,3] = zlocs[:]
loss_val = zeros(num)

deflection = [0.0,0.0]
for i = 1:num
    loss_val[i] = wake_model(locs[i,:], deflection, model, turbine)
end

println(loss_val)
# plot(ylocs,loss_val)



model = Multizone([-0.5,0.3,1.0],0.05,[0.5,1.0,5.5],12.0,1.3)
for i = 1:num
    loss_val[i] = wake_model(locs[i,:], deflection, model, turbine)
end
println(loss_val)
plot(ylocs,loss_val)
