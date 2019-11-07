include("turbines.jl")
include("wake_models.jl")
using Plots

model = Jensen(0.1)

coord = Coord(0.0,0.0,0.0)
turbine = Turbine(coord,100.,90.,2.0/3.0)

xlocs = ones(50).*100.0
ylocs = range(-100,stop=100,length=50)
zlocs = ones(50).*turbine.hub_height
locs = zeros(length(xlocs),3)
locs[:,1] = xlocs[:]
locs[:,2] = ylocs[:]
locs[:,3] = zlocs[:]
loss_val = zeros(length(xlocs))

deflection = [0.0,0.0]
for i = 1:length(xlocs)
    loss_val[i] = wake_model(locs[i,:], deflection, model, turbine)
end

println(loss_val)
plot(ylocs,loss_val)
