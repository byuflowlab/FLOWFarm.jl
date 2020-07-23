#=== Set everything up for FlowFarm ===#
cd("/Users/nbaker/Documents/GitHub/FlowFarm.jl/test/")
include("iea37_specific_functions.jl")
using FlowFarm; const ff = FlowFarm
import YAML
using PyPlot
using CSV

#--- Read in windfarm boundary data ---#
# Which case study we're doing. 'cs3' or 'cs4'
str_case = "4"
#- Rip the boundary coordinates from the .yaml file -#
file_dir = "./inputfiles/"
bnry_file_name_orig = "iea37-boundary-cs" * str_case * ".yaml"
bnry_file_name = string(file_dir,bnry_file_name_orig)
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
bndry_x_clsd, bndry_y_clsd = ff.closeBndryLists(bndry_x, bndry_y)

# Place the turbines for the first area
layout_file_name = "./inputfiles/iea37-ex-opt4.yaml"
~, ~, fname_turb, ~ = ff.get_turb_loc_YAML(layout_file_name)
turbine_file_name = string("./inputfiles/",fname_turb)
~, ~, ~, ~, turb_diam, ~ = ff.get_turb_atrbt_YAML(turbine_file_name)
region = 1

# Plot boundary
bndry_x = bndry_x_clsd[region]
bndry_y = bndry_y_clsd[region]

plot(bndry_x, bndry_y)
num_bndry_pts = length(bndry_x)
println(num_bndry_pts)
#scatter(bndry_x[20],bndry_y[20])#, turb_diam/2.0, fill=true,color="black")
# Plot boundary verticies
for i = 1:num_bndry_pts
     scatter(bndry_x[i],bndry_y[i])#, turb_diam/2.0, fill=true,color="black")
     plt.text(bndry_x[i]+turb_diam, bndry_y[i]+turb_diam, string(i))
end

# Formatting
axis("square")
axis("off")
plt.show()